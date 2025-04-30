from pathlib import Path
from urllib.parse import urlparse
from pathlib import PurePosixPath
from subprocess import run
import logging
import pandas as pd
from bs4 import BeautifulSoup

from src.config.cnfg_database import TAXONKIT_DB, ASSEMBLY_SUMMARY_GENBANK, ASSEMBLY_SUMMARY_REFSEQ
from src.config.cnfg_software import MAMBA
from src.config.cnfg_taxonomy import BELOW_FAMILY_RANKS
from src.utils.util_command import execute_command
from src.utils.util_write import list2txt


def get_taxonomy_id_from_sciname(sciname: str, infodir: Path) -> list:
    """
    使用 taxonkit 将 scientific name 转为 taxonomy id.
    仅允许输入科以下的分类级别。

    :param sciname: scientific name
    :return: scientific name 在当前分类级别下的 taxonomy id 列表
    """
    logging.info(f"获取 {sciname} 的 taxonomy id 列表")
    # 获取 scientific name 的 taxonomy id 和 rank
    cmd_name2taxid = (
        f"echo '{sciname}' | {MAMBA} -n basic run taxonkit name2taxid --data-dir {TAXONKIT_DB} --show-rank --sci-name")
    output = execute_command(cmd_name2taxid)
    try:
        _, txid, rank = output.split("\t")
    except ValueError:
        raise ValueError(f"输出格式错误或为 homotypic synonym: {cmd_name2taxid}\n输出: {output}")
    # 检查 rank 是否在允许范围内, 科及一下级别
    if rank not in BELOW_FAMILY_RANKS:
        raise ValueError(f"当前 scientific name 的 rank 不在允许范围内(科及以下级别): {rank}")
    # 获取当前分类级别下的 taxonomy id 列表
    cmd_list_taxids = (
        f"{MAMBA} -n basic run taxonkit list --data-dir {TAXONKIT_DB} --indent '' --ids {txid}")
    taxid_output = execute_command(cmd_list_taxids)
    taxids = taxid_output.split("\n")
    if not taxids or taxids == [""]:
        raise RuntimeError(
            f"当前 scientific name 的 taxid 列表为空: {cmd_list_taxids}\n输出: {taxid_output}")
    # 结果写入到文件
    list2txt(taxids, infodir.joinpath("taxids.txt"))
    return taxids


def get_assembly_summary_by_taxids(taxids: list, infodir: Path) -> pd.DataFrame:
    """
    筛选 assembly summary 中 taxids 的条目, 返回 dataframe
    :param taxids: taxonomy id 列表
    :param infodir: 输出目录
    :return: 当前 taxids 的 assembly summary dataframe
    """
    logging.info(f"获取 {taxids} 的 refseq&genbank assembly summary")
    # 读入 RefSeq 的 assembly_summary_refseq.txt 文件
    usecols = ["#assembly_accession", "refseq_category", "taxid", "species_taxid", "organism_name",
               "assembly_level", "gbrs_paired_asm", "ftp_path", "genome_size", "scaffold_count", "contig_count"]
    df_rf = pd.read_csv(ASSEMBLY_SUMMARY_REFSEQ, sep="\t", skiprows=1,
                        na_values=["na", ""], dtype=object, usecols=usecols)
    df_gb = pd.read_csv(ASSEMBLY_SUMMARY_GENBANK, sep="\t", skiprows=1,
                        na_values=["na", ""], dtype=object, usecols=usecols)
    df_rs_cur_tax = df_rf.loc[df_rf["taxid"].isin(taxids),]
    # refseq 中包含的条目不要重复下载
    gbrs_paired_asms = df_rs_cur_tax["gbrs_paired_asm"].dropna().unique()
    df_gb_cur_tax = df_gb[
        (df_gb["taxid"].isin(taxids)) & (~df_gb["#assembly_accession"].isin(gbrs_paired_asms))]
    df_rsgb_cnct = pd.concat([df_rs_cur_tax, df_gb_cur_tax], axis=0, ignore_index=True)
    # 输出到文件
    df_rsgb_cnct.to_csv(Path(infodir).joinpath(
        "assembly_summary_concat_refseq_genbank.tsv"), sep="\t", index=False)
    return df_rsgb_cnct


def download_genome_files(df_asmb_smry: pd.DataFrame, alldir: Path) -> None:
    """
    下载 genome 目录下面的指定文件, fna/gff/gtf/faa
    :param df_rsgb_cnct: 当前物种 assembly summary dataframe
    :param alldir: 当前物种存放 genome 的目录
    """
    logging.info(f"下载 {df_asmb_smry.shape[0]} 个基因组的文件")
    # 迭代每行, 每行为一个基因组
    for row in df_asmb_smry.iterrows():
        asmb_acc = row[1]["#assembly_accession"]
        ftp_path = row[1]["ftp_path"]
        prfx = PurePosixPath(urlparse(ftp_path).path).name
        # 创建当前基因组的目录
        dir_cur_gnm = Path(alldir).joinpath(asmb_acc)
        dir_cur_gnm.mkdir(parents=True, exist_ok=True)
        # 获取当前基因组的 ftp 目录页面, 要把 '/' 加上
        res = run(f"curl -l {ftp_path}/", shell=True, capture_output=True, text=True, check=True)
        html_content = res.stdout.strip()
        # 查看 fna, gtf, gff, faa 哪些文件可以下载
        target_files = [prfx + kw for kw in ["_genomic.fna.gz",
                                             "_genomic.gff.gz", "_genomic.gtf.gz", "_protein.faa.gz"]]
        existed_links = []
        soup = BeautifulSoup(html_content, "html.parser")
        for a_tag in soup.find_all("a", href=True):
            href = a_tag["href"]
            if href in target_files:
                existed_links.append(href)
        # 开始下载, 写入到文件, 方便后续没成功的文件继续下载
        with open(f"{dir_cur_gnm}/dwnld.sh", "w") as f:
            f.write(f"wget -c {ftp_path}/md5checksums.txt -O {dir_cur_gnm}/md5checksums.txt\n")
            for link in existed_links:
                f.write(
                    f"wget -c {ftp_path}/{link} -O {dir_cur_gnm}/{PurePosixPath(urlparse(link).path).name}\n")
        run(f"bash {dir_cur_gnm}/dwnld.sh", shell=True, check=True, timeout=600)
        # MD5 校验
        res = run(f"cd {dir_cur_gnm} && md5sum -c md5checksums.txt | grep -c OK",
                  shell=True, check=True, capture_output=True, text=True)
        if int(res.stdout.strip()) == len(existed_links):
            run(f"touch {dir_cur_gnm}/md5checksums.OK", shell=True, check=True)
        else:
            run(f"touch {dir_cur_gnm}/md5checksums.FAILED", shell=True, check=True)

import logging
import subprocess
import sys
from pathlib import Path
from urllib.parse import urlparse
import click
import pandas as pd

from src.config.software import ASCP, MAMBA
from src.config.configure import PRIVATE_KEY, TAXONKIT_DB, PROJECT_DIR

logger = logging.getLogger(__file__)
logger.addHandler(logging.NullHandler())

logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s - %(levelname)s - %(message)s')


def df_dup_remove(refseq, genbank):
    """
    去除重复的基因组，保留refseq中的基因组
    :param refseq: refseq summary 数据框
    :param genbank: genbank summary 数据库
    :return: 返回去重后的数据框
    """
    genbank = genbank[~genbank["gbrs_paired_asm"].isin(refseq["assembly_accession"])]
    res = pd.concat([refseq, genbank]).reset_index(drop=True)
    return res


def query_summary(dataframe, taxonid):
    """
    通过taxonid查询summary文件，获取下载路径
    :param dataframe: summary 数据框
    :param taxonid: taxonid
    :return: 返回包含taxonid的summary数据框
    """
    res = []
    index = set()
    for i in taxonid:
        index.update(dataframe[dataframe["taxid"] == int(i)].index)
        index.update(dataframe[dataframe["species_taxid"] == int(i)].index)
    if len(index) > 0:
        return dataframe.loc[list(index)]
    else:
        return None


def ascp_download(url, d_out: Path):
    """
    使用ascp下载文件
    :param url: ftp路径
    :param d_out: 输出路径
    :return: 返回下载的文件路径
    """
    ori_path = urlparse(url)
    stem_path = ori_path.path
    name = stem_path.strip().split('/')[-1]
    cmd = f"{ASCP} -i {PRIVATE_KEY} -k 1 -T -l100m anonftp@ftp.ncbi.nlm.nih.gov:{url} {d_out}"
    logging.debug(f"运行命令: {cmd}")
    ret = subprocess.run(cmd, shell=True)
    if ret.returncode == 0:
        return d_out.joinpath(name)
    else:
        sys.exit(logger.error(f"{ret.returncode}: {cmd}"))


def expand_taxon(taxon):
    """
    将taxonid展开，获取其子id
    :param taxon: taxonid
    :return: 返回展开后的taxonid集合
    """
    res = set()
    cmd = f"{MAMBA} run -n basic taxonkit list --data-dir {TAXONKIT_DB} -i {taxon}"
    ret = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE)
    for i in ret.stdout.decode().strip().split("\n"):
        if i.strip():
            res.add(i.strip())
    return res


def get_download_cmd(ori_path, d_out):
    """
    从 ftp 获取下载路径
    :param ori_path: ftp路径
    :param d_out: 输出路径
    :return: 返回下载命令
    """
    ori_path = urlparse(ori_path)
    stem_path = ori_path.path
    name = stem_path.strip().split('/')[-1]
    cmd = f"{ASCP} -i {PRIVATE_KEY} -k 1 -T -l100m anonftp@ftp.ncbi.nlm.nih.gov:{stem_path}/{name}_genomic.fna.gz {d_out}"
    return cmd


def file_exist(cmd):
    """
    过滤已下载完成的基因组的下载命令
    :param cmd: 下载命令
    :return: 返回是否存在下载完成的基因组
    """
    download_path, d_out = cmd.strip().split(':')[-1].strip().split()
    f_name = download_path.strip().split('/')[-1]
    f_out = Path(f"{d_out}/{f_name}")
    if f_out.exists():
        return True
    else:
        return False


# Main
@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('--taxonid',
              required=True,
              type=click.Path(),
              help="taxonid, 使用逗号分隔")
@click.option('--ptype',
              required=True,
              type=click.Choice(['Bacteria', "Viruses", "Fungi"]),
              help="物种类型")
@click.option('--db',
              required=True,
              type=click.Path(),
              help="数据库目录")
@click.option('--name',
              required=True,
              type=click.Path(),
              help="物种数据库目录名称")
@click.option('--latin',
              type=click.Path(),
              required=True,
              help="拉丁名")
@click.option("--summary_refseq",
              type=click.Path(),
              required=True,
              help="assembly_summary_refseq.txt 文件, 下载地址 https://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt")
@click.option("--summary_genbank",
              type=click.Path(),
              required=True,
              help="assembly_summary_genbank.txt 文件, 下载地址 https://ftp.ncbi.nlm.nih.gov/genomes/genbank/assembly_summary_genbank.txt")
@click.option("--thread",
              default=16,
              show_default=True,
              type=int,
              help="下载线程数")
@click.option("--timeout",
              default=45,
              type=int,
              show_default=True,
              help="超时等待时间")
@click.option("--run/--no-run",
              default=False,
              show_default=True,
              help="是否运行下载脚本")
@click.option("--force/--no-force",
              default=False,
              show_default=True,
              help="是否强制下载")
@click.option("--filter/--no-filter",
              default=True,
              show_default=True,
              help="是否通过 excluded_from_refseq 筛选 summary 基因组")
def main(taxonid, ptype, db, name, latin, summary_refseq, summary_genbank, thread, timeout, run, force, filter):
    """
    从 NCBI 的 refseq 和 genbank 数据库下载基因组数据
    """
    d_database = Path(db).absolute()
    d_database.mkdir(exist_ok=True, parents=True)

    # 数据下载
    d_out = d_database.joinpath(name)
    d_out.mkdir(exist_ok=True, parents=True)

    l_taxonid = set()
    logger.info(f"Expand TaxonID {taxonid}")
    for i in taxonid.strip().split(','):
        l_taxonid.update(expand_taxon(i))

    filter_keyword = ["contaminated",
                      "chimeric",
                      "low quality sequence",
                      "misassembled",
                      "mixed culture",
                      "derived from environmental source",
                      "derived from metagenome",
                      "metagenome"]
    filter_keyword = '|'.join(filter_keyword)

    logger.info(f"Read the summary file: {summary_refseq}")
    df_refseq = pd.read_csv(summary_refseq, sep="\t", skiprows=1,
                            na_values=['na', ''], low_memory=False)
    df_refseq.rename(columns={df_refseq.columns[0]: "assembly_accession"}, inplace=True)
    df_refseq = df_refseq[df_refseq["ftp_path"].notna()]
    if filter:
        df_refseq = df_refseq[~df_refseq["excluded_from_refseq"].str.contains(
            filter_keyword, na=False)]

    logger.info(f"Read the summary file: {summary_genbank}")
    df_genbank = pd.read_csv(summary_genbank, sep="\t", skiprows=1,
                             na_values=['na', ''], low_memory=False)
    df_genbank.rename(columns={df_genbank.columns[0]: "assembly_accession"}, inplace=True)
    df_genbank = df_genbank[df_genbank["ftp_path"].notna()]
    if filter:
        df_genbank = df_genbank[~df_genbank["excluded_from_refseq"].str.contains(
            filter_keyword, na=False)]

    logger.info(f"Remove duplicate assem")
    df_summary = df_dup_remove(df_refseq, df_genbank)

    logger.info(f"Generate download script")
    d_seq = d_out.joinpath("seq")
    d_seq.mkdir(exist_ok=True, parents=True)
    f_download = d_seq.joinpath("download.sh")
    f_metadata = d_out.joinpath("metadata.tsv")
    df_query = query_summary(df_summary, l_taxonid)
    if df_query is not None:
        df_query.to_csv(f_metadata, sep="\t", index=False)
    try:
        cmd = df_query["ftp_path"].apply(get_download_cmd, d_out=d_seq)
    except:
        sys.exit(logger.error(df_query.loc[:, ["asm_name", "ftp_path"]]))
    with open(f_download, 'w') as OUT:
        for i in cmd:
            if force:
                print(i, file=OUT)
            else:
                if file_exist(i):
                    pass
                else:
                    print(i, file=OUT)

    f_script = d_out.joinpath("download.sh")
    with open(f_script, 'w') as OUT:
        print(f"{MAMBA} run -n basic parallel --timeout {timeout} -j {thread} < {f_download}", file=OUT)
        if ptype == "Viruses":
            print(f"""
{MAMBA} run -n python3.12 poetry -C {PROJECT_DIR} run \
    python -m src.kml_qpcr.entrez_search \
    '({latin})' '(complete genome)' -o {d_seq}/list.txt""",
                file=OUT)
            print(f"""
{MAMBA} run -n python3.12 poetry -C {PROJECT_DIR} run \
    python src.kml_qpcr.entrez_download -t {d_seq}/list.txt -p {d_seq}/nucleotide""",
                file=OUT)

    if run:
        logger.info(f"Start to download...")
        subprocess.run(f"bash {f_script} > {f_script}.log 2>&1", shell=True)
        logger.info(f"Finish download")
    else:
        logger.info(f"Run the script {f_script} manual")


if __name__ == "__main__":
    main()

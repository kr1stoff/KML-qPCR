from pathlib import Path
from subprocess import run
import logging
from functools import reduce
import pandas as pd

from src.config.cnfg_software import CSVTK, PARALLEL, ACTIVATE
from src.config.cnfg_database import CHECKV_DB


def assess_genome_quality(sci_name: str, genome_set_dir: str, threads: int, pthgn_type: str, force: bool = False):
    """
    基因组质量评估
    :gnmdir: 基因组目录, 例 /data/mengxf/Project/KML250416_chinacdc_pcr/genomes/Ehrlichia_chaffeensis
    :threads: 线程数
    :pthgn_type: 病原类型, 例 Bacteria 或 Viruses
    :force: 是否强制重新运行 checkM/checkV, 默认识别到结果文件就跳过
    :return: None

    目录结构:
    genomes/Ehrlichia_chaffeensis
    ├── all
    │   ├── GCF_000013145.1
    │   ├── GCF_000632815.1
    │   ├── ...
    │   └── GCF_000632965.10
    └── ...
    """
    logging.info(f"开始基因组质量评估: {sci_name}, 病原类型: {pthgn_type}, 线程数: {threads}, 强制重新运行: {force}")
    gnmdir = Path(genome_set_dir).joinpath(sci_name.replace(" ", "_"))
    if pthgn_type == "Bacteria":
        # 细菌
        run_checkm(gnmdir, threads, force)
        # filter_by_checkm(gnmdir)
    else:
        # 病毒
        run_checkv(gnmdir, threads, force)
        filter_by_checkv(gnmdir)


def run_checkm(gnmdir: Path, threads: int, force: bool):
    """
    运行 checkM
    checkM 评估 bins 多个基因组
    gnmdir: 基因组目录, 例 /data/mengxf/Project/KML250416_chinacdc_pcr/genomes/Ehrlichia_chaffeensis
    threads: 线程数
    """
    # 创建 checkM 输入 bins 目录，复制并解压所有基因组
    all_dir = gnmdir.joinpath("all")
    bins_dir = gnmdir.joinpath("genome_assess/checkm_input")
    bins_dir.mkdir(parents=True, exist_ok=True)
    run(f"cp {all_dir}/*/*_genomic.fna {bins_dir}", shell=True, check=True)
    run(f"gunzip --force {bins_dir}/*gz", shell=True, check=True)
    # 运行 checkM
    checkm_dir = gnmdir.joinpath("genome_assess/checkm")
    checkm_dir.mkdir(parents=True, exist_ok=True)
    result_file = checkm_dir.joinpath("result.tsv")
    # * 如果结果文件已存在且不强制运行，则跳过
    if result_file.exists() and not force:
        logging.info(f"checkM 结果文件 {result_file} 已存在, 跳过运行.")
        return
    # 并行太多 mamba run 会报文件锁的错误
    checkm_cmd = f"""
    source {ACTIVATE} qpcr
    checkm lineage_wf \
        -x fna --tab_table -t {threads} --pplacer_threads {threads} \
        -f {result_file} \
        {bins_dir} \
        {checkm_dir}
    conda deactivate
    """
    logging.debug(f"运行 checkM: {checkm_cmd}")
    run(checkm_cmd, shell=True, check=True, executable="/bin/bash")
    run(f"{CSVTK} -t csv2xlsx {result_file}", shell=True, check=True)


def get_genome_anno_quality(gnmdir: Path):
    """
    获取基因组注释质量
    :param gnmdir: 基因组目录
    """
    # 输入 Prokka 注释结果的目录
    prk_dir = gnmdir / "genome_annotate"
    # 获取去所有 prokka.tsv 特征文件中 rRNA, tRNA 数量
    rna_mtrx = []
    for feat_file in prk_dir.glob("*/prokka.tsv"):
        # /data/mengxf/Project/KML250416_chinacdc_pcr/genomes/Coxiella_Burnetii/genome_annotate/GCA_000300315.1/prokka.tsv
        gnm_num = feat_file.parent.name
        df = pd.read_csv(feat_file, sep="\t")
        # tRNA 数量
        trna_count = df[df["ftype"] == "tRNA"].shape[0]
        # rRNA 数量, 分为 16S, 23S, 5S
        rrna_df = df[df["ftype"] == "rRNA"]
        rrna_16s_count = rrna_df[rrna_df["product"].str.contains("16S")].shape[0]
        rrna_23s_count = rrna_df[rrna_df["product"].str.contains("23S")].shape[0]
        rrna_5s_count = rrna_df[rrna_df["product"].str.contains("5S")].shape[0]
        rna_mtrx.append([gnm_num, rrna_23s_count, rrna_16s_count, rrna_5s_count, trna_count])
    # 转换为 DataFrame
    rna_df = pd.DataFrame(rna_mtrx, columns=["genome_id", "23SrRNA", "16SrRNA", "5SrRNA", "tRNA"])
    # 保存为 CSV
    rna_df.to_csv(gnmdir / "genome_assess/rna_quality.csv", index=False)


def merge_checkm_rna_datafram(assess_dir: Path):
    """
    """


# TODO 加入了注释质量，需要重新过滤
# def filter_by_checkm(gnmdir: Path):
#     """
#     根据 checkM 结果过滤基因组, 生成 high_quality_genomes.txt 文件
#     :param checkm_dir: checkM 结果目录
#     :return: None
#     """
#     out_gnm_file = gnmdir / "genome_assess/high_quality_genomes.txt"
#     checkm_dir = gnmdir.joinpath("genome_assess/checkm")
#     restab = f"{checkm_dir}/result.tsv"
#     df = pd.read_csv(restab, sep="\t", usecols=["Bin Id", "Completeness", "Contamination"])
#     # ! [250526 FJH] 过滤条件
#     # 1.基因组完整度 (Completeness) ≥ 90%
#     # 2.污染度 (Contamination) ≤ 5%
#     fltr_df = df[(df["Completeness"] >= 90) & (df["Contamination"] <= 5)].copy()
#     fltr_df["Genome"] = fltr_df["Bin Id"].apply(lambda x: "_".join(x.split("_")[:2]))
#     fltr_df["Genome"].drop_duplicates().to_csv(out_gnm_file, index=False, header=False)


def run_checkv(gnmdir: Path, threads: int, force: bool):
    """
    运行 checkV
    checkV 评估单个基因组，不能像 checkM 那样评估 bins 多个基因组. 运行单个基因组的 checkV 后合并表格
    gnmdir: 基因组目录, 例 /data/mengxf/Project/KML250416_chinacdc_pcr/genomes/Bandavirus_dabieense
    threads: 线程数
    """
    # 配置输出目录
    all_dir = gnmdir.joinpath("all")
    checkv_dir = gnmdir.joinpath("genome_assess/checkv")
    checkv_bins_dir = checkv_dir / "bins"
    checkv_bins_dir.mkdir(parents=True, exist_ok=True)
    # * 如果结果文件已存在且不强制运行，则跳过
    result_file = checkv_dir.joinpath("checkv_summary.tsv")
    if result_file.exists() and not force:
        logging.info(f"checkV 结果文件 {result_file} 已存在, 跳过运行.")
        return
    # 并行和线程
    prl_num = min(4, threads)
    sgl_thrd = threads // prl_num
    # 搜索所有的 fna 文件, 写入批量运行脚本
    with open(checkv_dir / "checkv_batch.sh", "w") as f:
        for fna in all_dir.glob("*/*.fna"):
            checkv_cmd = f"checkv end_to_end -t {sgl_thrd} -d {CHECKV_DB} {fna} {checkv_bins_dir}/{fna.parent.name}\n"
            f.write(checkv_cmd)
    # 运行失败就中断
    checkv_cmd = f"""
    source {ACTIVATE} qpcr
    cat {checkv_dir}/checkv_batch.sh | {PARALLEL} -j {prl_num}
    conda deactivate
    """
    res = run(checkv_cmd, shell=True, check=True, executable="/bin/bash")
    logging.debug(f"运行 checkV: {checkv_cmd}")
    if res.returncode != 0:
        raise RuntimeError("CheckV 批量运行失败.")
    # 合并结果
    qlt_smrys = list(checkv_bins_dir.glob("*/quality_summary.tsv"))
    dfs = []
    for altsmr in qlt_smrys:
        df = pd.read_csv(altsmr, sep="\t")
        df.insert(0, "genome_id", altsmr.parent.name)
        dfs.append(df)
    dfmrg = reduce(lambda x, y: pd.concat([x, y]), dfs)
    dfmrg.to_excel(checkv_dir / "checkv_summary.xlsx", index=False)
    dfmrg.to_csv(result_file, sep="\t", index=False)


def filter_by_checkv(gnmdir: Path):
    """
    根据 checkV 结果过滤基因组, 生成 high_quality_genomes.txt 文件
    :param checkv_dir: checkV 结果目录
    :return: None
    """
    checkv_dir = gnmdir / "genome_assess/checkv"
    out_gnm_file = gnmdir / "genome_assess/high_quality_genomes.txt"
    restab = checkv_dir / "checkv_summary.tsv"
    # 病毒 checkv
    df = pd.read_csv(restab, sep="\t", usecols=[
                     "genome_id", "checkv_quality", "warnings", "completeness", "contamination"])
    # ! [250526 FJH] 过滤条件
    # 1.基因组质量(Checkv_Quality) 为高/中/低质量
    # 2.无warnings
    # 3.基因组完整度(Completeness) 为100. PS.如为基因组为分段形式则以总和计算.
    # 4.污染度(Contamination) ≤5
    # checkv_quality 1: high quality, 2: medium quality, 3: low quality
    to_move_genomes = df[~df['checkv_quality'].isin(
        ['Low-quality', 'Medium-quality', 'High-quality'])]['genome_id'].unique().tolist()
    df = df[~df['genome_id'].isin(to_move_genomes)]
    # 无 warnings 和 contamination <= 5
    df = df[(df['warnings'].isna()) & (df['contamination'] <= 5)][["genome_id", "completeness"]]
    # 基因组完整度 (Completeness) 为100 (四舍五入). PS.如为基因组为分段形式则以总和计算
    df = round(df.groupby('genome_id').sum())
    df[df['completeness'] == 100].index.to_series().to_csv(out_gnm_file, index=False, header=False)

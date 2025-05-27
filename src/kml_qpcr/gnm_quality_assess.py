from pathlib import Path
from subprocess import run
import logging
from functools import reduce
import pandas as pd

from src.config.cnfg_software import CSVTK, PARALLEL, ACTIVATE
from src.config.cnfg_database import CHECKV_DB


def assess_genome_quality(genome_dir: Path, threads: int, pthgn_type):
    """
    基因组质量评估
    genome_dir: 基因组目录, 例 /data/mengxf/Project/KML250416_chinacdc_pcr/genomes/Ehrlichia_chaffeensis
    threads: 线程数
    目录结构:
    genomes/Ehrlichia_chaffeensis
    ├── all
    │   ├── GCF_000013145.1
    │   ├── GCF_000632815.1
    │   ├── ...
    │   └── GCF_000632965.10
    └── ...
    """
    if pthgn_type == "Bacteria":
        # 细菌
        run_checkm(genome_dir, threads)
        filter_by_checkm(genome_dir)
    else:
        # 病毒
        run_checkv(genome_dir, threads)
        filter_by_checkv(genome_dir)


def run_checkm(genome_dir: Path, threads: int):
    """
    运行 checkM
    checkM 评估 bins 多个基因组
    genome_dir: 基因组目录, 例 /data/mengxf/Project/KML250416_chinacdc_pcr/genomes/Ehrlichia_chaffeensis
    threads: 线程数
    """
    # 创建 checkM 输入 bins 目录，复制并解压所有基因组
    all_dir = genome_dir.joinpath("all")
    bins_dir = all_dir.joinpath("bins")
    bins_dir.mkdir(parents=True, exist_ok=True)
    run(f"cp {all_dir}/*/*_genomic.fna.gz {bins_dir}", shell=True, check=True)
    run(f"gunzip --force {bins_dir}/*gz", shell=True, check=True)
    # 运行 checkM
    checkm_dir = genome_dir.joinpath("genome_assess/checkm")
    checkm_dir.mkdir(parents=True, exist_ok=True)
    # 并行太多 mamba run 会报文件锁的错误
    checkm_cmd = f"""
    source {ACTIVATE} qpcr
    checkm lineage_wf \
        -x fna --tab_table -t {threads} --pplacer_threads {threads} \
        -f {checkm_dir}/result.tsv \
        {bins_dir} \
        {checkm_dir}
    conda deactivate
    """
    logging.debug(f"运行 checkM: {checkm_cmd}")
    run(checkm_cmd, shell=True, check=True, executable="/bin/bash")
    run(f"{CSVTK} -t csv2xlsx {checkm_dir}/result.tsv", shell=True, check=True)


def filter_by_checkm(genome_dir):
    """
    根据 checkM 结果过滤基因组, 在 checkm 结果目录下生成 high_quality_genomes.txt 文件
    :param checkm_dir: checkM 结果目录
    :return: None
    """
    checkm_dir = genome_dir.joinpath("genome_assess/checkm")
    restab = f"{checkm_dir}/result.tsv"
    out_gnm_file = f"{checkm_dir}/high_quality_genomes.txt"
    df = pd.read_csv(restab, sep="\t", usecols=["Bin Id", "Completeness", "Contamination"])
    # ! [250526 FJH] 过滤条件
    # 1.基因组完整度 (Completeness) ≥ 90%
    # 2.污染度 (Contamination) ≤ 5%
    df[(df["Completeness"] >= 90) & (df["Contamination"] <= 5)
       ]["Bin Id"].to_csv(out_gnm_file, index=False, header=False)


def run_checkv(genome_dir: Path, threads: int):
    """
    运行 checkV
    checkV 评估单个基因组，不能像 checkM 那样评估 bins 多个基因组. 运行单个基因组的 checkV 后合并表格
    genome_dir: 基因组目录, 例 /data/mengxf/Project/KML250416_chinacdc_pcr/genomes/Bandavirus_dabieense
    threads: 线程数
    """
    # 配置输出目录
    all_dir = genome_dir.joinpath("all")
    checkv_dir = genome_dir.joinpath("genome_assess/checkv")
    checkv_bins_dir = checkv_dir / "bins"
    checkv_bins_dir.mkdir(parents=True, exist_ok=True)
    # 并行和线程
    prl_num = 4
    sgl_thrd = threads // prl_num
    # 搜索所有的 fna.gz 文件, 写入批量运行脚本
    with open(checkv_dir / "checkv_batch.sh", "w") as f:
        for fna in all_dir.glob("GC*/*_genomic.fna.gz"):
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
    qlt_smrys = list(checkv_bins_dir.glob("GC*/quality_summary.tsv"))
    dfs = []
    for altsmr in qlt_smrys:
        df = pd.read_csv(altsmr, sep="\t")
        df.insert(0, "genome_id", altsmr.parent.name)
        dfs.append(df)
    dfmrg = reduce(lambda x, y: pd.concat([x, y]), dfs)
    dfmrg.to_excel(checkv_dir / "checkv_summary.xlsx", index=False)
    dfmrg.to_csv(checkv_dir / "checkv_summary.tsv", sep="\t", index=False)


def filter_by_checkv(genome_dir):
    """
    根据 checkV 结果过滤基因组, 在 checkv 结果目录下生成 high_quality_genomes.txt 文件
    :param checkv_dir: checkV 结果目录
    :return: None
    """
    checkv_dir = genome_dir.joinpath("genome_assess/checkv")
    restab = f"{checkv_dir}/checkv_summary.tsv"
    out_gnm_file = f"{checkv_dir}/high_quality_genomes.txt"
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

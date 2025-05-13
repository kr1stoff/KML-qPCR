from pathlib import Path
from subprocess import run
import logging
from functools import reduce
import pandas as pd

from src.config.cnfg_software import MAMBA
from src.config.cnfg_database import CHECKV_DB


def assess_genome_quality(dir_genome: Path, threads: int, pthgn_type):
    """
    基因组质量评估
    dir_genome: 基因组目录, 例 /data/mengxf/Project/KML250416_chinacdc_pcr/genomes/Ehrlichia_chaffeensis
    threads: 线程数
    目录结构:
    genomes/Ehrlichia_chaffeensis
    ├── all
    │   ├── GCF_000013145.1
    │   ├── GCF_000632815.1
    │   ├── ...
    │   └── GCF_000632965.1
    └── ...
    """
    if pthgn_type == "Bacteria":
        # 细菌
        run_checkm(dir_genome, threads)
    else:
        # 病毒
        run_checkv(dir_genome, threads)

    # TODO 后入按需添加 BUSCO, QUAST, etc.


def run_checkm(dir_genome: Path, threads: int):
    """
    运行 checkM
    checkM 评估 bins 多个基因组
    dir_genome: 基因组目录, 例 /data/mengxf/Project/KML250416_chinacdc_pcr/genomes/Ehrlichia_chaffeensis
    threads: 线程数
    """
    # 创建 checkM 输入 bins 目录，复制并解压所有基因组
    dir_all = dir_genome.joinpath("all")
    dir_bins = dir_all.joinpath("bins")
    dir_bins.mkdir(parents=True, exist_ok=True)
    run(f"cp {dir_all}/*/*_genomic.fna.gz {dir_bins}", shell=True, check=True)
    run(f"gunzip --force {dir_bins}/*gz", shell=True, check=True)
    # 运行 checkM
    dir_checkm = dir_genome.joinpath("genome_assess/checkm")
    dir_checkm.mkdir(parents=True, exist_ok=True)
    cmd_checkm = f"""
    {MAMBA} -n qpcr run checkm lineage_wf \
        -x fna --tab_table -t {threads} --pplacer_threads {threads} \
        -f {dir_checkm}/result.tsv \
        {dir_bins} \
        {dir_checkm}
    """
    logging.debug(f"运行 checkM: {cmd_checkm}")
    run(cmd_checkm, shell=True, check=True)
    run(f"{MAMBA} -n basic run csvtk -t csv2xlsx {dir_checkm}/result.tsv", shell=True, check=True)


def run_checkv(dir_genome: Path, threads: int):
    """
    运行 checkV
    checkV 评估单个基因组，不能像 checkM 那样评估 bins 多个基因组。运行单个基因组的 checkV 后合并表格
    dir_genome: 基因组目录, 例 /data/mengxf/Project/KML250416_chinacdc_pcr/genomes/Bandavirus_dabieense
    threads: 线程数
    """
    # 配置输出目录
    dir_all = dir_genome.joinpath("all")
    dir_checkv = dir_genome.joinpath("genome_assess/checkv")
    dir_checkv_bins = dir_checkv / "bins"
    dir_checkv_bins.mkdir(parents=True, exist_ok=True)
    # 并行和线程
    prl_num = 4
    sgl_thrd = threads // prl_num
    # 搜索所有的 fna.gz 文件, 写入批量运行脚本
    with open(dir_checkv / "checkv_batch.sh", "w") as f:
        for fna in dir_all.glob("GC*/*_genomic.fna.gz"):
            cmd_checkv = f"{MAMBA} -n qpcr run checkv end_to_end -t {sgl_thrd} -d {CHECKV_DB} {fna} {dir_checkv_bins}/{fna.parent.name}\n"
            f.write(cmd_checkv)
    # 运行失败就中断
    res = run(
        f"cat {dir_checkv}/checkv_batch.sh | {MAMBA} -n basic run parallel -j {prl_num}", shell=True, check=True)
    logging.debug(f"运行 checkV: {res.args}")
    if res.returncode != 0:
        raise RuntimeError("CheckV 批量运行失败.")
    # 合并结果
    qlt_smrys = list(dir_checkv_bins.glob("GC*/quality_summary.tsv"))
    dfs = []
    for altsmr in qlt_smrys:
        df = pd.read_csv(altsmr, sep="\t")
        df.insert(0, "genome_id", altsmr.parent.name)
        dfs.append(df)
    dfmrg = reduce(lambda x, y: pd.concat([x, y]), dfs)
    dfmrg.to_excel(dir_checkv / "checkv_summary.xlsx", index=False)
    dfmrg.to_csv(dir_checkv / "checkv_summary.tsv", sep="\t", index=False)

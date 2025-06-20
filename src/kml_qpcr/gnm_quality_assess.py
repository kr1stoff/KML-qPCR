from pathlib import Path
from subprocess import run
import logging
from functools import reduce
import pandas as pd

from src.kml_qpcr.base import BaseQPCR
from src.config.cnfg_software import CSVTK, ACTIVATE
from src.config.cnfg_database import CHECKV_DB
from src.utils.util_command import multi_run_command


class GenomeQualityAssessor(BaseQPCR):
    def __init__(self, sci_name: str, genome_set_dir: str, threads: int, force: bool = False):
        """
        初始化基因组质量评估器
        :sci_name: 物种名称, 例 Ehrlichia_chaffeensis
        :gnm_dgenome_set_dirir: 基因组目录, 例 KML250416_chinacdc_pcr/genomes/Ehrlichia_chaffeensis
        :threads: 线程数
        :force: 是否强制重新运行 checkM/checkV, 默认识别到结果文件就跳过
        """
        super().__init__(sci_name, genome_set_dir, threads, force)
        # 分离株基因组评估目录
        self.assess_dir = self.gnm_dir / "genome_assess"

    def run(self):
        """基因组质量评估流程"""
        logging.info(f"开始基因组质量评估: {self.gnm_dir}, 线程数: {self.threads}, 强制重新运行: {self.force}")
        self.run_checkm()
        self.get_genome_anno_quality()
        self.merge_checkm_rna_stats()
        self.filter_by_merge_df_bacteria()

    def run_checkm(self) -> None:
        """运行 checkM. checkM 评估 bins 多个基因组"""
        # 创建 checkM 输入 bins 目录，复制并解压所有基因组
        all_dir = self.gnm_dir.joinpath("all")
        bins_dir = self.gnm_dir.joinpath("genome_assess/checkm_input")
        bins_dir.mkdir(parents=True, exist_ok=True)
        run(f"cp {all_dir}/*/*.fna {bins_dir}", shell=True, check=True)
        # 运行 checkM
        checkm_dir = self.gnm_dir.joinpath("genome_assess/checkm")
        checkm_dir.mkdir(parents=True, exist_ok=True)
        result_file = checkm_dir.joinpath("result.tsv")
        # * 如果结果文件已存在且不强制运行，则跳过
        if result_file.exists() and not self.force:
            logging.warning(f"checkM 结果文件 {result_file} 已存在, 跳过运行.")
            return
        # 并行太多 mamba run 会报文件锁的错误
        checkm_cmd = f"""
        source {ACTIVATE} qpcr
        checkm lineage_wf \
            -x fna --tab_table -t {self.threads} --pplacer_threads {self.threads} \
            -f {result_file} \
            {bins_dir} \
            {checkm_dir}
        conda deactivate
        """
        logging.debug(f"运行 checkM: {checkm_cmd}")
        run(checkm_cmd, shell=True, check=True, executable="/bin/bash")
        run(f"{CSVTK} -t csv2xlsx {result_file}", shell=True, check=True)

    def get_genome_anno_quality(self) -> None:
        """获取基因组注释质量"""
        # 输入 Prokka 注释结果的目录
        prk_dir = self.gnm_dir / "genome_annotate"
        # 获取去所有 GCF.tsv 特征文件中 rRNA, tRNA 数量
        rna_mtrx = []
        for feat_file in prk_dir.glob("*/*.tsv"):
            # genomes/Coxiella_Burnetii/genome_annotate/GCA_000300315.1/GCA_000300315.1.tsv
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
        rna_df = pd.DataFrame(
            rna_mtrx, columns=["genome_id", "23SrRNA", "16SrRNA", "5SrRNA", "tRNA"])
        # 保存为 CSV
        rna_df.to_csv(self.assess_dir / "rna_quality.csv", index=False)

    def merge_checkm_rna_stats(self) -> None:
        """
        合并 checkM 和 RNA 统计表
        :param assess_dir: 评估目录
        """
        # checkm 结果, 仅保留必须三列
        checkm_df = pd.read_csv(f"{self.assess_dir}/checkm/result.tsv", sep="\t",
                                usecols=["Bin Id", "Completeness", "Contamination"])
        checkm_df["genome_id"] = checkm_df.apply(
            lambda row: "_".join(row["Bin Id"].split("_")[:2]), axis=1)
        # 重排列表顺序
        checkm_df = checkm_df[["genome_id", "Completeness", "Contamination"]]
        # 重新读一下 RNA 统计表
        rna_df = pd.read_csv(self.assess_dir / "rna_quality.csv")
        # 合并
        merge_df = pd.merge(checkm_df, rna_df, on="genome_id")
        merge_df.to_excel(self.assess_dir / "checkm_rna_statistics.xlsx", index=False)
        merge_df.to_csv(self.assess_dir / "checkm_rna_statistics.csv", index=False)

    def filter_by_merge_df_bacteria(self):
        """
        根据 checkM 和 RNA 统计表合并结果过滤基因组, 生成 high_quality_genomes.txt 文件
        :param assess_dir: 评估目录
        """
        df = pd.read_csv(self.assess_dir / "checkm_rna_statistics.csv")
        # ! [250526 FJH] 过滤条件
        # 1. 基因组完整度 (Completeness) ≥ 90%
        # 2. 污染度 (Contamination) ≤ 5%
        # 3. 注释到23S rRNA基因、16S rRNA基因、5S rRNA基因和至少18个tRNA基因
        fltr_df = (
            df[df["Completeness"] > 90]
            .loc[df["Contamination"] < 5]
            .loc[df["23SrRNA"] > 0]
            .loc[df["16SrRNA"] > 0]
            .loc[df["5SrRNA"] > 0]
            .loc[df["tRNA"] >= 18]
        )
        fltr_df["genome_id"].to_csv(
            self.assess_dir / "high_quality_genomes.txt", index=False, header=False)


class GenomeQualityAssessorViruses(GenomeQualityAssessor):
    def __init__(self, sci_name: str, genome_set_dir: str, threads: int, force: bool = False):
        super().__init__(sci_name, genome_set_dir, threads, force)

    def run(self):
        """病毒基因组评估流程"""
        logging.info(f"开始病毒基因组质量评估: {self.gnm_dir}, 线程数: {self.threads}, 强制重新运行: {self.force}")
        self.run_checkv()
        self.filter_by_checkv()

    def run_checkv(self):
        """运行 checkV"""
        # 配置输出目录
        all_dir = self.gnm_dir.joinpath("all")
        checkv_dir = self.gnm_dir.joinpath("genome_assess/checkv")
        checkv_bins_dir = checkv_dir / "bins"
        checkv_bins_dir.mkdir(parents=True, exist_ok=True)
        # * 如果结果文件已存在且不强制运行，则跳过
        result_file = checkv_dir.joinpath("checkv_summary.tsv")
        if result_file.exists() and not self.force:
            logging.warning(f"checkV 结果文件 {result_file} 已存在, 跳过运行.")
            return
        # 搜索所有的 fna 文件, 写入批量运行脚本
        checkv_cmds = []
        for fna in all_dir.glob("*/*.fna"):
            checkv_cmd = f"source {ACTIVATE} qpcr && checkv end_to_end -t 1 -d {CHECKV_DB} {fna} {checkv_bins_dir}/{fna.parent.name} && conda deactivate"
            checkv_cmds.append(checkv_cmd)
        multi_run_command(checkv_cmds, self.threads)
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

    def filter_by_checkv(self):
        """根据 checkV 结果过滤基因组, 生成 high_quality_genomes.txt 文件"""
        checkv_dir = self.assess_dir / "checkv"
        out_gnm_file = self.assess_dir / "high_quality_genomes.txt"
        restab = checkv_dir / "checkv_summary.tsv"
        # 病毒 checkv
        df = pd.read_csv(restab, sep="\t", usecols=[
            "genome_id", "checkv_quality", "warnings", "completeness", "contamination"])
        # ! [250603 FJH] 过滤条件
        # 1.基因组质量(Checkv_Quality) 为高/中/低质量
        # 2.无warnings
        # 3.基因组完整度(Completeness) > 90%. PS.如为基因组为分段形式则以总和计算.
        # checkv_quality 1: high quality, 2: medium quality, 3: low quality
        to_move_genomes = df[~df['checkv_quality'].isin(
            ['Low-quality', 'Medium-quality', 'High-quality'])]['genome_id'].unique().tolist()
        df = df[~df['genome_id'].isin(to_move_genomes)]
        # warnings 和 contamination
        df = df[df['warnings'].isna()][["genome_id", "completeness"]]
        # 基因组完整度
        df = df.groupby('genome_id').sum()
        df[df['completeness'] > 90].index.to_series().to_csv(
            out_gnm_file, index=False, header=False)

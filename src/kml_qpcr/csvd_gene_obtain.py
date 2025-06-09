from pathlib import Path
from subprocess import run
import logging
import pandas as pd

from src.config.cnfg_software import ACTIVATE, SEQKIT, SEQTK
from src.utils.util_command import multi_run_command


class ConservedGenePredictor:
    def __init__(self, sci_name: str, genome_set_dir: str, threads: int, core_islt_perc: int, core_blastp_idnt, force: bool):
        """
        初始化保守基因预测器
        :sci_name: 目标微生物科学名
        :genome_set_dir: 基因组集合目录
        :threads: 线程数
        :core_islt_perc: 核心基因覆盖分离株百分比
        :core_blastp_idnt: 核心基因在分离株间相似度
        :force: 是否强制执行
        """
        self.sci_name = sci_name
        self.genome_set_dir = genome_set_dir
        self.threads = threads
        self.core_islt_perc = core_islt_perc
        self.core_blastp_idnt = core_blastp_idnt
        self.force = force
        # 当前微生物基因组目录
        self.gnm_dir = Path(genome_set_dir).joinpath(sci_name.replace(" ", "_"))
        # 保守基因目录
        self.csvd_dir = self.gnm_dir / "conserved_gene"
        self.csvd_dir.mkdir(exist_ok=True, parents=True)
        # 高质量分离株列表
        self.hq_gnms = self.get_high_quality_genomes()

    def run(self):
        """保守基因预测"""
        logging.info("获取保守基因")
        # 准备 Roary 输入文件
        self.prepare_roary_input()
        # 运行 Roary
        self.run_roary()
        # 筛选核心单拷贝基因
        core_sglcp_genes = self.filter_core_single_copy_gene()
        # 拆分所有高质量分离株的基因序列
        self.seqkit_split_isolates_ffn()
        # 输出所有保守基因序列合集
        self.output_conserved_gene_set(core_sglcp_genes)
        # 统计保守基因长度
        self.calc_gene_length(core_sglcp_genes)

    def get_high_quality_genomes(self) -> list[str]:
        """获取高质量分离株列表"""
        with open(self.gnm_dir / "genome_assess/high_quality_genomes.txt") as f:
            return f.read().splitlines()

    def prepare_roary_input(self) -> None:
        """准备 Roary 输入文件, 使用筛选后的分离株"""
        # 复制 Roary 需要的 .gff 文件到 Roary 输入目录
        roary_indir = self.csvd_dir / "roary_input"
        roary_indir.mkdir(exist_ok=True, parents=True)
        for gnm in self.hq_gnms:
            run(f"cp {self.gnm_dir}/genome_annotate/{gnm}/{gnm}.gff {roary_indir}/{gnm}.gff",
                shell=True, check=True)

    def run_roary(self) -> None:
        """
        运行 Roary 预测保守基因, 核心文件:
        - clustered_proteins: 基因覆盖的基因组, 计算单拷贝基因. 
        - gene_presence_absence.Rtab: 基因和基因组的对照表, 筛选出核心基因
        """
        roary_dir = self.csvd_dir / "roary"
        # 如果已经存在或非强制, 跳过
        if (roary_dir.joinpath("gene_presence_absence.Rtab").exists() and
            roary_dir.joinpath("clustered_proteins").exists() and
                (not self.force)):
            logging.warning(f"Roary 已经运行过 {roary_dir}, 跳过.")
            return
        # 如果存在 Roary 结果目录无论是否完整, Roary 都会认为是已经运行过. 结果会生成到其他文件夹, 需要删除
        if roary_dir.exists():
            run(f"rm -r {roary_dir}", shell=True, check=True)
        # * [250603 FJH] 核心基因样本比例调整至 100%; 相似度也调整至 100%
        cmd = f"""
        source {ACTIVATE} meta
        roary -p {self.threads} -cd {self.core_islt_perc} -i {self.core_blastp_idnt} {self.csvd_dir}/roary_input/*.gff -f {roary_dir}
        conda deactivate
        """
        logging.debug(cmd)
        run(cmd, shell=True, check=True, executable="/bin/bash")

    def filter_core_single_copy_gene(self) -> list[str]:
        """
        过滤核心单拷贝基因. 1. 覆盖超过阈值的核心基因; 2. 单拷贝基因; 3. 有名称的基因 (删掉group_*)
        :return: 核心单拷贝基因列表
        """
        # 基因存在缺失矩阵
        df = pd.read_table(self.csvd_dir / "roary/gene_presence_absence.Rtab",
                           sep='\t', index_col=0)
        # 1. 覆盖超过阈值的核心基因
        core_genes = df[df.apply(lambda x: x.sum() / len(x) * 100, axis=1)
                        >= self.core_islt_perc].index.tolist()
        # 分离株总数
        islt_count = df.shape[1]
        # 2. 单拷贝基因
        sglcp_genes = []
        # 核心基因明细文件
        with open(self.csvd_dir / "roary/clustered_proteins") as f:
            for line in f:
                gene, gene_cntt = line.strip().split(': ')
                islts = ["_".join(gid.split("_")[:2]) for gid in gene_cntt.split("\t")]
                # 基因数大于分离株总数 或 基因占比小于核心基因占比
                if (len(islts) > islt_count) or (len(islts) / islt_count * 100 < self.core_islt_perc):
                    continue
                # 基因数等于分离株总数, 但不是单拷贝
                if len(islts) != len(set(islts)):
                    continue
                sglcp_genes.append(gene)
        # 1&2 核心单拷贝基因列表
        core_sglcp_itsct = list(set(core_genes).intersection(set(sglcp_genes)))
        # 3. 有名称的基因 (删掉group_*)
        # ! 基因名有路径分隔符号 '/', 引起报错. 例如 pgk/tpi.
        core_sglcp_genes = []
        for gene in core_sglcp_itsct:
            if gene.startswith("group_") or ("/" in gene):
                continue
            core_sglcp_genes.append(gene)
        with open(self.csvd_dir / "core_single_copy_genes.txt", "w") as f:
            f.write("\n".join(core_sglcp_genes))
        return core_sglcp_genes

    def seqkit_split_isolates_ffn(self):
        """Seqkit 按照序列 ID 拆分所有高质量基因组注释基因文件 .ffn"""
        all_genes_dir = self.csvd_dir / "all_genes"
        all_genes_dir.mkdir(exist_ok=True, parents=True)
        # 是否强制执行
        if len(list(all_genes_dir.iterdir())) >= len(self.hq_gnms) and not self.force:
            logging.warning(f"Seqkit 已拆分完所有分离株基因 {all_genes_dir}, 跳过.")
            return
        cmds = []
        for gnm in self.hq_gnms:
            cmd = f"{SEQKIT} split --force --by-id --by-id-prefix '' {self.gnm_dir}/genome_annotate/{gnm}/{gnm}.ffn -O {all_genes_dir}/{gnm}"
            cmds.append(cmd)
        multi_run_command(cmds, self.threads)

    def output_conserved_gene_set(self, core_sglcp_genes: list[str]):
        """输出保守基因序列合集"""
        # 保守基因序列合集目录
        csvd_gene_set_dir = self.csvd_dir / "csvd_gene_seq_set"
        csvd_gene_set_dir.mkdir(exist_ok=True, parents=True)
        cat_cmds = []
        with open(self.csvd_dir / "roary/clustered_proteins") as f:
            for line in f:
                gene, gene_cntt = line.strip().split(': ')
                # 是否为单拷贝核心基因
                if gene in core_sglcp_genes:
                    cat_files = []
                    for gnid in gene_cntt.split("\t"):
                        isltid = "_".join(gnid.split("_")[:2])
                        cat_files.append(f"{self.csvd_dir}/all_genes/{isltid}/{gnid}.ffn")
                    cmd = f"cat {' '.join(cat_files)} > {csvd_gene_set_dir}/{gene}.ffn && " \
                        f"{SEQTK} comp {csvd_gene_set_dir}/{gene}.ffn > {csvd_gene_set_dir}/{gene}.ffn.comp.tsv"
                    cat_cmds.append(cmd)
        multi_run_command(cat_cmds, self.threads)

    def calc_gene_length(self, core_sglcp_genes: list[str]):
        """统计保守基因长度. 中位数,极差,标准差,变异系数"""
        gene_length_stats_mat = []
        for gene in core_sglcp_genes:
            comp_file = self.csvd_dir / "csvd_gene_seq_set" / f"{gene}.ffn.comp.tsv"
            df = pd.read_table(comp_file, sep="\t", header=None, usecols=[1])
            gl = df.iloc[:, 0]
            gene_length_stats_mat.append(
                # 基因名, 中位数, 极差, 标准差, 变异系数
                [gene, gl.median(), gl.max() - gl.min(), gl.std(), gl.std() / gl.mean()])
        odf = pd.DataFrame(gene_length_stats_mat, columns=["gene", "median", "range", "std", "cv"])
        odf.to_csv(self.csvd_dir / "gene_length_stat.csv", index=False)
        odf.to_excel(self.csvd_dir / "gene_length_stat.xlsx", index=False)

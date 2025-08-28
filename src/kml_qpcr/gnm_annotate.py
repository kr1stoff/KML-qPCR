from pathlib import Path
import logging

from src.kml_qpcr.base import BaseQPCR
from src.config.cnfg_software import ACTIVATE
from src.utils.util_command import multi_run_command


class GenomeAnnotator(BaseQPCR):
    def __init__(self, sci_name: str, genome_set_dir: str, threads: int, force: bool):
        """
        注释基因组.
        :param sci_name: 物种学名.
        :param genome_set_dir: 基因组集目录.
        :param threads: 全局线程数.
        :param force: 是否强制重新运行 prokka, 默认识别到结果文件就跳过.
        """
        super().__init__(sci_name, genome_set_dir, threads, force)
        # 存放注释结果的文件夹
        self.gnm_annt_dir = self.gnm_dir / "genome_annotate"
        self.gnm_annt_dir.mkdir(parents=True, exist_ok=True)

    def run(self) -> None:
        logging.info(f"开始注释基因组: {self.gnm_dir}, 线程数: {self.threads}, 强制重新运行: {self.force}")
        # 运行 prokka
        self.run_prokka()

    def run_prokka(self):
        """使用 Prokka 进行基因组注释."""
        # 当前物种目录下 all 文件夹
        all_dir = self.gnm_dir / "all"
        # 如果不是强制重新运行, 并且 prokka 已经运行完成, 则跳过注释
        gff_count = len(list(self.gnm_annt_dir.glob("*/*.gff")))
        fna_count = len(list(all_dir.glob("*/*.fna")))
        if (not self.force) and (gff_count >= fna_count):
            logging.warning("Prokka 已经运行完成, 跳过注释.")
            return
        # 批量运行 prokka
        # 和后面统一不用 PARALELL 用 multiprocessing.Pool 替代
        prk_cmds = []
        for fna in all_dir.glob("*/*.fna"):
            gnm_id = fna.parent.name
            # ! Roary 需要每个 GFF 文件 basename 不同. Error: GFF files must have unique basenames
            prk_cmd = f"source {ACTIVATE} meta && prokka --cpu 1 --force --prefix {gnm_id} --outdir {self.gnm_annt_dir}/{gnm_id} --kingdom Bacteria --addgenes --quiet --locustag {gnm_id} {fna} && conda deactivate"
            prk_cmds.append(prk_cmd)
        multi_run_command(prk_cmds, self.threads)

# todo 病毒注释

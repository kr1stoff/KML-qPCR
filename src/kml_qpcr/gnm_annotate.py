from pathlib import Path
from subprocess import run
import logging

from src.config.cnfg_software import ACTIVATE, PARALLEL
from src.utils.util_file import count_file_lines


def annotate_genomes(sci_name: str, genome_set_dir: str, threads: int, force: bool) -> None:
    """
    注释基因组.
    :param sci_name: 物种学名.
    :param genome_set_dir: 基因组集目录.
    :param threads: 全局线程数.
    :param force: 是否强制重新运行 prokka, 默认识别到结果文件就跳过.
    :return: None
    """
    logging.info(f"开始注释基因组: {sci_name}, 线程数: {threads}, 强制重新运行: {force}")
    # TODO 病毒
    sci_name = sci_name.replace(" ", "_")
    # 当前物种目录
    gnm_dir = Path(genome_set_dir) / sci_name
    # 当前物种目录下 all 文件夹
    all_dir = gnm_dir / "all"
    # 存放注释结果的文件夹
    gnm_annt_dir = gnm_dir / "genome_annotate"
    gnm_annt_dir.mkdir(parents=True, exist_ok=True)
    # 运行 prokka
    run_prokka(all_dir, gnm_annt_dir, threads, force)


def run_prokka(all_dir: Path, gnm_annt_dir: Path, threads: int, force: bool):
    """
    使用 Prokka 进行基因组注释.
    :param all_dir: 包含所有基因组的目录.
    :param gnm_annt_dir: 注释结果输出目录.
    :param threads: 使用的线程数.
    """
    # 并行数看情况修改
    prl_num = min(4, threads)
    sgl_thrds = threads // prl_num
    # ! all 目录里面只放基因组文件夹, 不要放别的
    with open(f"{gnm_annt_dir}/prokka_batch.sh", "w") as f:
        for cur_gnm_dir in all_dir.iterdir():
            gnm_num = cur_gnm_dir.name
            fna = list(cur_gnm_dir.glob("*.fna"))[0]
            prk_cmd = f"prokka --cpu {sgl_thrds} --force --prefix prokka --outdir {gnm_annt_dir}/{gnm_num} --kingdom Bacteria --addgenes --quiet --locustag {gnm_num} {fna}"
            f.write(prk_cmd + "\n")
    # 进入 meta 环境运行 prokka_batch.sh
    prk_batch_cmd = f"""
source {ACTIVATE} meta
cat {gnm_annt_dir}/prokka_batch.sh | {PARALLEL} -j {prl_num}
conda deactivate
"""
    if force or not prokka_is_complete_run(gnm_annt_dir):
        run(prk_batch_cmd, shell=True, check=True, executable="/bin/bash")
    else:
        logging.info("Prokka 已经运行完成, 跳过注释.")


def prokka_is_complete_run(gnm_annt_dir: Path) -> bool:
    """
    检查 prokka 是否已经运行完成.
    :param gnm_annt_dir: 注释结果输出目录.
    :return: 如果 prokka 结果文件存在则返回 True, 否则返回 False.
    """
    gnm_count = count_file_lines(f"{gnm_annt_dir}/prokka_batch.sh")
    result_files = list(gnm_annt_dir.glob("*/prokka.gff"))
    return len(result_files) == gnm_count

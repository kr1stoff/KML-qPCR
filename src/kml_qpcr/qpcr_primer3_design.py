import click
import logging
from pathlib import Path
from Bio import SeqIO
from subprocess import run
from multiprocessing import Pool

from src.config.cnfg_software import MAMBA


logging.basicConfig(level=logging.DEBUG,
                    format="%(asctime)s - %(levelname)s - %(filename)s - %(message)s",
                    datefmt="%Y-%m-%d %H:%M:%S")


@click.command()
@click.option("--fasta", required=True, help="输入用于设计引物的FASTA序列文件.")
@click.option("--workdir", default="qpcr_analysis", show_default=True, help="引物分析工作目录.")
@click.option("--threads", default=8, type=int, show_default=True, help="最大线程数.")
@click.help_option(help="显示帮助信息.")
def main(fasta, workdir, threads):
    """
    批量设计 qPCR 引探
    输入多区域 fasta 序列文件, 拆分序列写入到多个配置文件并行跑 primer3, 输出 primer3 结果
    """
    logging.info("开始运行 qPCR 引物设计")
    home = Path(__file__).resolve().parents[1]
    template = home.joinpath("config/template.p3")
    shell_dir = Path(workdir).resolve().joinpath("shell")
    shell_dir.mkdir(parents=True, exist_ok=True)

    split_and_write_primer_inputs(fasta, workdir, template, shell_dir, threads)
    parallel_run_primer3(shell_dir, threads)
    logging.info("完成 qPCR 引物设计")


def split_and_write_primer_inputs(fasta, workdir, template, shell_dir, threads):
    """
    拆分序列写入到多个配置文件, 多个区域 primer3 命令行脚本写入到 batch 文件中.

    :param fasta: 输入用于设计引物的FASTA序列文件.
    :param workdir: 引物分析工作目录.
    :param template: primer3 配置模板文件.
    :param shell_dir: shell 脚本目录.
    :param threads: 使用的 CPU 核心数.
    """
    primer3_indir = Path(workdir).resolve().joinpath("primer3_in")
    primer3_indir.mkdir(parents=True, exist_ok=True)
    primer3_outdir = Path(workdir).resolve().joinpath("primer3_out")
    primer3_outdir.mkdir(parents=True, exist_ok=True)

    records = SeqIO.parse(fasta, "fasta")
    # Pool.map 并行写入 primer3 输入文件的列表
    write_primer3_input_tuple = []
    # 并行运行 primer3 的 shell 列表, 写入到 batch.sh 给 parallel 运行
    para_shells = []
    # primer input 模板, 不在 for 中使用, 频繁 IO 速度慢
    primer3_template = open(template).read()

    for record in records:
        name = record.name.replace(":", "_")    # : 替换为 _, 避免文件名冲突
        seq = record.seq
        write_primer3_input_tuple.append((primer3_indir, name, seq, primer3_template))
        para_shells.append(
            f"{MAMBA} -n qpcr run primer3_core < {primer3_indir}/{name}.p3 > {primer3_outdir}/{name}.out")

    # 并行写 primer3 输入文件
    with Pool(threads) as pool:
        pool.map(parallel_write_primer3_input, write_primer3_input_tuple)

    # * 写个多个区域 primer3 命令行到 batch 文件中
    with open(f"{shell_dir}/primer3_batch.sh", "w") as f:
        for ps in para_shells:
            f.write(ps + "\n")


def parallel_write_primer3_input(input_tuple):
    """
    并行写入 primer3 输入文件的方法
    """
    primer3_indir, name, seq, primer3_template = input_tuple
    # 写个多个区域配置文件
    with open(f"{primer3_indir}/{name}.p3", "w") as g:
        g.write(f"SEQUENCE_ID={name}\nSEQUENCE_TEMPLATE={seq}\n{primer3_template}")


def parallel_run_primer3(shell_dir, threads):
    """
    并行运行 primer3.

    :param shell_dir: shell 脚本目录.
    :param threads: 并行数.
    """
    cmd = f"cat {shell_dir}/primer3_batch.sh | {MAMBA} -n basic run parallel -j {threads}"
    logging.debug(f"运行 primer3 命令: {cmd}")
    res = run(cmd, shell=True, executable="/bin/bash", capture_output=True, encoding="utf-8")

    with open(f"{shell_dir}/primer3_batch.sh.log", "w") as f:
        f.write(f"[stdout]\n{res.stdout}\n[stderr]\n{res.stderr}]\n")


if __name__ == "__main__":
    main()

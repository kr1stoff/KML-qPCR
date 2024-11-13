import yaml
import click
import logging
from pathlib import Path
from Bio import SeqIO
from subprocess import run
from multiprocessing import Pool


logging.basicConfig(level=logging.DEBUG,
                    format="%(asctime)s - %(levelname)s - %(filename)s - %(message)s",
                    datefmt="%Y-%m-%d %H:%M:%S")


@click.command()
@click.option('-i', '--fasta', type=click.Path(exists=True), help='输入用于设计引物的FASTA序列文件.')
@click.option('-c', '--cores', default=8, type=int, help='使用的 CPU 核心数.')
@click.option('-w', '--workdir', default='qpcr_analysis', help='引物分析工作目录.')
@click.help_option('-h', '--help')
def main(fasta, workdir, cores):
    """
    批量设计 qPCR 引探
    输入多区域 fasta 序列文件, 拆分序列写入到多个配置文件并行跑 primer3, 输出 primer3 结果
    """
    logging.info('开始运行 qPCR 引物设计')
    logging.info(f'输入文件: {fasta}')
    logging.info(f'工作目录: {workdir}')
    logging.info(f'并行数: {cores}')
    home = Path(__file__).resolve().parents[1]
    template = home.joinpath('config/template.p3')
    primer3 = yaml.safe_load(open(home.joinpath('config/software.yaml')).read())['primer3']
    parallel = yaml.safe_load(open(home.joinpath('config/software.yaml')).read())['parallel']
    shell_dir = Path(workdir).resolve().joinpath('shell')
    shell_dir.mkdir(parents=True, exist_ok=True)

    split_and_write_primer_inputs(fasta, workdir, template, primer3, shell_dir, cores)
    parallel_run_primer3(shell_dir, workdir, cores, parallel)
    logging.info('完成 qPCR 引物设计')


def split_and_write_primer_inputs(fasta, workdir, template, primer3, shell_dir, cores):
    """
    拆分序列写入到多个配置文件, 多个区域 primer3 命令行脚本写入到 batch 文件中.
    该步骤 primer3_in.DONE 文件存在则跳过.

    :param fasta: 输入用于设计引物的FASTA序列文件.
    :param workdir: 引物分析工作目录.
    :param template: primer3 配置模板文件.
    :param primer3: primer3 命令行.
    :param shell_dir: shell 脚本目录.
    :param cores: 使用的 CPU 核心数.
    """
    logging.info('拆分序列写入到多个配置文件')
    primer3_indir = Path(workdir).resolve().joinpath('primer3_in')
    primer3_indir.mkdir(parents=True, exist_ok=True)
    primer3_outdir = Path(workdir).resolve().joinpath('primer3_out')
    primer3_outdir.mkdir(parents=True, exist_ok=True)
    primer3_indir_done_flag = primer3_indir.resolve().joinpath('primer3_in.DONE')

    if primer3_indir_done_flag.exists():
        logging.info('primer3_in 目录已存在, 跳过拆分序列写入到多个配置文件')
        return

    records = SeqIO.parse(fasta, 'fasta')
    # Pool.map 并行写入 primer3 输入文件的列表
    write_primer3_input_tuple = []
    # 并行运行 primer3 的 shell 列表, 写入到 batch.sh 给 parallel 运行
    para_shells = []
    # primer input 模板, 不在 for 中使用, 频繁 IO 速度慢
    primer3_template = open(template).read()

    for record in records:
        name = record.name.replace(':', '_')    # : 替换为 _, 避免文件名冲突
        seq = record.seq
        write_primer3_input_tuple.append((primer3_indir, name, seq, primer3_template))
        para_shells.append(f'{primer3} < {primer3_indir}/{name}.p3 > {primer3_outdir}/{name}.out')

    # 并行写 primer3 输入文件
    with Pool(cores) as pool:
        pool.map(parallel_write_primer3_input, write_primer3_input_tuple)

    # * 写个多个区域 primer3 命令行到 batch 文件中
    with open(f'{shell_dir}/primer3_batch.sh', 'w') as f:
        for ps in para_shells:
            f.write(ps + '\n')

    # 创建 primer3_in.DONE 文件
    primer3_indir_done_flag.touch()


def parallel_write_primer3_input(input_tuple):
    """
    并行写入 primer3 输入文件的方法
    """
    primer3_indir, name, seq, primer3_template = input_tuple
    # 写个多个区域配置文件
    with open(f'{primer3_indir}/{name}.p3', 'w') as g:
        g.write(f'SEQUENCE_ID={name}\nSEQUENCE_TEMPLATE={seq}\n{primer3_template}')


def parallel_run_primer3(shell_dir, workdir, cores, parallel):
    """
    并行运行 primer3.
    primer3_out.DONE 文件存在则跳过.

    :param shell_dir: shell 脚本目录.
    :param cores: 并行数.
    :param parallel: parallel 命令行.
    """
    logging.info('并行运行 primer3')
    primer3_outdir_done_flag = Path(workdir).resolve().joinpath('primer3_out/primer3_out.DONE')
    if primer3_outdir_done_flag.exists():
        logging.info('primer3_out 目录已存在, 跳过并行运行 primer3')
        return

    cmd = f'cat {shell_dir}/primer3_batch.sh | {parallel} -j {cores}'
    res = run(cmd, shell=True, executable='/bin/bash', capture_output=True, encoding='utf-8')

    with open(f'{shell_dir}/primer3_batch.sh.log', 'w') as f:
        f.write(f'[stdout]\n{res.stdout}\n[stderr]\n{res.stderr}]\n')

    # 创建 primer3_out.DONE 文件
    primer3_outdir_done_flag.touch()


if __name__ == '__main__':
    main()

import re
import click
import logging
import sys
from pathlib import Path
from multiprocessing import Pool


logging.basicConfig(level=logging.DEBUG,
                    format="%(asctime)s - %(levelname)s - %(filename)s - %(message)s",
                    datefmt="%Y-%m-%d %H:%M:%S")


@click.command()
@click.option('-w', '--workdir', default='qpcr_analysis', help='引物分析工作目录.')
@click.option('-c', '--cores', default=8, type=int, help='使用的 CPU 核心数.')
@click.help_option('-h', '--help')
def main(workdir, cores):
    """
    输入输出同目录，接上 qpcr_primer3_design.py 的工作目录，解析 primer3 输出引物表格
    """
    logging.info('开始解析 primer3 输出引物表格')
    logging.info(f'工作目录: {workdir}')
    primer3_outdir = Path(workdir).resolve().joinpath('primer3_out')
    primer3_pardir = Path(workdir).resolve().joinpath('primer3_parse')
    primer3_pardir.mkdir(parents=True, exist_ok=True)
    primer3_outdir_done_flag = primer3_pardir.resolve().joinpath('primer3_parse.DONE')

    # 检查输入
    check_input(primer3_outdir, primer3_outdir_done_flag)

    out_par_tuple = []
    for p3out in Path(primer3_outdir).glob('*.out'):
        p3par = primer3_pardir.joinpath(p3out.name.replace('.out', '.tsv'))
        out_par_tuple.append((p3out, p3par))

    # * 多进程解析 primer3 输出表格, tqdm 显示进度条
    with Pool(cores) as pool:
        pool.map(parse_primer3_out, out_par_tuple)

    primer3_outdir_done_flag.touch()
    logging.info('解析 primer3 输出表格完成')


def check_input(primer3_outdir, primer3_outdir_done_flag):
    """
    检查输入, 1. 输入不存在； 2. 已解析完成

    :param primer3_outdir: primer3 输出目录
    :param primer3_outdir_done_flag: primer3 输出目录解析完成标志
    """
    if not primer3_outdir.exists():
        logging.error(f'primer3_out 目录不存在, 无法解析 primer3 结果')
        sys.exit(1)

    if primer3_outdir_done_flag.exists():
        logging.info(f'primer3_parse 目录已存在, 跳过解析 primer3 输出引物表格')
        sys.exit(0)


def parse_sequence_id(seqid) -> tuple:
    """
    解析sequence_id, 符合格式`chrom:start-end`则返回 [chrom, start, end]
    例子: SEQUENCE_ID=D90400.1_sliding_4261-4560

    :param seqid: sequence_id
    :return: [chrom, start, end]
    """
    ptn = re.compile(r'(.*?)_sliding_(\d+)-(\d+)')
    if re.match(ptn, seqid):
        res = re.findall(ptn, seqid)[0]
        # ! seqkit sliding : 1-bases
        return (res[0], str(int(res[1]) - 1), res[2])
    else:
        return (seqid, '-', '-')


def calc_ampl_gc(amplicon_seq) -> float:
    """
    计算引物的 GC 含量
    """
    ampl_gc = 0
    for base in amplicon_seq:
        if base in ['G', 'C']:
            ampl_gc += 1
    return ampl_gc / len(amplicon_seq)


def parse_primer3_out(input_tuple):
    """
    解析 primer3 输出表格

    :param primer3_out: primer3 输出表格
    :param parse_tab: 解析后的表格
    """
    primer3_out, parse_tab = input_tuple
    content = open(primer3_out, 'r').read()
    sequence_id = re.findall(r"SEQUENCE_ID=(\S+)", content)[0]
    sequence_template = re.findall(r"SEQUENCE_TEMPLATE=(\S+)", content)[0]
    chrom, start, end = parse_sequence_id(sequence_id)
    # ! primer3 : 0-based
    forward_pos = re.findall(r"PRIMER_LEFT_\d+=(\d+),(\d+)", content)
    forward_seq = re.findall(r"PRIMER_LEFT_\d+_SEQUENCE=(\S+)", content)
    forward_tm = re.findall(r"PRIMER_LEFT_\d+_TM=(\S+)", content)
    forward_gc = re.findall(r"PRIMER_LEFT_\d+_GC_PERCENT=(\S+)", content)
    reverse_pos = re.findall(r"PRIMER_RIGHT_\d+=(\d+),(\d+)", content)
    reverse_seq = re.findall(r"PRIMER_RIGHT_\d+_SEQUENCE=(\S+)", content)
    reverse_tm = re.findall(r"PRIMER_RIGHT_\d+_TM=(\S+)", content)
    reverse_gc = re.findall(r"PRIMER_RIGHT_\d+_GC_PERCENT=(\S+)", content)
    amplicon_len = re.findall(r"PRIMER_PAIR_\d+_PRODUCT_SIZE=(\S+)", content)
    amplicon_tm = re.findall(r"PRIMER_PAIR_\d+_PRODUCT_TM=(\S+)", content)
    probe_pos = re.findall(r"PRIMER_INTERNAL_\d+=(\d+),(\d+)", content)
    probe_seq = re.findall(r"PRIMER_INTERNAL_\d+_SEQUENCE=(\S+)", content)
    probe_tm = re.findall(r"PRIMER_INTERNAL_\d+_TM=(\S+)", content)
    probe_gc = re.findall(r"PRIMER_INTERNAL_\d+_GC_PERCENT=(\S+)", content)

    with open(parse_tab, 'w') as f:
        header = 'index\tchromosome\tforward_start\tforward_end\tforward_length\tforward_tm' \
            '\tforward_gc\tforward_sequence\treverse_start\treverse_end\treverse_length\treverse_tm' \
            '\treverse_gc\treverse_sequence\tprobe_start\tprobe_end\tprobe_length\tprobe_tm\tprobe_gc' \
            '\tprobe_sequence\tamplicon_tm\tamplicon_gc\tamplicon_length\tamplicon_sequence\n'
        f.write(header)

        for i in range(len(forward_pos)):
            amplicon_seq_curr = sequence_template[int(forward_pos[i][0])-1: int(reverse_pos[i][0])]
            amplicon_gc_curr = calc_ampl_gc(amplicon_seq_curr)

            # * 输出结果使用 samtools faidx 1-based 坐标
            outlist = [
                str(i),
                chrom,
                str(int(start) + int(forward_pos[i][0])),
                # ! 根据 samtools 1-based 矫正坐标
                str(int(start) + int(forward_pos[i][0]) + int(forward_pos[i][1]) - 1),
                forward_pos[i][1],
                forward_tm[i],
                forward_gc[i],
                forward_seq[i],
                str(int(start) + int(reverse_pos[i][0])),
                # ! Reverse 起始位置在右端，需要减去长度
                str(int(start) + int(reverse_pos[i][0]) - int(reverse_pos[i][1]) + 1),
                reverse_pos[i][1],
                reverse_tm[i],
                reverse_gc[i],
                reverse_seq[i],
                str(int(start) + int(probe_pos[i][0])),
                str(int(start) + int(probe_pos[i][0]) + int(probe_pos[i][1]) - 1),
                probe_pos[i][1],
                probe_tm[i],
                probe_gc[i],
                probe_seq[i],
                amplicon_tm[i],
                str(amplicon_gc_curr),
                amplicon_len[i],
                amplicon_seq_curr
            ]

            f.write('\t'.join(outlist) + '\n')


if __name__ == '__main__':
    main()

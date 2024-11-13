import csv
import json
import click
import logging
import sys
from pathlib import Path
from collections import namedtuple
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
    包容性评估, 引探匹配的基因组放在最前面

    :param workdir: 引物分析工作目录.
    :param cores: 使用的 CPU 核心数.
    """
    logging.info('开始评估引物包容性')
    blast_out = Path(workdir).resolve().joinpath('blast/primers.out')
    primer_json = Path(workdir).resolve().joinpath('blast/primers.json')
    primer3_pardir = Path(workdir).resolve().joinpath('primer3_parse')
    incl_dir = Path(workdir).resolve().joinpath('inclusivity')
    incl_dir.mkdir(exist_ok=True)
    incl_dir_done_flag = incl_dir.resolve().joinpath('inclusivity.DONE')

    # * 检查输入
    check_input(primer3_pardir, incl_dir_done_flag)

    # 引物 名称和序列 的映射
    prim_name_dict = json.load(open(primer_json))

    inclu_input_tuple = []
    for parse_tab in primer3_pardir.glob('*.tsv'):
        inclu_input_tuple.append((parse_tab, blast_out, prim_name_dict,
                                 incl_dir.joinpath(parse_tab.name)))
    # * 并行分析
    with Pool(cores) as pool:
        pool.map(evaluate_inclusivity, inclu_input_tuple)

    # * 合并引物列表, 去除空文件
    merge_primer_table(workdir)

    incl_dir_done_flag.touch()
    logging.info('评估引物包容性完成')


def merge_primer_table(workdir):
    """
    合并引物列表, 去除空文件
    """
    outlist = []
    for incl_tab in Path(workdir).resolve().joinpath('inclusivity').glob('*.tsv'):
        with open(incl_tab) as f:
            reader = csv.DictReader(f, delimiter='\t')
            # ! 直接 len(list(reader)) 会迭代完条目, 后面就写不进 outlist 了, 所以先转成 list
            reader_list = list(reader)
            if len(reader_list) == 0:
                continue
            for row in reader_list:
                outlist.append(row.values())

    with open(Path(workdir).resolve().joinpath('inclusivity.tsv'), 'w') as f:
        f.write('\t'.join(reader.fieldnames) + '\n')
        for it in outlist:
            f.write('\t'.join(it) + '\n')


def check_input(primer3_pardir, incl_dir_done_flag):
    """
    检查输入
    """
    if not primer3_pardir.exists():
        logging.error(f'primer3_parse 目录不存在, 无法解析 primer3 结果')
        sys.exit(1)

    if incl_dir_done_flag.exists():
        logging.info(f'引物包容性评估目录已存在, 跳过评估引物包容性')
        sys.exit(0)


def get_namedtuple_index_header():
    """
    解析 primer3_parse 输出的 primer 表格，返回 namedtuple index-header
    """
    # qseqid sseqid length qlen slen sstart send qstart qend qcovs pident nident evalue bitscore sstrand
    headers = ['qseqid', 'sseqid', 'length', 'qlen', 'slen', 'sstart', 'send',
               'qstart', 'qend', 'qcovs', 'pident', 'nident', 'evalue', 'bitscore', 'sstrand']
    dict_index_column = {col: idx for idx, col in enumerate(headers)}
    nih = namedtuple('NIH', dict_index_column.keys())(**dict_index_column)
    return nih


def evaluate_inclusivity(input_tuple):
    parse_tab, blast_out, prim_name_dict, incl_out = input_tuple
    nih = get_namedtuple_index_header()
    # 阈值
    min_prod_len = 70
    max_prod_len = 1000
    # 所有引物，输出引物名称
    primer_names = []
    with open(parse_tab) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            primer_names.append(prim_name_dict[row['forward_sequence']])
            primer_names.append(prim_name_dict[row['reverse_sequence']])
            primer_names.append(prim_name_dict[row['probe_sequence']])

    # 在 blast 里面把上面的引物比对信息提取出来
    blast_dict = {}
    with open(blast_out) as f:
        for line in f:
            lns = line.strip().split('\t')
            # 只保留当前 parse table 的比对结果
            if lns[nih.qseqid] not in primer_names:
                continue
            # ! 过滤 3 端错配超过 2bp
            if int(lns[nih.qlen]) - int(lns[nih.qend]) > 2:
                continue
            blast_dict.setdefault(lns[nih.qseqid], [])
            blast_dict[lns[nih.qseqid]].append(lns)

    # 过滤并输出结果
    with open(parse_tab) as f, open(incl_out, 'w') as g:
        reader = csv.DictReader(f, delimiter='\t')
        g.write('\t'.join(['hits'] + reader.fieldnames) + '\n')

        for row in reader:
            inclu_nucles = []
            forward_hits = blast_dict[prim_name_dict[row['forward_sequence']]]
            reverse_hits = blast_dict[prim_name_dict[row['reverse_sequence']]]
            probe_hits = blast_dict[prim_name_dict[row['probe_sequence']]]
            # 满足的引物对
            primer_pairs = []
            for fh in forward_hits:
                sseqid_f, sstart_f, sstrand_f = fh[nih.sseqid], int(fh[nih.sstart]), fh[nih.sstrand]
                # ! 保留上游引物正向比对
                if sstrand_f != 'plus':
                    continue
                # * blast minus, 起始位置和终止位置按照 5->3, sstart > send
                for rh in reverse_hits:
                    sseqid_r, send_r, sstrand_r = rh[nih.sseqid], int(rh[nih.send]), rh[nih.sstrand]
                    # ! 保留下游引物反向比对
                    if sstrand_r != 'minus':
                        continue
                    # ! 保留上下游引物在同一条染色体上，并且长度在阈值内
                    if (sseqid_f == sseqid_r) and (min_prod_len < (send_r - sstart_f) < max_prod_len):
                        primer_pairs.append([sseqid_f, int(fh[nih.send]), send_r])

            for pp in primer_pairs:
                for ph in probe_hits:
                    pp_chrom, pp_start, pp_end = pp
                    sseqid_p, sstart_p, send_p, sstrand_p = ph[nih.sseqid], int(
                        ph[nih.sstart]), int(ph[nih.send]), ph[nih.sstrand]
                    #! 保留探针正向比对
                    if sstrand_p != 'plus':
                        continue
                    # ! 染色体相同, F引物右端 --- 探针 --- R引物左端
                    if (sseqid_p == pp_chrom) and (sstart_p > pp_start) and (send_p < pp_end):
                        inclu_nucles.append(pp_chrom)

            if len(inclu_nucles) > 0:
                # * inclu_nucles 去一下重复
                g.write('\t'.join([','.join(set(inclu_nucles))] + list(row.values())) + '\n')


if __name__ == '__main__':
    main()

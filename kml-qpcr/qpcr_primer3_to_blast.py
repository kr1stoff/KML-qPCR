from pathlib import Path
import sys
import click
import csv
import logging
import json


logging.basicConfig(level=logging.DEBUG,
                    format="%(asctime)s - %(levelname)s - %(filename)s - %(message)s",
                    datefmt="%Y-%m-%d %H:%M:%S")


@click.command()
@click.option('-w', '--workdir', default='qpcr_analysis', help='引物分析工作目录.')
@click.help_option('-h', '--help')
def main(workdir):
    """
    去除序列重复的引物. 分别导出 fasta 和 json 文件
    """
    logging.info('开始准备 blast 分析输入')
    blast_dir = Path(workdir).joinpath('blast')
    blast_dir.mkdir(exist_ok=True)
    primer3_parse_dir = Path(workdir).joinpath('primer3_parse')

    if not primer3_parse_dir.exists():
        logging.error(f'目录 {workdir} 下没有 primer3_parse 文件夹, 请检查.')
        sys.exit(1)

    # 引物字典
    prim_dict = {}
    counter = 0
    for ptsv in primer3_parse_dir.glob('*.tsv'):
        with open(ptsv) as file:
            reader = csv.DictReader(file, delimiter='\t')
            for row in reader:
                prim_dict, counter = insert_seq_to_dict(row['forward_sequence'], prim_dict, counter)
                prim_dict, counter = insert_seq_to_dict(row['reverse_sequence'], prim_dict, counter)
                prim_dict, counter = insert_seq_to_dict(row['probe_sequence'], prim_dict, counter)

    # 写 fasta
    with open(blast_dir.joinpath('primers.fasta'), 'w') as f:
        for seq, name in prim_dict.items():
            f.write(f'>{name}\n{seq}\n')

    # 写 json
    with open(blast_dir.joinpath('primers.json'), 'w') as f:
        json.dump(prim_dict, f)

    logging.info('完成')


def insert_seq_to_dict(seq, seq_dict, counter):
    """
    向字典中插入序列
    """
    if seq not in seq_dict:
        seq_dict[seq] = 'p' + str(counter)
        counter += 1
    return seq_dict, counter


if __name__ == '__main__':
    main()

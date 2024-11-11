#!/usr/bin/env python3


import yaml
import click
from pathlib import Path
from Bio import SeqIO


@click.command()
@click.option('-i', '--fasta', type=click.Path(exists=True), help='输入用于设计引物的FASTA序列文件.')
@click.option('-w', '--workdir', help='引物分析工作目录.')
def main(fasta, workdir):
    """
    批量设计 qPCR 引探
    输入多区域 fasta 序列文件, 拆分序列写入到多个配置文件并行跑 primer3, 输出 primer3 结果
    在工作目录中创建 primer3 目录, 每个区域一个配置文件, 每个区域一个输出文件
    """
    records = SeqIO.parse(fasta, 'fasta')
    for record in records:
        name = record.name.replace(':', '_')    # : 替换为 _, 避免文件名冲突
        seq = record.seq
        home = Path(__file__).resolve().parents[1]
        template = home.joinpath('etc/template.p3')
        primer3 = yaml.safe_load(open(home.joinpath('etc/software.yml')).read())['primer3']
        primer3_workdir = Path(workdir).resolve().joinpath('primer3')

        # 写个多个区域配置文件
        with open(template) as f, open(f'{primer3_workdir}/{name}.p3', 'w') as g:
            g.write(f'SEQUENCE_ID={name}\nSEQUENCE_TEMPLATE={seq}\n{f.read()}')


if __name__ == '__main__':
    main()

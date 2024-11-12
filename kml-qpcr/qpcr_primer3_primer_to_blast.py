from pathlib import Path
import sys
import click
import csv
import logging

logging.basicConfig(level=logging.DEBUG,
                    format="%(asctime)s - %(levelname)s - %(filename)s - %(message)s",
                    datefmt="%Y-%m-%d %H:%M:%S")


@click.command()
@click.option('-w', '--workdir', default='qpcr_analysis', help='引物分析工作目录.')
@click.option('-c', '--cores', default=8, type=int, help='使用的 CPU 核心数.')
@click.help_option('-h', '--help')
def main(workdir, cores):
    """
    去除序列重复的引物. 分别导出 fasta 和 table 文件
    """
    primer3_parse_dir = Path(workdir).joinpath('primer3_parse')

    if not primer3_parse_dir.exists():
        logging.error(f'目录 {workdir} 下没有 primer3_parse 文件夹, 请检查.')
        sys.exit(1)


if __name__ == '__main__':
    main()

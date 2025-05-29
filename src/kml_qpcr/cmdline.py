import click

from src.kml_qpcr.gnm_download import download_genome_files
from src.kml_qpcr.gnm_quality_assess import assess_genome_quality
from src.kml_qpcr.cstm_gnms_load import load_customer_genomes
from src.kml_qpcr.gnm_annotate import annotate_genomes


@click.group()
@click.version_option(version="0.1.0", help="显示版本信息.", prog_name="KML-qPCR")
@click.help_option(help="显示帮助信息.")
def cli():
    """病原微生物定量PCR分析工具"""


def common_options(func):
    """通用参数装饰器"""
    func = click.option("--threads", default=4, type=int, show_default=True, help="全局线程数.")(func)
    func = click.option("--genome-set-dir", default="kml_qpcr_genomes",
                        show_default=True, help="项目基因组集目录, 输出结果在该目录下.")(func)
    func = click.option("--sci-name", required=True, help="输入物种学名.")(func)
    func = click.help_option(help="显示帮助信息.")(func)
    return func


@cli.command()
@common_options
def download(sci_name, genome_set_dir, threads):
    """下载参考数据库. 线程参数用于解压, 不适用于下载."""
    download_genome_files(sci_name, genome_set_dir, threads)


@cli.command()
@common_options
@click.option("--customer-genome-dir", required=True, help="输入客户基因组目录.")
def load(customer_genome_dir, sci_name, genome_set_dir, threads):
    """接入客户基因组集, 格式化成符合项目结构的目录结构."""
    load_customer_genomes(customer_genome_dir, sci_name, genome_set_dir, threads)


@cli.command()
@common_options
@click.option("--force", is_flag=True, help="强制重新运行 prokka, 默认识别到结果文件就跳过.")
def annotate(sci_name, genome_set_dir, threads, force):
    """注释基因组"""
    annotate_genomes(sci_name, genome_set_dir, threads, force)


@cli.command()
@common_options
@click.option("--pathogen-type", type=click.Choice(["Bacteria", "Viruses"]), default="Bacteria", show_default=True, help="输入病原类型.")
@click.option("--force", is_flag=True, help="强制重新运行 checkM/checkV, 默认识别到结果文件就跳过.")
def assess(sci_name, genome_set_dir, pathogen_type, threads, force):
    """质控评估"""
    assess_genome_quality(sci_name, genome_set_dir, threads, pathogen_type, force)

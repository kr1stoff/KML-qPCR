import click
from pathlib import Path
from src.kml_qpcr.gnm_download import download_genome_files
from src.kml_qpcr.gnm_quality_assess import assess_genome_quality


@click.group()
@click.version_option(version="0.1.0", help="显示版本信息.", prog_name="KML-qPCR")
@click.help_option(help="显示帮助信息.")
def mycli():
    """病原微生物定量PCR分析工具"""


@mycli.command()
@click.option("--sci-name", required=True, help="输入物种学名.")
@click.option("--genome-set-dir", default="kml_qpcr_genomes", show_default=True, help="项目基因组集目录, 输出结果在该目录下.")
@click.help_option(help="显示帮助信息.")
def download(sci_name, genome_set_dir):
    """下载参考数据库"""
    download_genome_files(sci_name, genome_set_dir)


@mycli.command()
@click.option("--sci-name", required=True, help="输入物种学名.")
@click.option("--genome-set-dir", default="kml_qpcr_genomes", show_default=True, help="项目基因组集目录, 输出结果在该目录下.")
@click.option("--pathogen-type", type=click.Choice(["Bacteria", "Viruses"]), default="Bacteria", show_default=True, help="输入病原类型.")
@click.option("--threads", default=4, type=int, show_default=True, help="全局线程数.")
@click.option("--force", is_flag=True, help="强制重新运行 checkM/checkV, 默认识别到结果文件就跳过.")
@click.help_option(help="显示帮助信息.")
def assess(sci_name, genome_set_dir, pathogen_type, threads, force):
    """质控评估"""
    assess_genome_quality(sci_name, genome_set_dir, threads, pathogen_type, force)

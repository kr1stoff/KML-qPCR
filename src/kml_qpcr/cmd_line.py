import click
from pathlib import Path
from src.kml_qpcr.gnm_download import (
    get_taxonomy_id_from_sciname,
    get_assembly_summary_by_taxids,
    download_genome_files
)
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
    outdir = Path(genome_set_dir).joinpath(sci_name.replace(" ", "_"))
    # 物种基因组集 info 目录
    infodir = Path(outdir).joinpath("info")
    infodir.mkdir(parents=True, exist_ok=True)
    taxids = get_taxonomy_id_from_sciname(sci_name, infodir)
    df_rsgb = get_assembly_summary_by_taxids(taxids, infodir)
    # 下载基因组文件
    alldir = Path(outdir).joinpath("all")
    alldir.mkdir(parents=True, exist_ok=True)
    download_genome_files(df_rsgb, alldir)


@mycli.command()
@click.option("--sci-name", required=True, help="输入物种学名.")
@click.option("--genome-set-dir", default="kml_qpcr_genomes", show_default=True, help="项目基因组集目录, 输出结果在该目录下.")
@click.option("--pathogen-type", type=click.Choice(["Bacteria", "Viruses"]), default="Bacteria", show_default=True, help="输入病原类型.")
@click.option("--threads", default=4, type=int, show_default=True, help="全局线程数.")
@click.help_option(help="显示帮助信息.")
def assess(sci_name, genome_set_dir, pathogen_type, threads):
    """质控评估"""
    gnmdir = Path(genome_set_dir).joinpath(sci_name.replace(" ", "_"))
    assess_genome_quality(gnmdir, threads, pathogen_type)

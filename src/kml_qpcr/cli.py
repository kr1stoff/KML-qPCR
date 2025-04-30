import click
from pathlib import Path
from src.kml_qpcr.download import get_taxonomy_id_from_sciname, get_assembly_summary_by_taxids, download_genome_files


@click.group()
@click.version_option(version="0.1.0", help="显示版本信息.", prog_name="KML-qPCR")
@click.help_option(help="显示帮助信息.")
@click.option("--genome-dir", default="kml_qpcr_genomes", help="项目基因组集目录, 输出结果在.")
@click.pass_context
def mycli(ctx, genome_dir):
    """病原微生物定量PCR分析工具"""
    ctx.ensure_object(dict)
    ctx.obj["gnmdir"] = genome_dir


@mycli.command()
@click.option("--sci-name", required=True, help="输入 scientific name.")
@click.pass_context
def download(ctx, sci_name):
    """下载参考数据库"""
    outdir = Path(ctx.obj["gnmdir"]).joinpath(sci_name.replace(" ", "_"))
    # 物种基因组集 info 目录
    infodir = Path(outdir).joinpath("info")
    infodir.mkdir(parents=True, exist_ok=True)
    taxids = get_taxonomy_id_from_sciname(sci_name, infodir)
    df_rsgb = get_assembly_summary_by_taxids(taxids, infodir)
    # 下载基因组文件
    alldir = Path(outdir).joinpath("all")
    alldir.mkdir(parents=True, exist_ok=True)
    download_genome_files(df_rsgb, alldir)

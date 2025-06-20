import click

from src.kml_qpcr.gnm_download import download_genome_files
from src.kml_qpcr.gnm_quality_assess import GenomeQualityAssessor, GenomeQualityAssessorViruses
from src.kml_qpcr.cstm_gnms_load import load_customer_genomes
from src.kml_qpcr.gnm_annotate import GenomeAnnotator
from src.kml_qpcr.csvd_gene_obtain import ConservedGenePredictor
from src.kml_qpcr.spec_gene_obtain import SpeciticityGeneObtainer


@click.group()
@click.version_option(version="0.1.0", help="显示版本信息.", prog_name="KML-qPCR")
@click.help_option(help="显示帮助信息.")
def cli():
    """病原微生物定量PCR分析工具"""


def common_options(func):
    """通用参数装饰器"""
    func = click.option("--force", is_flag=True, help="强制重新运行第三方软件 默认识别到结果文件就跳过.")(func)
    func = click.option("--threads", default=4, type=int, show_default=True, help="全局线程数.")(func)
    func = click.option("--genome-set-dir", default="kml_qpcr_genomes",
                        show_default=True, help="项目微生物总目录, 输出结果在该目录下.")(func)
    func = click.option("--sci-name", required=True, help="输入物种学名. 例如 'Coxiella Burnetii'.")(func)
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
def annotate(sci_name, genome_set_dir, threads, force):
    """注释基因组"""
    ga = GenomeAnnotator(
        sci_name=sci_name,
        genome_set_dir=genome_set_dir,
        threads=threads,
        force=force
    )
    ga.run()


@cli.command()
@common_options
@click.option("--pathogen-type", type=click.Choice(["Bacteria", "Viruses"]), default="Bacteria", show_default=True, help="输入病原类型.")
def assess(sci_name, genome_set_dir, pathogen_type, threads, force):
    """质控评估"""
    assessor_class = (
        GenomeQualityAssessor
        if pathogen_type == "Bacteria"
        else GenomeQualityAssessorViruses
    )
    gqa = assessor_class(
        sci_name=sci_name,
        genome_set_dir=genome_set_dir,
        threads=threads,
        force=force
    )
    gqa.run()


@cli.command()
@common_options
@click.option("--core-isolates-percent", type=int, default=100, show_default=True, help="核心基因覆盖分离株百分比阈值.")
@click.option("--blastp-identity", type=int, default=100, show_default=True, help="BlastP 相似度阈值")
def conserved(sci_name, genome_set_dir, threads, force, core_isolates_percent, blastp_identity):
    """保守区域预测"""
    cgp = ConservedGenePredictor(
        sci_name=sci_name,
        genome_set_dir=genome_set_dir,
        threads=threads,
        core_islt_perc=core_isolates_percent,
        core_blastp_idnt=blastp_identity,
        force=force
    )
    cgp.run()

@cli.command()
@common_options
def specificity(sci_name, genome_set_dir, threads, force):
    """特异性基因预测"""
    sgo = SpeciticityGeneObtainer(
            sci_name=sci_name,
            genome_set_dir=genome_set_dir,
            threads=threads,
            force=force
    )
    sgo.run()

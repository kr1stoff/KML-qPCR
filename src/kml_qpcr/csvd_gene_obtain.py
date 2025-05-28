from pathlib import Path
from subprocess import run
import logging


def obtain_conserved_genes():
    """"""
    gnmsfile = "/data/mengxf/Project/KML250416_chinacdc_pcr/genomes/Ehrlichia_chaffeensis/genome_assess/checkm/high_quality_genomes.txt"
    alldir = "/data/mengxf/Project/KML250416_chinacdc_pcr/genomes/Ehrlichia_chaffeensis/all"
    gffdir = "/data/mengxf/Project/KML250416_chinacdc_pcr/genomes/Ehrlichia_chaffeensis/conserved_gene/genbank2gff3"
    # 准备 bp_genbank2gff3 的输入
    prepare_bp_genbank2gff3_input(gnmsfile, alldir, gffdir)
    # TODO 后续 roary 流程
    # mamba run -n basic bp_genbank2gff3 --dir conserved_gene/genbank2gff3
    # mamba run -n meta roary -p 32 conserved_gene/genbank2gff3/*gff -f conserved_gene/roary


def prepare_bp_genbank2gff3_input(gnmsfile, alldir, gffdir):
    """
    准备 bp_genbank2gff3 的输入
    :param gnmsfile: 高质量基因组的名称文件
    :param alldir: 所有基因组的目录
    :param gffdir: 输出目录
    :return:
    """
    high_qlt_gnms = open(gnmsfile, "r").read().strip().split("\n")
    for gnm in high_qlt_gnms:
        try:
            gbff = list(Path(alldir).glob(f"{gnm}/*.gbff.gz"))[0]
        except IndexError:
            # logging.error(f"Cannot find gbff file for {gnm} in {alldir}")
            logging.error(f"基因组 {alldir}/{gnm} 没有找到对应的 gbff 文件")
            continue
        run(f"cp {gbff} {gffdir}/{gnm}.gbff.gz", shell=True, check=True)

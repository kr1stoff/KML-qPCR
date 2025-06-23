import logging
from pathlib import Path
from subprocess import run
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from src.kml_qpcr.base import BaseQPCR
from src.config.cnfg_database import BLAST_CORE_NT
from src.config.cnfg_software import BLASTN


class SpeciticityGeneObtainer(BaseQPCR):
    def __init__(self, sci_name: str, genome_set_dir: str, threads: int, force: bool):
        """
        获取特异性基因
        :param sci_name: 物种名称
        :param genome_set_dir: 基因组集合目录
        :param threads: 线程数
        :param force: 是否强制重新计算
        """
        super().__init__(sci_name, genome_set_dir, threads, force)
        self.csvd_gene_dir = self.gnm_dir / "conserved_gene"
        self.blast_dir = self.gnm_dir / "specific_gene/blast"
        self.blast_dir.mkdir(exist_ok=True, parents=True)
        self.blast_query = self.blast_dir / "csvd_gene_seq_set.fasta"
        self.blast_out = self.blast_dir / "blastn.tsv"

    def run(self):
        self.merge_csvd_gene_for_blast()
        # 是否强制重新运行
        if not self.force and self.blast_out.exists():
            logging.warning("BLAST 结果文件已存在, 跳过")
            return
        self.run_blast()

    def merge_csvd_gene_for_blast(self):
        """合并保守基因形成 fasta, 用于 blast 比对"""
        # * 多次运行没有删除 csvd_gene_seq_set 目录里面预测的基因可能是乱的. 读预测到的保守基因列表, 再做合并
        with open(self.csvd_gene_dir / "core_single_copy_genes.txt") as f:
            csvd_genes = [line.strip() for line in f]
        # 读取每个保守基因的 fasta 文件, 并合并成一个 fasta 文件
        genes = []
        for gene in csvd_genes:
            ffn = f"{self.csvd_gene_dir}/csvd_gene_seq_set/{gene}.ffn"
            parser = SeqIO.parse(ffn, "fasta")
            first_gene = next(parser)
            seqrcd = SeqRecord(first_gene.seq, id=gene, description="")
            genes.append(seqrcd)
        SeqIO.write(genes, self.blast_query, "fasta")

    def run_blast(self):
        """运行 blast 命令"""
        logging.info("运行 blast core nt")
        blastdb = Path(BLAST_CORE_NT).resolve().parent
        cmd = fr"""
        export BLASTDB={blastdb}
        {BLASTN} -num_threads {self.threads} \
            -query {self.blast_query} \
            -db {BLAST_CORE_NT} \
            -out {self.blast_out} \
            -perc_identity 60 \
            -qcov_hsp_perc 60 \
            -outfmt "6 qseqid sseqid ssciname staxid pident qcovs length nident evalue bitscore"
        sed -i '1i\qseqid\tsseqid\tssciname\tstaxid\tpident\tqcovs\tlength\tnident\tevalue\tbitscore' {self.blast_out}
        """
        logging.debug(cmd)
        run(cmd, shell=True, check=True, executable="/bin/bash")

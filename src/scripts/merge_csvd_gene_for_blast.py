from pathlib import Path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


gnm_dir = "/data/mengxf/Project/KML250416_chinacdc_pcr/genomes/Coxiella_Burnetii"
gnm_dir = Path(gnm_dir)
csvd_gene_dir = gnm_dir / "conserved_gene/csvd_gene_seq_set"
blast_dir = gnm_dir / "specific_gene/blast"
blast_dir.mkdir(exist_ok=True, parents=True)


genes = []
for ffn in csvd_gene_dir.glob("*.ffn"):
    gene = ffn.stem
    parser = SeqIO.parse(ffn, "fasta")
    first_gene = next(parser)
    seqrcd = SeqRecord(first_gene.seq, id=gene, description="")
    genes.append(seqrcd)

SeqIO.write(genes, blast_dir / "csvd_gene_seq_set.fasta", "fasta")

# blastout
"""
mamba run -n basic blastn -num_threads 16 \
    -query csvd_gene_seq_set.fasta \
    -db /data/mengxf/Database/NCBI/blast/db/core_nt/core_nt \
    -out blastn.out2 \
    -perc_identity 60 \
    -qcov_hsp_perc 60 \
    -outfmt "6 qseqid sseqid staxid pident qcovs length qlen slen sstart send qstart qend nident evalue bitscore"
"""

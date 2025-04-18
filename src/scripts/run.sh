# 窗口300 步长30 不限行宽度
seqkit sliding -s 30 -W 300 -w 0 HPV23genomes.fa -o genomes_primer/sliding.fa

# primer3 批量引物设计
poetry run python kml-qpcr/qpcr_primer3_design.py -c 32 \
    -i /data/mengxf/Project/KML241108_HPV_qPCR_primer_design/genomes_primer/sliding.fa \
    -w /data/mengxf/Project/KML241108_HPV_qPCR_primer_design/genomes_primer

# 解析 primer3 引物到 表格
poetry run python kml-qpcr/qpcr_primer3_parse.py \
    -w /data/mengxf/Project/KML241108_HPV_qPCR_primer_design/genomes_primer -c 32

# blast 输入
poetry run python kml-qpcr/qpcr_primer3_to_blast.py \
    -w /data/mengxf/Project/KML241108_HPV_qPCR_primer_design/genomes_primer

# blast
time /home/mengxf/miniforge3/envs/basic/bin/blastn -task blastn-short \
    -word_size 7 -num_threads 32 -perc_identity 70 -qcov_hsp_perc 50 -max_hsps 100 -max_target_seqs 100 \
    -query /data/mengxf/Project/KML241108_HPV_qPCR_primer_design/genomes_primer/blast/primers.fasta \
    -db /data/mengxf/Project/KML241108_HPV_qPCR_primer_design/HPV23genomes.fa \
    -out /data/mengxf/Project/KML241108_HPV_qPCR_primer_design/genomes_primer/blast/primers.out \
    -outfmt "6 qseqid sseqid length qlen slen sstart send qstart qend qcovs pident nident evalue bitscore sstrand"

# 包容性评估

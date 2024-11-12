# 窗口300 步长30 不限行宽度
seqkit sliding -s 30 -W 300 -w 0 HPV23genomes.fa -o genomes_primer/sliding.fa

# primer3 批量引物设计
poetry run python kml-qpcr/qpcr_primer3_design.py -c 32 -i /data/mengxf/Project/KML241108_HPV_qPCR_primer_design/genomes_primer/sliding.fa -w /data/mengxf/Project/KML241108_HPV_qPCR_primer_design/genomes_primer

# 解析 primer3 引物到 表格
poetry run python kml-qpcr/qpcr_primer3_parse.py -w /data/mengxf/Project/KML241108_HPV_qPCR_primer_design/genomes_primer/ -c 32

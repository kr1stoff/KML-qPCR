# 列出病原列表的 taxid, rank. assembly_summary 没有物种名这一列, 但是有 species_taxid 这一列
# csvtk cut -f NCBI命名 pathogen_list.csv | sed '1d' | mamba -n basic run taxonkit name2taxid --data-dir /data/mengxf/Database/NCBI/taxonomy --show-rank --sci-name | csvtk -t pretty -H
echo "Borrelia burgdorferi" | mamba -n basic run taxonkit name2taxid --data-dir /data/mengxf/Database/NCBI/taxonomy --show-rank --sci-name
# species group rank - "spotted fever group"
mamba -n basic run taxonkit list --data-dir /data/mengxf/Database/NCBI/taxonomy --indent '' --ids 114277
# species rank - Ehrlichia chaffeensis
mamba -n basic run taxonkit list --data-dir /data/mengxf/Database/NCBI/taxonomy --indent '' --ids 945

# 出血热 S 基因引物设计
cd /data/mengxf/Project/KML250416_chinacdc_pcr/genomes/Orthonairovirus_haemorrhagiae
# S 基因
grep '>' nucleotide/sequence.fasta | grep -Ev 'UNVERIFIED|MAG' | grep complete | grep -E 'segment S|nucleocapsid protein' | sed 's/>//g' >nucleotide/S.headers.txt
seqkit grep --threads 32 --line-width 0 --by-name -f nucleotide/S.headers.txt -o nucleotide/S.fna nucleotide/sequence.fasta
mkdir -p genome_assess/checkv
mamba -n qpcr run checkv end_to_end -t 8 -d /data/mengxf/Database/checkV/checkv-db-v1.5 nucleotide/S.fna genome_assess/checkv/S
cp genome_assess/checkv/S/quality_summary.tsv genome_assess/checkv/checkv_summary.tsv
csvtk csv2xlsx -t genome_assess/checkv/checkv_summary.tsv
# * 根据 checkv 做质量筛选
csvtk -t filter2 -f '$checkv_quality!="Not-determined"' genome_assess/checkv/checkv_summary.tsv | csvtk -t cut -f contig_id | sed '1d' >genome_assess/checkv/S.remain_list.txt
seqkit grep --threads 32 --line-width 0 --by-name -r -f genome_assess/checkv/S.remain_list.txt nucleotide/S.fna -o genome_assess/checkv/S.filtered.fna
# 滑窗获取分段序列
seqkit sliding --step 100 --window 200 --line-width 0 genome_assess/checkv/S.filtered.fna -o primer_design/S.sliding.fna
# 去重
seqkit rmdup --by-seq --line-width 0 primer_design/S.sliding.fna -o primer_design/S.sliding.rmdup.fna
# primer3 引物设计
# ! 输入输出用绝对路径
mamba -n python3.12 run poetry -C /data/mengxf/GitHub/KML-qPCR run python -m src.kml_qpcr.primer3_design --threads 32 \
    --fasta /data/mengxf/Project/KML250416_chinacdc_pcr/genomes/Orthonairovirus_haemorrhagiae/primer_design/S.sliding.rmdup.fna \
    --workdir /data/mengxf/Project/KML250416_chinacdc_pcr/genomes/Orthonairovirus_haemorrhagiae/primer_design
#  primer3 引物解析
mamba -n python3.12 run poetry -C /data/mengxf/GitHub/KML-qPCR run python -m src.kml_qpcr.primer3_parse --threads 32 \
    --workdir /data/mengxf/Project/KML250416_chinacdc_pcr/genomes/Orthonairovirus_haemorrhagiae/primer_design
#3 计算包容性和输出简并引物
mamba run -n python3.12 poetry -C /data/mengxf/GitHub/KML-qPCR run python -m src.kml_qpcr.inclusivity_and_degenerate \
	--ref-seqs /data/mengxf/Project/KML250416_chinacdc_pcr/genomes/Orthonairovirus_haemorrhagiae/genome_assess/checkv/S.filtered.fna \
	--workdir /data/mengxf/Project/KML250416_chinacdc_pcr/genomes/Orthonairovirus_haemorrhagiae/primer_design \
	--threads 32

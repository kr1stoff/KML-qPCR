# 列出病原列表的 taxid, rank. assembly_summary 没有物种名这一列, 但是有 species_taxid 这一列
# csvtk cut -f NCBI命名 pathogen_list.csv | sed '1d' | mamba -n basic run taxonkit name2taxid --data-dir /data/mengxf/Database/NCBI/taxonomy --show-rank --sci-name | csvtk -t pretty -H
echo "Borrelia burgdorferi" | mamba -n basic run taxonkit name2taxid --data-dir /data/mengxf/Database/NCBI/taxonomy --show-rank --sci-name
# species group rank - "spotted fever group"
mamba -n basic run taxonkit list --data-dir /data/mengxf/Database/NCBI/taxonomy --indent '' --ids 114277
# species rank - Ehrlichia chaffeensis
mamba -n basic run taxonkit list --data-dir /data/mengxf/Database/NCBI/taxonomy --indent '' --ids 945
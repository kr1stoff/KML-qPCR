from pathlib import Path
from src.kml_qpcr.gnm_quality_assess import get_genome_anno_quality, merge_checkm_rna_stats, filter_by_merge_df_bacteria

# gnmdir = "/data/mengxf/Project/KML250416_chinacdc_pcr/genomes/Borreliella_burgdorferi"
# gnmdir = "/data/mengxf/Project/KML250416_chinacdc_pcr/genomes/Coxiella_Burnetii"
# get_genome_anno_quality(Path(gnmdir))


assdir = "/data/mengxf/Project/KML250416_chinacdc_pcr/genomes/Coxiella_Burnetii/genome_assess"
# merge_checkm_rna_stats(Path(assdir))
filter_by_merge_df_bacteria(Path(assdir))

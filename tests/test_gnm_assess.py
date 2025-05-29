from pathlib import Path
from src.kml_qpcr.gnm_quality_assess import get_genome_anno_quality

gnmdir = "/data/mengxf/Project/KML250416_chinacdc_pcr/genomes/Borreliella_burgdorferi"
# gnmdir = "/data/mengxf/Project/KML250416_chinacdc_pcr/genomes/Coxiella_Burnetii"

get_genome_anno_quality(Path(gnmdir))

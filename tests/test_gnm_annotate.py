from src.kml_qpcr.gnm_annotate import prokka_is_complete_run
from pathlib import Path

anno_dir = "/data/mengxf/Project/KML250416_chinacdc_pcr/genomes/Coxiella_Burnetii/genome_annotate"
# anno_dir = "/data/mengxf/Project/KML250416_chinacdc_pcr/genomes/Borreliella_burgdorferi/genome_annotate"

print(prokka_is_complete_run(Path(anno_dir)))

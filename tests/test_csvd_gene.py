from pathlib import Path
from src.kml_qpcr.csvd_gene_obtain import ConservedGenePredictor


gnm_dir = "/data/mengxf/Project/KML250416_chinacdc_pcr/genomes/Coxiella_Burnetii"

cgp = ConservedGenePredictor(
    "Coxiella Burnetii",
    "/data/mengxf/Project/KML250416_chinacdc_pcr/genomes",
    32, 100, 100, False)
# core_sglcp_genes = cgp.filter_core_single_copy_gene()
# cgp.seqkit_split_isolates_ffn()
# cgp.output_conserved_gene_set(core_sglcp_genes)
cgp.run()

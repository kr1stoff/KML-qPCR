from src.kml_qpcr.spec_gene_obtain import SpeciticityGeneObtainer


# sgo = SpeciticityGeneObtainer(
#     "Coxiella Burnetii",
#     "/data/mengxf/Project/KML250416-chinacdc-pcr/genomes",
#     32, False
# )
sgo = SpeciticityGeneObtainer(
    "Ehrlichia chaffeensis",
    "/data/mengxf/Project/KML250416-chinacdc-pcr/genomes",
    32, False
)
sgo.merge_csvd_gene_for_blast()

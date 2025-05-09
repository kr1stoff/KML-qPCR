from src.kml_qpcr.gnm_download import get_taxonomy_id_from_sciname


# sciname = "Escherichia coli"
# sciname = "spotted fever group"
sciname = "Ehrlichia chaffeensis"

taxids = get_taxonomy_id_from_sciname(sciname)
print('\n'.join(taxids))

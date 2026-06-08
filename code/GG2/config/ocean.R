list(
  dataset = "ocean",
  input = list(
    format = "mitag",
    table = "./ocean/data/downloaded_data/16S/miTAG.taxonomic.profiles.release.tsv",
    sequences = "./ocean/data/downloaded_data/16S/16S.OTU.SILVA.reference.sequences.fna"
  ),
  output = list(
    # Base directory where GG2 working folders (input/, PI/) will be created
    directory = "./ocean/data/processed_data/16S/GG2"
  ),
  greengenes = list(
    backbone_qza = "./greengenes/2024.09.backbone.full-length.fna.qza",
    tree_qza     = "./greengenes/2024.09.phylogeny.asv.nwk.qza",
    taxonomy_qza = "./greengenes/2024.09.taxonomy.id.nwk.qza"
  ),
  vsearch = list(
    default_perc_identity = 0.99
  ),
  feature_id_column = "OTU.rep", # Column name for feature IDs in the input table
  sample_mapping = list(
    # For TARA: map sample_label -> analysis_id for traceability to environmental data
    lookup_file         = "./ocean/data/processed_data/environmental/filtered/sample_id_lookup.csv",
    sample_label_column = "sample_label",
    analysis_id_column  = "analysis_id"
  )
)


list(
  dataset = "soil",
  input = list(
    format = "dada2",
    table = "./soil/data/processed_data/16S/DADA/merged/merged_seqtab_asv.csv",
    sequences = "./soil/data/processed_data/16S/DADA/merged/merged_asv_sequences.fna"
  ),
  output = list(
    # Base directory where GG2 working folders (input/, PI/) will be created
    directory = "./soil/data/processed_data/16S/GG2"
  ),
  greengenes = list(
    backbone_qza = "./greengenes/2024.09.backbone.full-length.fna.qza",
    tree_qza     = "./greengenes/2024.09.phylogeny.asv.nwk.qza",
    taxonomy_qza = "./greengenes/2024.09.taxonomy.id.nwk.qza"
  ),
  vsearch = list(
    default_perc_identity = 0.99
  ),
  feature_id_column = "ASV_ID", # Name used in the original DADA2 table (columns after sample_id)
  sample_mapping = NULL # Soil sample_id already matches environmental data
)


# CCA config for ocean.
# Env: envdata.csv from read_data (data_path/environmental/filtered/envdata.csv or env_file).
# Composition: GG2 exports under data_path/16S/GG2/<perc_identity>/final.
# Paths are relative to repo root.

list(
  dataset = "ocean",
  data_path = "./ocean/data/processed_data",
  results_path = "./ocean/results/CCA",
  env_file = NULL,
  hyperparam = list(
    k = 15,
    lambda1_range = 10^seq(-2.5, 5, length.out = 21),
    lambda2_range = 0,
    seed = 1
  ),
  null = list(
    n_shuffles = 101,
    seeds = 0:100
  ),
  plotting_params = list(
    env_var_labels = c( # no units because measurements are scaled anyway so the loadings aren't in those units
      depth = "depth",
      temperature = "temperature",
      salinity = "salinity",
      oxygen = "oxygen",
      no2 = "NO2",
      po4 = "PO4",
      hco3 = "HCO3",
      co3 = "CO3"
    ),
    tax_palette <- c(
      Phylum  = "#decfb2",
      Class   = "#ddc68e",
      Order   = "#dcbd6c",
      Family  = "#dab449",
      Genus   = "#bd863e",
      Species = "#a05532",
      OTU     = "#812727"
    )
  )
)

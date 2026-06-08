# CCA config for soil.
# Env: envdata.csv from read_data (data_path/envdata.csv or env_file).
# Composition: GG2 exports under data_path/16S/GG2/<perc_identity>/final.
# Paths are relative to repo root.

list(
  dataset = "soil",
  data_path = "./soil/data/processed_data",
  results_path = "./soil/results/CCA",
  # Env: default is data_path/environmental/filtered/envdata.csv; override if envdata.csv lives elsewhere
  env_file = NULL,
  # Composition: will look in data_path/16S/GG2/<perc_identity>/final for seqtab.tsv, taxonomy.csv
  hyperparam = list(
    k = 15,
    lambda1_range = 10^seq(-5, 5, length.out = 21),
    lambda2_range = 0,
    seed = 1
  ),
  null = list(
    n_shuffles = 101,
    seeds = 0:100
  ),
  # Optional: plotting (env var order/labels for env loadings in step 4 and 5; tax_palette for step 5 crosstax compare)
  # env_var_labels: named vector (names = var names in data, matching read_data base_vars; values = display labels)
  plotting_params = list(
    env_var_labels = c(
      ph = "pH",
      soil_c = "C",
      soil_n = "N",
      soil_p = "P",
      soil_c_n_ratio = "C:N",
      clay_silt = "clay-silt"
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

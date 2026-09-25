# CCA config for soil restricted to one hemisphere (latitude from site_locations.csv).
# Same settings as soil.R; only `hemisphere` and `results_path` differ.
# Paths are relative to repo root.

cfg <- source("code/CCA/config/soil.R", local = TRUE)$value
cfg$hemisphere <- "south"
cfg$results_path <- "./soil/results/CCA_hemisphere/south"
cfg

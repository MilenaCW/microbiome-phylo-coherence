# soil (from env_metadata.xlsx)
# - geographic: latitude, longitude
# - measured: clay_silt, soil_c, soil_n, soil_p, soil_c_n_ratio, ph
# - EO/meteorological: aridity_index, mdr, max_t, min_t, psea, uv_light, npp
# - classifiers: continent, forest, grassland, shrubland, ecosystem_type
# Base variables to keep (ordered): ph, soil_c, soil_n, soil_c_n_ratio, soil_p, clay_silt

list(
  dataset = "soil",
  input = list(
    # NOTE: this file is not present in the repo snapshot; update if needed.
    path = "soil/data/downloaded_data/env_metadata.xlsx"
  ),
  output = list(
    dir = "soil/data/processed_data/environmental/filtered"
  ),
  # Base variables to keep (order relevant)
  base_vars = c("ph",
    "soil_c",
    "soil_n",
    "soil_p",
    "clay_silt"),
  # Optional extras (kept unless missingness/redundancy filters drop them)
  extra_vars = c(),
  rename_map = c(
    "ID_sequencing" = "sample_id",
    "Latitude" = "latitude",
    "Longitude" = "longitude",
    "pH" = "ph",
    "Clay_silt" = "clay_silt",
    "Soil_C" = "soil_c",
    "Soil_N" = "soil_n",
    "Soil_P" = "soil_p",
    "Soil_C_N_ratio" = "soil_c_n_ratio",
    "Aridity_Index" = "aridity_index",
    "MDR" = "mdr",
    "MAXT" = "max_t",
    "MINT" = "min_t",
    "PSEA" = "psea",
    "UV_Light" = "uv_light",
    "NPP2003_2015" = "npp",
    "Continent" = "continent",
    "Forest" = "forest",
    "Grassland" = "grassland",
    "Shrubland" = "shrubland",
    "Ecosystem_type" = "ecosystem_type"
  ),
  filters = list(
    # Optional numeric ranges for outlier trimming (set empty list() for none)
    ranges = list( # as estimated by eye
      ph = list(keep_na = FALSE),
      soil_c = list(max = 15, keep_na = FALSE),
      soil_n = list(max = 0.6, keep_na = FALSE),
      soil_p = list(max = 1500, keep_na = FALSE),
      clay_silt = list(keep_na = FALSE)
    )
  ),
  var_catalog = tibble::tribble(
    ~dataset, ~var, ~type, ~units, ~description, ~source_column,
    "soil", "latitude", "geographic", "deg", "Latitude", "Latitude",
    "soil", "longitude", "geographic", "deg", "Longitude", "Longitude",
    "soil", "clay_silt", "measured", "percent", "Texture (% clay + silt)", "Clay_silt",
    "soil", "soil_c", "measured", "percent", "Soil carbon", "Soil_C",
    "soil", "soil_n", "measured", "percent", "Soil nitrogen", "Soil_N",
    "soil", "soil_p", "measured", "mg P/kg soil", "Soil phosphorus", "Soil_P",
    "soil", "soil_c_n_ratio", "measured", "unitless", "Soil C:N ratio", "Soil_C_N_ratio",
    "soil", "ph", "measured", "unitless", "Soil pH", "pH",
    "soil", "aridity_index", "eo_meteorological", "Precipitation / Evapotranspiration", "Aridity index", "Aridity_Index",
    "soil", "mdr", "eo_meteorological", "degC", "Mean diurnal temperature range", "MDR",
    "soil", "max_t", "eo_meteorological", "degC", "Max temperature", "MAXT",
    "soil", "min_t", "eo_meteorological", "degC", "Min temperature", "MINT",
    "soil", "psea", "eo_meteorological", "CV, %", "Precipitation seasonality", "PSEA",
    "soil", "uv_light", "eo_meteorological", "unitless", "UV intensity (0-16)", "UV_Light",
    "soil", "npp", "eo_meteorological", "NDVI index", "Net primary productivity proxy", "NPP2003_2015",
    "soil", "continent", "classifier", NA, "Continent", "Continent",
    "soil", "forest", "classifier", NA, "Forest indicator", "Forest",
    "soil", "grassland", "classifier", NA, "Grassland indicator", "Grassland",
    "soil", "shrubland", "classifier", NA, "Shrubland indicator", "Shrubland",
    "soil", "ecosystem_type", "classifier", NA, "Ecosystem type", "Ecosystem_type"
  )
)

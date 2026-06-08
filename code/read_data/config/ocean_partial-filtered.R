# ocean (from OM.CompanionTables.xlsx Table W1/W8 + miTAG table header)
# - joins: sample_label -> analysis_id (from miTAG headers) and sample_label -> pangea_id (W1)
# - measured (core): depth, temperature, po4, no2, no2no3, si, salinity, oxygen
# - geographic: latitude, longitude (from W8 Mean_Lat/Mean_Long)
# - optional: carbonate chemistry from TARA_SAMPLES_CONTEXT_ENV-DEPTH-CARB (median/q2 only)
# Base variables to keep (ordered): depth, po4, si, salinity, oxygen
# Extra variables (removed for analysis): temperature, n02, n02n03; carbonate vars in extra_vars

list(
  dataset = "ocean",
  input = list(
    mitag_tax_profiles_tsv = "ocean/data/downloaded_data/16S/miTAG.taxonomic.profiles.release.tsv",
    companion_tables_xlsx = "ocean/data/downloaded_data/environmental/OM.CompanionTables.xlsx",
    carbonate_xlsx = "ocean/data/downloaded_data/environmental/TARA_SAMPLES_CONTEXT_ENV-DEPTH-CARB_20170515.xlsx"
  ),
  carbonate = list(
    target_measurements = c(
      "pH", "Carbon dioxide", "Carbon dioxide, partial pressure",
      "Fugacity of carbon dioxide in seawater", "Bicarbonate ion", "Carbonate ion",
      "Carbon, total", "Alkalinity, total", "Calcite saturation state", "Aragonite saturation state"
    )
  ),
  output = list(
    dir = "ocean/data/processed_data/environmental/partial-filtered"
  ),
  # Base variables to keep (order relevant)
  base_vars = c("depth",
                "temperature",
                "po4",
                "no2no3",
                "si",
                "no2",
                "salinity",
                "oxygen"),
  # Optional extras (kept unless missingness/redundancy filters drop them)
  extra_vars = c(),
  # Rename map applied to Table W8 and to carbonate data (raw col names -> canonical)
  rename_map = c(
    "PANGAEA Sample ID" = "pangea_id",
    "Mean_Lat*" = "latitude",
    "Mean_Long*" = "longitude",
    "Mean_Depth [m]*" = "depth",
    "Mean_Temperature [deg C]*" = "temperature",
    "Mean_Salinity [PSU]*" = "salinity",
    "Mean_Oxygen [umol/kg]*" = "oxygen",
    "NO2 [umol/L]**" = "no2",
    "PO4 [umol/L]**" = "po4",
    "NO2NO3 [umol/L]**" = "no2no3",
    "SI [umol/L]**" = "si",
    # Carbonate (CARB XLSX): raw = lowercased measurement + _q2
    "ph_q2" = "pH",
    "carbon dioxide_q2" = "co2",
    "carbon dioxide, partial pressure_q2" = "co2_pp",
    "fugacity of carbon dioxide in seawater_q2" = "co2_f",
    "bicarbonate ion_q2" = "hco3",
    "carbonate ion_q2" = "co3",
    "carbon, total_q2" = "carbon_total",
    "alkalinity, total_q2" = "alkalinity_total",
    "calcite saturation state_q2" = "calcite_saturation_state",
    "aragonite saturation state_q2" = "aragonite_saturation_state"
  ),
  filters = list(
    # Optional numeric ranges for outlier trimming (set empty list() for none)
    ranges = list( # as estimated by eye
      depth = list(keep_na = FALSE),
      temperature = list(keep_na = FALSE),
      po4 = list(keep_na = FALSE),
      no2 = list(keep_na = FALSE),
      no2no3 = list(keep_na = FALSE),
      si = list(keep_na = FALSE),
      salinity = list(keep_na = FALSE),
      oxygen = list(keep_na = FALSE)
    )
  ),
  var_catalog = tibble::tribble(
    ~dataset, ~var, ~type, ~units, ~description, ~source_column, ~table_source,
    "ocean", "latitude", "geographic", "deg", "Mean latitude", "Mean_Lat*", "Table W8",
    "ocean", "longitude", "geographic", "deg", "Mean longitude", "Mean_Long*", "Table W8",
    "ocean", "depth", "measured", "m", "Mean depth", "Mean_Depth [m]*", "Table W8",
    "ocean", "temperature", "measured", "degC", "Mean temperature", "Mean_Temperature [deg C]*", "Table W8",
    "ocean", "salinity", "measured", "PSU", "Mean salinity", "Mean_Salinity [PSU]*", "Table W8",
    "ocean", "oxygen", "measured", "umol/kg", "Mean oxygen", "Mean_Oxygen [umol/kg]*", "Table W8",
    "ocean", "no2", "measured", "umol/L", "NO2", "NO2 [umol/L]**", "Table W8",
    "ocean", "po4", "measured", "umol/L", "PO4", "PO4 [umol/L]**", "Table W8",
    "ocean", "no2no3", "measured", "umol/L", "NO2 + NO3", "NO2NO3 [umol/L]**", "Table W8",
    "ocean", "si", "measured", "umol/L", "Silica", "SI [umol/L]**", "Table W8",
    "ocean", "pH", "calculated", "", "pH median; calculated quantity (median) from measured values", "pH (median)", "CARB",
    "ocean", "co2", "calculated", "", "Carbon dioxide median; calculated quantity (median) from measured values", "Carbon dioxide (median)", "CARB",
    "ocean", "co2_pp", "calculated", "", "CO2 partial pressure median; calculated quantity (median) from measured values", "Carbon dioxide, partial pressure (median)", "CARB",
    "ocean", "co2_f", "calculated", "", "Fugacity CO2 seawater median; calculated quantity (median) from measured values", "Fugacity of carbon dioxide in seawater (median)", "CARB",
    "ocean", "hco3", "calculated", "", "Bicarbonate ion median; calculated quantity (median) from measured values", "Bicarbonate ion (median)", "CARB",
    "ocean", "co3", "calculated", "", "Carbonate ion median; calculated quantity (median) from measured values", "Carbonate ion (median)", "CARB",
    "ocean", "carbon_total", "calculated", "", "Total carbon median; calculated quantity (median) from measured values", "Carbon, total (median)", "CARB",
    "ocean", "alkalinity_total", "calculated", "", "Total alkalinity median; calculated quantity (median) from measured values", "Alkalinity, total (median)", "CARB",
    "ocean", "calcite_saturation_state", "calculated", "", "Calcite saturation state median; calculated quantity (median) from measured values", "Calcite saturation state (median)", "CARB",
    "ocean", "aragonite_saturation_state", "calculated", "", "Aragonite saturation state median; calculated quantity (median) from measured values", "Aragonite saturation state (median)", "CARB"
  )
)

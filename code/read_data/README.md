# Environmental data reader (code/read_data)

Loads raw environmental metadata from source files, standardizes column names, applies outlier filters, and writes `envdata.csv` (and supporting files) for use by the CCA pipeline.

---

## Scripts

| Script | Purpose |
|---|---|
| `read_envdata.R` | Central reader: loads source files, renames columns, filters outliers, writes outputs |
| `env_diagnostic.R` | Generates missingness UpSet plots and correlation correlograms (PDFs) |
| `read_all.sh` | Runs reader + diagnostics for soil and ocean with a given mode (`filtered` or `full`) |

---

## Configuration

Config files live in `config/` and must evaluate to an R `list`. Pre-built configs cover the datasets used in the manuscript:

| Config | Dataset | Variables |
|---|---|---|
| `soil_filtered.R` | Soil | pH, soil C/N/P, clay+silt (outlier-filtered) |
| `soil_full.R` | Soil | All available soil environmental variables |
| `ocean_filtered.R` | Ocean | Depth, temperature, PO4, NO2, salinity, oxygen |
| `ocean_partial-filtered.R` | Ocean | Full core variable set before high-correlation variables were removed |
| `ocean_full.R` | Ocean | All available ocean environmental variables |

### Key config fields

```r
list(
  dataset = "soil",   # "soil" or "ocean"

  input = list(
    # Dataset-specific paths to raw source files (Excel / XLSX / TSV).
    # Soil: path = "soil/data/downloaded_data/env_metadata.xlsx"
    # Ocean: mitag_tax_profiles_tsv, companion_tables_xlsx (+ optionally carbonate_xlsx)
  ),

  output = list(
    dir = "soil/data/processed_data/environmental/filtered"
  ),

  base_vars = c("ph", "soil_c", "soil_n", "soil_p", "clay_silt"),  # ordered; all written to envdata.csv
  extra_vars = c(),                                                   # optional; appended after base_vars

  rename_map = c(
    "pH"       = "ph",         # source column name -> canonical name
    "Soil_C"   = "soil_c",
    # ...
  ),

  filters = list(
    ranges = list(
      ph     = list(keep_na = FALSE),                  # drop NA only, no numeric bound
      soil_c = list(max = 15, keep_na = FALSE),        # drop rows where soil_c > 15 or NA
      soil_p = list(min = 0, max = 1500, keep_na = FALSE)
      # ...
    )
  ),

  var_catalog = tibble::tribble(
    ~dataset,  ~var,      ~type,       ~units,    ~description,    ~source_column,
    "soil",    "ph",      "measured",  "unitless", "Soil pH",       "pH",
    # ...
  )
)
```

`filters$ranges` entries support `min`, `max`, and `keep_na` (logical). Rows that fall outside the range or are NA (when `keep_na = FALSE`) are dropped; a summary is printed to the console.

---

## Running

Run reader + diagnostics for soil and ocean with the `filtered` or `full` variable set:

```bash
bash code/read_data/read_all.sh filtered
bash code/read_data/read_all.sh full
```

Run a single config manually (required for `ocean_partial-filtered.R` and any custom config):

```bash
Rscript code/read_data/read_envdata.R  --config code/read_data/config/ocean_partial-filtered.R
Rscript code/read_data/env_diagnostic.R --config code/read_data/config/ocean_partial-filtered.R
```

---

## Outputs

Written to the `output$dir` specified in the config:

| File | Contents |
|---|---|
| `envdata.csv` | sample_id + selected environmental variables; rows = samples |
| `var_catalog.csv` | Variable documentation: dataset, var, type, units, description, source column |
| `site_locations.csv` | sample_id, latitude, longitude (for mapping) |
| `sample_id_lookup.csv` | Ocean only — maps original sample labels to analysis IDs |
| `env_upset_missingness.pdf` | UpSet plot: which variables are missing together |
| `env_correlogram.pdf` | Pairwise correlation matrix with density diagonals |

---

## Integration with CCA

CCA step 0 (`code/CCA/scripts/00_read_data.R`) reads `envdata.csv` from the path specified in the CCA config. Run this module for each dataset/config combination before running the CCA pipeline.

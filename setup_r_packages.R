# Run once after creating the project environment:
#   Rscript setup_r_packages.R
# Installs R packages not available via conda.

install.packages("phylolm", repos = "https://cloud.r-project.org")

if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes", repos = "https://cloud.r-project.org")
}
remotes::install_github("ElenaTuzhilina/RCCA")

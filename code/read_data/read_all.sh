# Usage: bash read_all.sh <filtered|full>
MODE=${1:?Usage: $0 <filtered|full>}

# soil
Rscript read_envdata.R \
 --config ./code/read_data/config/soil_${MODE}.R
Rscript env_diagnostic.R \
 --config ./code/read_data/config/soil_${MODE}.R

# ocean
Rscript read_envdata.R \
 --config ./code/read_data/config/ocean_${MODE}.R
Rscript env_diagnostic.R \
 --config ./code/read_data/config/ocean_${MODE}.R

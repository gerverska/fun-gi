# Update the conda environment and R package space according to state files ####

# conda ####
conda env update -p ./env -f conda.yml

# R ####
R -e 'renv::restore(prompt = F)'

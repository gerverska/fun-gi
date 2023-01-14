#!/bin/bash

# This script allows for the creation of environments on different platforms ####
# This script will install the specified versions of each requested program ####

# Setup ####
conda remove -p ./env --all
> .Rprofile
rm renv -fr

# Configure conda before environment creation ####
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority flexible

# Create a conda environment from the exported .yml ####
conda env create -f conda.yml -p ./env

# Download all R packages specified in renv.lock ####
conda run -p ./env Rscript code/00-build.r

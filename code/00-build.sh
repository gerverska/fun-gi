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

# Create an environment and install packages from the bioconda channel ####
conda create -y -p ./env -c bioconda fastqc=0.11.9 multiqc=1.13 \
pheniqs=2.1.0 cutadapt=4.2 itsx=1.1.3 vsearch=2.22.1 

# Install packages from the conda-forge channel ####
conda install -y -p ./env -c conda-forge r-base=4.2.2 r-renv=0.16.0 \
pkg-config=0.29.2 libgit2=1.5.0 git=2.39.0 gh=2.21.2

# Download all R packages specified in renv.lock ####
conda run -p ./env Rscript code/00-build.r

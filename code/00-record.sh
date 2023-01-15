#!/bin/bash

# Record the current states of the conda environment and R package space ####

# conda ####
conda env export --from-history > conda.yml

# R ####
R -e 'renv::snapshot()'

#!/bin/bash

# This script allows for the creation of environments on different platforms ####

# Setup ####
conda remove -p ./env --all

# Create a conda environment from the exported .yml ####
conda env create -c -f env.yml -p ./env

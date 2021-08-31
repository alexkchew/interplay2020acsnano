#!/bin/bash

# generate_cluster_of_points.sh
# The purpose of this script is to generate clusters of groups. 
# Then, using these groups, we will compute the hydrophobicity. 

# Written by: Alex K. Chew (12/8/2019)

## VARIABLES:


################################
### LOADING GLOBAL VARIABLES ###
################################
dir_loc="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
source "${dir_loc}/bin/bashrc.sh"

#########################
### INITIAL VARIABLES ###
#########################

## DEFINING WORKING DIRECTORY
path_sims="/home/shared/np_hydrophobicity_project/simulations/NP_SPHERICAL"
# "$1"

## DEFINING SPECIFIC PATH
specific_sim="MostlikelynpNVT_EAM_300.00_K_2_nmDIAM_C11COOH_CHARMM36jul2017_Trial_1_likelyindex_1"

## DEFINING ANALYSIS DIRECTORY
analysis_dir="${3:-hyd_analysis}"

## DEFINING GRID DIRECTORY
grid_dir="${4:-grid-0_1000-0.1}"

## DEFINING GRID FILE
grid_file="${5:-out_willard_chandler.dat}"

#########################
### DEFAULT VARIABLES ###
#########################
num_cluster_within_group="5"
num_cores="20"

## PYTHON SCRIPT
python_script="${PYTHON_SCRIPTS}/cluster_grid_points.py"

##########################################
### MAIN SCRIPT 
##########################################

## DEFINING FULL PAHT
full_path_sim="${path_sims}/${specific_sim}"

#############################
### RUNNING PYTHON SCRIPT ###
#############################
python3.6 "${python_script}" --path "${full_path_sim}" \
                             --analysis "${analysis_dir}" \
                             --gridfolder "${grid_dir}" \
                             --griddat "${grid_file}" \
                             --nwithingroup "${num_cluster_within_group}" \
                             --num_cores "${num_cores}"
                             






#!/bin/bash

# full_run_np_analysis.sh
# The purpose of this script is to run full nanoparticle analysis
# 
# Written by: Alex K. Chew & Brad C. Dallin (10/20/2019)

################################
### LOADING GLOBAL VARIABLES ###
################################
dir_loc="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
source "${dir_loc}/bin/bashrc.sh"

#########################
### DEFAULT VARIABLES ###
#########################

## BASH FILE
bash_file="${PATH2SCRIPTS}/run_np_analysis.sh"

## PREFIX
folder_prefix="MostlikelynpNVT_EAM_300.00_K_2_nmDIAM"
folder_suffix="CHARMM36jul2017_Trial_1_likelyindex_1"

############################
### DEFINING INPUT ARRAY ###
############################

## MAIN SIMULATION FOLDER
main_sim_folder="NP_SPHERICAL"
# "190920-Most_likely_np_sims_NVT_mixed_solvent_ROT_ligs"

## DEFINING LIGAND ARRAY
declare -a ligand_names=("C11OH")
# "ROT001" "ROT002"

## LOOPING
for each_ligand in "${ligand_names[@]}"; do
	## CREATING SIM FOLDER
	current_sim_folder="${folder_prefix}_${each_ligand}_${folder_suffix}"

	## PRINTING
	echo "Running: bash "${bash_file}" "${main_sim_folder}" "${current_sim_folder}""; sleep 2

	## RUNNING BASH SCRIPT
	bash "${bash_file}" "${main_sim_folder}" "${current_sim_folder}"

done




#!/bin/bash

# prep_pymol_water_snapshot.sh
# The purpose of this script is to prepare a snapshot of the system with water.
#
# Written by: Alex K. Chew (04/08/2020)
#
# USAGE:
#    bash /home/shared/np_hydrophobicity_project/run_scripts/prep_pymol_scripts/prep_pymol_water_snapshot.sh

################################
### LOADING GLOBAL VARIABLES ###
################################
dir_loc="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
source "${dir_loc}/../bin/bashrc.sh"

## INPUTS

## DEFINING PATH TO SIMULATION
path_to_sim="/home/akchew/scratch/nanoparticle_project/simulations/20200618-GNP_COSOLVENT_MAPPING/comap_PRO_1_1_50-EAM_300.00_K_2_nmDIAM_C11OH_CHARMM36jul2017_Trial_1_likelyindex_1"

## DEFINING INPUT PREFIX
input_prefix="sam_prod"

## DEFINING OUTPUT
output_prefix="${input_prefix}-pymol_snapshot"

## DEFINING LAST FRAME
last_frame="12000"

## CHECK IF EXISTS
stop_if_does_not_exist "${path_to_sim}"

## GOING TO DIR
cd "${path_to_sim}"

## RUNNING TRJCONV
gmx trjconv -f "${input_prefix}.xtc" \
            -s "${input_prefix}.tpr" \
            -o "${output_prefix}.pdb" \
            -dump "${last_frame}" \
            -pbc mol << INPUTS
System
INPUTS

         
echo "Completed generation for: ${path_to_sim}"
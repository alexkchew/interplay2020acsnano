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


#### FOR GNP SIMS

## SIM DIR NAME
sim_dirname="20200618-Most_likely_GNP_water_sims_FINAL"
# "20200401-Renewed_GNP_sims_with_equil"
## SPECIFIC FILE
specific_sim="MostlikelynpNVTspr_50-EAM_300.00_K_2_nmDIAM_dodecanethiol_CHARMM36jul2017_Trial_1_likelyindex_1"
specific_sim="MostlikelynpNVTspr_50-EAM_300.00_K_2_nmDIAM_C11OH_CHARMM36jul2017_Trial_1_likelyindex_1"

#### FOR PLANAR SIMS

## SIM DIR NAME
# sim_dirname="20200403-Planar_SAMs-5nmbox_vac_with_npt_equil"
# ## SPECIFIC FILE
# specific_sim="NVTspr_50_Planar_300.00_K_dodecanethiol_10x10_CHARMM36jul2017_intffGold_Trial_1-5000_ps"
# 

## DEFINING PATH TO SIM
path_to_sim="${PATH_SIMULATIONS}/${sim_dirname}/${specific_sim}"

## DEFINING INPUT PREFIX
input_prefix="sam_prod"

## DEFINING OUTPUT
output_prefix="${input_prefix}-pymol_water_snapshot"

## DEFINING LAST FRAME
# last_frame="50000"

## CHECK IF EXISTS
stop_if_does_not_exist "${path_to_sim}"

## GOING TO DIR
cd "${path_to_sim}"

## RUNNING TRJCONV
gmx trjconv -f "${input_prefix}.gro" \
            -s "${input_prefix}.tpr" \
            -o "${output_prefix}.pdb" \
            -pbc mol << INPUTS
System
INPUTS

#             -dump "${last_frame}" \


echo "Completed generation for: ${path_to_sim}"
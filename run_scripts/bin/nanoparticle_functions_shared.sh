#!/bin/bash

# nanoparticle_functions_shared.sh
# This is shared functions between nanoparticle project
# and lig from dir
#
# FUNCTIONS:
#   - extract_lig_name_from_dir: extract lig name from directory

## FUNCTION TO EXTRACT LIGAND NAME FROM NOMENCLATURE
# The purpose of this function is to extract ligand name from 
# nomenclature 'EAM_300.00_K_2_nmDIAM_ROT005_CHARMM36jul2017_Trial_1'
# INPUTS:
#   $1: directory name
# OUTPUTS:
#   ligand name
# USAGE:
#   ligand_name=$(extract_lig_name 'EAM_300.00_K_2_nmDIAM_ROT005_CHARMM36jul2017_Trial_1')
# DEBUGGING:
#   cut -d'_' -f1 <<< 'switch_solvents-50000-dmso-EAM_300.00_K_2_nmDIAM_ROT011_CHARMM36jul2017_Trial_1'
function extract_lig_name_from_dir () {
    ## DEFINING INPUTS
    input_dir_name_="$1"
    
    ## DEFINING INITIAL TEXT
    initial_text=$(cut -d'_' -f"1" <<< "${input_dir_name_}")
    
    if [[ "${initial_text}" == "EAM" ]]; then
        ligand_index_="6" # 6th position based on splitting '_'
        # EAM_300.00_K_2_nmDIAM_ROT011_CHARMM36jul2017_Trial_1
    elif [[ "${initial_text}" == "switch" ]] || [[ "${initial_text}" == "Mostlikelynp"* ]]; then
        ligand_index_="7" # 6th position based on splitting '_'
        # switch_solvents-50000-dmso-EAM_300.00_K_2_nmDIAM_ROT011_CHARMM36jul2017_Trial_1
    elif [[ "${initial_text}" == "MostlikelynpNVTspr_" ]]; then
        ligand_index_="8"

    elif [[ "${initial_text}" == "FrozenPlanar"* ]] || [[ "${initial_text}" == "FrozenGoldPlanar"* ]]; then
        ligand_index_="4"
    elif [[ "${initial_text}" == "NVTspr"* ]]; then
        ligand_index_="6"
    else
        echo "Error! Ligand name is not well-defined"
        echo "Check extract_lig_name_from_dir function in nanoparticle_functions.sh"
        sleep 5
    fi
    ## FINDING LIGANDN AME
    ligand_name=$(cut -d'_' -f"${ligand_index_}" <<< "${input_dir_name_}")
    
    ## PRINTING
    echo "${ligand_name}"

}
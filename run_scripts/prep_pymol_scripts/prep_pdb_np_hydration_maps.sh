#!/bin/bash

# prep_pdb_np_hydration_maps.sh
# This script is designed to prepare nanoparticle PDB files for hydration maps.

# VARIABLES:
#   $1: simulation path
#   $2: job name:
#       - none: default
#       - planar: planar SIMs
################################
### LOADING GLOBAL VARIABLES ###
################################
dir_loc="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
source "${dir_loc}/../bin/bashrc.sh"

#########################
### INITIAL VARIABLES ###
#########################

## DEFINING FORCEFIELD
forcefield="${FORCEFIELD}"

## DEFINING PATH TO FOLDERS
folder_name="$1"
path_to_folders="${PATH_SIMULATIONS}/${folder_name}"

## DEFINING PRODUCTION TPR AND XTC
input_tpr="sam_prod.tpr"
input_xtc="sam_prod.xtc"

## DEFINING FRAME TO OUTPUT
frame="2000"

## DEFINING OUTPUT PDB
output_pdb="sam_prod_${frame}.pdb"

## DEFINING IF YOU WANT SYSTEM
want_system=true

#########################
### DEFAULT VARIABLES ###
#########################

## DEFINING INDEX
index_file="gold_ligand.ndx"

###################
### MAIN SCRIPT ###
###################

## LOOPING THROUGH EACH FOLDER
for path_sim in ${path_to_folders}/*/ ; do

    ## CHECKING IF PATH EXISTS
    stop_if_does_not_exist "${path_sim}"

    ## GOING INTO DIRECTORY
    cd "${path_sim}"

    ## GETTING BASENAME
    sim_basename=$(basename ${path_sim})
    ## CHECKING THE DIFFERENT NAMES
    if [[ "${sim_basename}" != "FrozenGoldPlanar"* ]] && [[ "${sim_basename}" != "FrozenPlanar"* ]] && [[ "${sim_basename}" != "NVTspr"* ]]; then
        echo "Running nanoparticle case extraction"
        ## EXTRACTING LIGAND NAME
        read -a extract_array <<< $(extract_output_NVT_hydration_maps ${sim_basename})
        lig_name=${extract_array[1]}  
        
        ## DEFINING GOLD RESIDUE NAME
        gold_residue_name="AUNP"

    else
        echo "Running planar case extraction"
        ## PLANAR CASE
        lig_name=$(extract_lig_name_from_dir ${sim_basename})
        
        ## DEFINING GOLD RESIDUE NAME
        gold_residue_name="AUI"
        
    fi

    ## DEFINING LIGAND NAME ARRAY
    read -a ligand_name_array <<< "$(str2array_by_delim "${lig_name}" ",")"

    ## STORING LIG RES NAME
    declare -a lig_res_name_array=()

    ## LOOPING THROUGH EACH LIGAND
    for lig_name in ${ligand_name_array[@]}; do

        ## DEFINING ITP FILE
        lig_itp_file="${lig_name}.itp"

        ## DEFINING ITP FILE
        if [[ "${sim_basename}" == "FrozenGoldPlanar"* ]] || [[ "${sim_basename}" == "FrozenPlanar"* ]] || [[ "${sim_basename}" == "NVTspr"* ]]; then
            lig_itp_file="${forcefield}/${lig_itp_file}"
        fi

        ## CHECKING IF PATH EXISTS
        stop_if_does_not_exist "${lig_itp_file}"

        ## FINDING LIGAND RESIDUE NAME
        lig_residue_name="$(itp_get_resname ${lig_itp_file})"

        ## STORING LIG RESIDUE NAME
        lig_res_name_array+=( "${lig_residue_name}" )

    done

    residue_name_spaces="$(join_array_to_string " " ${lig_res_name_array[@]})"
    residue_name_underscore="$(join_array_to_string "_" "${lig_res_name_array[@]}")"

    ## CREATING INDEX FILE
    make_ndx_gold_with_ligand "${path_sim}/${input_tpr}" \
                              "${path_sim}/${index_file}"   \
                              "${residue_name_spaces}"   \
                              "${gold_residue_name}"  

    ## DEFINING COMBINED NAME
    combined_name="${gold_residue_name}_${residue_name_underscore}"

## GMX TRJCONV
gmx trjconv -s ${input_tpr} -f "${input_xtc}" -o "${output_pdb}" -dump "${frame}" -n ${index_file} -pbc mol << INPUTS
${combined_name}
INPUTS

if [[ "${want_system}" == true ]]; then
    output_pdb_system="${output_pdb%.pdb}_system.pdb"
    echo "Since want_system is true, outputting PDB with system: ${output_pdb_system}"

## GMX TRJCONV
gmx trjconv -s ${input_tpr} -f "${input_xtc}" -o "${output_pdb_system}" -dump "${frame}" -n ${index_file} -pbc mol << INPUTS
System
INPUTS

fi

done
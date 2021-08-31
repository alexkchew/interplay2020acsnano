#!/bin/bash

# functions.sh
# This contains all functions
## FUNCTIONS:
#   python3.6: runs python 3.6. Use pip3.6 to install new modules here.
#   extract_output_NVT_hydration_maps: function to extract input names for NVT

## ALIASING PYTHON3.6
function python3.6 () {
    "/usr/bin/python3" "$@"
}

### FUNCTION TO EXTRACT NOMENCLATURE
# This extracts outputname for NVT hydration maps
# INPUTS:
#   $1: input name
#       e.g. MostlikelynpNVT_EAM_300.00_K_2_nmDIAM_C11CONH2_CHARMM36jul2017_Trial_1_likelyindex_1
# OUTPUTS:
#   Array with the following:
#       0: diameter
#       1: ligand name
#       2: likely index
# USAGE:
#   read -a extract_array <<< $(extract_output_NVT_hydration_maps ${name})
function extract_output_NVT_hydration_maps (){
    ## DEFINING INPUTS
    input_name_="$1"
    
    ## SPLITTING ARRAY
    my_array=($(echo ${input_name_} | tr "_" "\n")) #"_"
    # RETURNS: 8_nm 300_K 1_mf aceticacid_formate_methylammonium_propane
    
    ## RETURNING
    diameter_=${my_array[4]}
    lig_name_=${my_array[6]}
    likely_index_=${my_array[-1]}
    
    ## DECLARING OUTPUT ARRAY
    declare -a output_array=("${diameter_}" \
                             "${lig_name_}" \
                             "${likely_index_}"
                             )
    
    ## PRINTING
    echo "${output_array[@]}"

}
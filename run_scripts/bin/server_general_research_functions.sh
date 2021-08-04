#!/bin/bash
# server_general_research_functions.sh
# This script contains general code that can be sourced by any script.
# It's main purpose is to mainstream transferrable functions that can be used in multiple projects

### FUNCTION LIST ###
### SPECIALIZED FUNCTIONS
#       combined_alchemical_molecule_name: generates a combined alchemical molecule name

### GENERALIZED FUNCTIONS
#       print_script_name: prints name of the script
#       join_by: joins arrays and generates a string
#       str2array_by_delim: converts string to array by using a deliminator

## DIRECTORY FUNCTIONS
#       check_file_exist: checks if a file exists (True/False output)
#       create_dir: checks and creates directories
#       stop_if_does_not_exist: function to check if file exists and stops if the file does not exist

## indexing functions
#       index_read_list_: reads index file names
#       index_check_: checks to see if residue name is within index file
#       intersection_two_arrays: finds intersection between two arrays
#       index_make_ndx_residue_: creates index file for a specific residue
#       generate_split_array: can split array up

## GENERAL FUNCTIONS
#       compute_num_cores_required: computes the number of cores that are required

## TRAJECTORY FUNCTIONS
#       xtc_check_total_time_ps: checks for total frames in an xtc file (same as checkxtc)

## GRO FILE FUNCTIONS
#       read_gro_box_lengths: read gro box lengths
#       compute_volume_from_box_length: computs box length from reading gro file
#       center_gro_: center gro file
#       gro_measure_box_size: reads gro box size (same as previous)
#       gro_compute_volume: computes volume based on box size
#       gro_compute_num_density: computes number density based on volume
#       gro_count_total_atoms: counts total atoms for a gro file
#       gro_combine_two_files: combines two gro files, assuming temp file 2 should have the box size
#       gro_get_first_residue_name: get first residue name
#       truncate_select_nearby_gro_: select nearby items for g ro -- truncation
#       gro_expand_box: expands the box
#       gro_combine_two_files_by_z: combines gro file by z-dimension

## MDP FILE FUNCTIONS
#       mdp_find_option_values: finds mdp options
#       mdp_add_energy_groups: adds energy group to mdp file

## ITP FUNCTIONS
#       itp_fix_genrestr_numbering: fix itp files restraints
#       convert_itp_nocharge: takes a itp file, then turns off all charges on it
#       itp_get_resname: gets itp file residue name
#       itp_extract_prm: extracts parameter file from itp
#       itp_get_residue_name_within_atom: get residue name within atom directive
#       itp_find_total_atoms: finds total number of atoms based on the itp file

## TOPOLOGY FUNCTIONS
#       top_clear_molecules: clears topology file of molecules
#       top_read_molecules: reads topology molecules
#       add_include_topology: Adds an include statement at the end of a topology file
#       topology_add_include_specific: Includes specific text into topology based on findings
#       topology_add_b4_include_specific: Includes topology information before a specific string
#       top_update_residue_num: corrects number of residues for a topology file
#       top_correct_residue_name: corrects topology based on a gro file
#       top_extract_total_residue: extracts total residue from topology file
#       top_extract_multiple_residues: method to extract multiple residue values to a csv file
#       fix_posre_for_single_residue: fixes position restraints for a single resiude.

## MOLECULAR INFO FUNCTIONS
#       check_molecule_info: check if the molecule information exists
#       write_molecular_info: function that writes molecular information
#       calc_pdb_molecular_weight: computes pdb file molecular weight
#       check_xtc_time: faster version of checking xtc time

## MISCELLANEOUS
#       round_nearest_integer: rounds value to nearest integer
#       join_array_to_string: joins arrays to a string
#       extract_variable_name_from_file: looks for variable name in a file, then extracts it

## PRINTING
echo "*** LOADING GENERAL FUNCTIONS (server_general_research_functions.sh) ***"

## ALIASING PYTHON3
function python3 () {
    "/usr/bin/python3.4" "$@"
}


## ALIASING PYTHON3.6
function python3.6 () {
    "/usr/bin/python3" "$@"
}
## SERVER FUNCTIONS
# Changing work folder for different servers


### FUNCTION TO COMBINE MOLECULE NAME
# The purpose of this function is combine the molecule names for an alchemical-analysis
# INPUTS:
#   $1: molecule_1_name
#   $2: molecule_2_name
# OUTPUTS:
#   combined molecule name
# USAGE:
#   combined_molecule_name "12-propanediol" "propanal"
#       This will output: 12-propanediol~propanal
function combined_alchemical_molecule_name () {
    ## DEFINING INPUTS
    molecule_1_=$1
    molecule_2_=$2
    
    ## DEFINING OUTPUTS
    output_molecule_="${molecule_1_}~${molecule_2_}"
    
    ## PRINTING
    echo "${output_molecule_}"
}

### FUNCTION TO CONVERT STRING TO ARRAY
# The purpose of this function is to convert a string to an array
# INPUTS:
#   $1: input string
#   $2: deliminator type (e.g. ';')
# OUTPUTS:
#   array that you can use
# USAGE EXAMPLE:
    #   str2array "mixed_solvent_methanol_60-EAM_300.00_K_2_nmDIAM_C11CF3_CHARMM36jul2017_Trial_1"
#   As a variable:
#       read -a test_var <<< "$(str2array "mixed_solvent_methanol_60-EAM_300.00_K_2_nmDIAM_C11CF3_CHARMM36jul2017_Trial_1")"
function str2array_by_delim () {
    ## INPUTS
    input_string="$1"
    deliminator="${2-_}" # type of deliminator
    
    ## READING ARRAY
    IFS="${deliminator}" read -ra converted_array <<< "${input_string}"
    
    ## PRINTING ARRAY
    echo "${converted_array[@]}"
}


### FUNCTION TO FIND HOST NAME INFORMATION
# The purpose of this function is to find the hostname information (e.g. key word, wall time, etc.)
# INPUTS:
#   - simply the variables you want to store it in
# OUTPUTS:
#   - this function rewrites the variable names and outputs the correct host name information
## USAGE: find_hostname_info SERVERKEYWORD WALLTIME NUMOFCPUS
#   e.g. echo "
function find_hostname_info () {
    ## DEFINING INPUTS
    local __serverkeyword="${1-SERVERKEYWORD}"
    local __walltime="${2:-WALLTIME}"
    local __ncpus="${3:-NUMOFCPUS}"

if [[ $(hostname) == *"stampede"* ]]; then # Check if we are in the stampede server
    SERVERKEYWORD="STAMPEDE" # Server name (Used to find the folder with submission scripts)
    WALLTIME="48:00:00" # Change WALLTIME in days-hh:mm:ss
    NUMOFCPUS=32 # Change NUMBEROFCPUS

elif [[ $(hostname) == *"comet"* ]]; then
    SERVERKEYWORD="COMET" # Server name (Used to find the folder with submission scripts)
    WALLTIME="48:00:00" # Change WALLTIME in days-hh:mm:ss
    NUMOFCPUS=24 # Change NUMBEROFCPUS

elif [[ $(hostname) == *"swarm"* ]]; then
    SERVERKEYWORD="SWARM" # Server name (Used to find the folder with submission scripts)
    WALLTIME="7-00:00:00" # Change WALLTIME in days-hh:mm:ss
    NUMOFCPUS=28 # Change NUMBEROFCPUS

elif [[ $(echo $HOSTNAME) == *"chtc"* ]]; then # Check if we are at the CHTC Supercomputers
    SERVERKEYWORD="CHTC_HPC" # Server name (Used to find the folder with submission scripts)
    WALLTIME="5-00:00:00" # Change WALLTIME in days-hh:mm:ss
    NUMOFCPUS=20 # Change NUMBEROFCPUS
else
    SERVERKEYWORD=""
    WALLTIME=""
    NUMOFCPUS=""
fi
    ## EVALUATING
    eval ${__serverkeyword}="${SERVERKEYWORD}"
    eval ${__walltime}="${WALLTIME}"
    eval ${__ncpus}="${NUMOFCPUS}"
}

## FUNCTION TO PRINT SCRIPT NAME (No inputs)
# USAGE: print_script_name
function print_script_name () {
    script_name=`basename "$0"`
    echo "*** RUNNING SCRIPT: ${script_name} ***"
}

## FUNCTION TO JOIN ARRAYS
# USAGE: join_by , "${FOO[@]}" #a,b,c
# REFERENCE: https://stackoverflow.com/questions/1527049/how-can-i-join-elements-of-an-array-in-bash
function join_by { local IFS="$1"; shift; echo "$*"; }

### FUNCTION TO READ GRO FILE BOX
# The purpose of this function is to read the last line of the gro file and get the  box lengths
## INPUTS:
#   $1: gro file
## OUTPUTS:
#   box_length in the form of an array
## USAGE:
#   read -a test<<<$(read_gro_box_lengths sphere_gold.gro)
#   echo "${test[0]} <-- should get you the first box length, etc.
function read_gro_box_lengths () {
    # DEFINING VARIABLES
    gro_file="$1"
    # FINDING BOX LENGTH
    declare -a box_length=$(tail ${gro_file} -n1 | awk '{print $1, $2, $3}')
    ## PRINTING
    echo "${box_length[@]}"
}

### FUNCTION TO COMPUTE BOX VOLUME
# The purpose of this function is to compute the box volume given the box length
## INPUTS:
#   $1: box length 1
#   $2: box length 2
#   $3: box length 3
## OUTPUTS:
#   box_length in the form of an array
## USAGE:
#   molecular_volume="$(compute_volume_from_box_length "${box_length[@]}")"
function compute_volume_from_box_length () {
    ## DEFINING VARIABLES
    box_length_1="$1"
    box_length_2="$2"
    box_length_3="$3"
    
    ## COMPUTING BOX VOLUME
    box_volume=$(echo "${box_length_1}*${box_length_2}*${box_length_3}" | bc -l )

    echo "${box_volume}"
}



### FUNCTION TO CALCULATE BOX VOLUME
# The purpose of this function is to get the volume of a gro file based on the box lengths
## INPUTS:
#   $1: gro file
## OUTPUTS:
#   box volume (float)
## USAGE:
#   GENERAL USAGE: -- simply prints out the volume
#       read_gro_box_volumes sphere_gold.gro
#   VARIABLE USAGE: 
#       volume=$(read_gro_box_volumes sphere_gold.gro)
function read_gro_box_volumes () {
    # DEFINING VARIABLES
    gro_file="$1"
    ## FINDING THE BOX LENGTHS
    read -a box_length <<< $(read_gro_box_lengths ${gro_file})
    ## FINDING BOX VOLUME
    box_volume=$(awk -v x_dim=${box_length[0]} -v y_dim=${box_length[1]} -v z_dim=${box_length[2]} 'BEGIN{ printf "%.3f", x_dim * y_dim * z_dim }')
    ## PRINTING
    echo "${box_volume}"
}

### FUNCTION TO SELECT SPECIFIC GROUP AND OUTPUT A NEW GRO FILE
#   The purpose of this function is to select specific gro based on groups that are nearby, then creates a new gro file based on that.
#   Assumptions:
#       - You have a AUNP system with ligands
#       - You need to find all molecules around this system
#   INPUTS:
#       $1: input tpr file
#       $2: input gro file
#       $3: output gro file
#       $4: residue names (e.g. AUNP, or 'AUNP RO1')
#       $5: truncate distance in nm
#       $6: index file
#   OUTPUTS:
#       new output gro file that has the selected gro information
function truncate_select_nearby_gro_ () {
    ## DEFINING INPUTS
    input_tpr_file_="$1"
    input_gro_file_="$2"
    output_gro_file_="$3"
    residue_names_="$4"
    truncate_distance_="$5"
    index_file_name="${6:-truncate.ndx}" # Default name: truncate.ndx
    
    ## PRINTING
    echo -e "\n--------- truncate_select_nearby_gro_ ---------"
    echo "  Selecting the following residues: ${residue_names_}"
    echo "  Truncation distance: ${truncate_distance_} nm"
    echo "  Working on gro file: ${input_gro_file_} --> ${output_gro_file_}"
    echo "  TPR file: ${input_tpr_file_} // index file: ${index_file_name}"
    
    ## RUNNING GMX SELECT
    gmx select -f "${input_gro_file_}" -s "${input_tpr_file_}" -on "${index_file_name}" -select "same residue as within ${truncate_distance_} of resname ${residue_names_} or resname ${residue_names_}" >/dev/null 2>&1
    
    ## EXPORTING VIA TRJCONV
    gmx trjconv -f "${input_gro_file_}" -n "${index_file_name}" -o ${output_gro_file_} -s "${input_tpr_file_}" >/dev/null 2>&1
}

### FUNCTION TO CENTER MOLECULE
#   The purpose of this function is to center the box for a specific molecule
#   Assumptions:
#       - You will be inputting a gro file
#       - You may have a case where multiple of the same residue shows up in your TPR file
#   INPUTS:
#       $1: input tpr file
#       $2: input gro file
#       $3: output gro file
#       $4: centering residue name
#       $5: index file
#   OUTPUTS:
#       output gro file with molecule centered
function center_gro_ () {
    ## DEFINING INPUTS
    input_tpr_file_="$1"
    input_gro_file_="$2"
    output_gro_file_="$3"
    center_residue_name_="$4"
    index_file_name="${5:-center.ndx}" # Default name: index.ndx
    
    ## PRINTING
    echo -e "\n--------- center_gro_ "---------
    echo "CREATING INDEX FILE FOR: ${center_residue_name_}"
    ## CREATING AN INDEX FILE
gmx make_ndx -f "${input_tpr_file_}" -o "${index_file_name}" >/dev/null 2>&1 << INPUTS
keep 0
r ${center_residue_name_}
q
INPUTS
# keep 1
    echo "CENTERING GRO FILE: ${input_gro_file_} --> ${output_gro_file_}"
    echo "   TPR FILE: ${input_tpr_file_}"
    echo "   INDEX FILE: ${index_file_name}"
    ## CENTERING THE MOLECULE
gmx trjconv -f "${input_gro_file_}" -s "${input_tpr_file_}" -o "${output_gro_file_}" -n "${index_file_name}" -center -pbc mol >/dev/null 2>&1 << INPUTS
${center_residue_name_}
System
INPUTS
}



### DEFINING FUNCTION TO READ ALL INDEX FILE VARIABLE NAMES
# The purpose of this function is to load the index file and read all the variables starting with "[ ]". This script is useful in selection of indexes that you care about.
## INPUTS:
#   $1: index file that you want to read
## OUTPUTS:
#   list of the index file variable names
## USAGE: read -a index_list <<< $(index_read_list_ ${index_file})
function index_read_list_ () 
{
    ## DEFINING INPUTS
    index_file_="$1"
    ## READING INDEX LIST, STORING AS AN ARRAY
    declare -a index_list=($(grep -E "\[" ${index_file_}  | sed "s/\[ //g" | sed "s/ \]//g"))
    ## STORING RESULTS
    echo "${index_list[@]}"
}

### FUNCTION TO SEE IF SOLVENT IS WITHIN THE INDEX
# INPUTS:
#   $1: residue name
#   $2: index file
function index_check_ () {
    l_res_name="$1"
    l_index_file="$2"
    if [ -z "$(grep ${l_res_name} ${l_index_file})" ]; then 
        echo "False"
    else
        echo "True"
    fi
}

### FUNCTION TO CREATE A SPECIFIC INDEX FOR A RESIDUE
## The purpose of this function is to make an index for a specific residue. This function is useful if your residue is multiply defined --- requiring you to create an index that is correct.
# INPUTS:
#       $1: input tpr file
#       $2: output ndx file
#       $3: atom selection (residue name)
# OUTPUTS:
#       index file with only the residue
# USAGE: index_make_ndx_residue_ test.tpr index.ndx SOL
function index_make_ndx_residue_ () {
    ## DEFINING INPUTS:
    input_tpr_file_="$1"
    ndx_file_="$2"    
    atom_selection_="$3"
    
    ## PRINTING
    echo -e "\n----- index_make_ndx_residue_ -----"
    echo "  Creating index file for residue: ${atom_selection}"
    echo "  Using TPR file: ${input_tpr_file_} // index file: ${ndx_file_}"
    
### CREATING SPECIFIC INDEX FILE
gmx make_ndx -f "${input_tpr_file_}" -o "${ndx_file_}" &> /dev/null << INPUTS
keep 0
keep 1
r ${atom_selection_}
q
INPUTS
}

### FUNCTION TO CREATE A SPECIFIC INDEX FOR A RESIDUE WITH THE SYSTEM
## The purpose of this function is to make an index for a specific residue. This function is useful if your residue is multiply defined --- requiring you to create an index that is correct.
# INPUTS:
#       $1: input tpr file
#       $2: output ndx file
#       $3: atom selection (residue name)
# OUTPUTS:
#       index file with only the residue
# USAGE: index_make_ndx_residue_ test.tpr index.ndx SOL
function index_make_ndx_residue_with_system () {
    ## DEFINING INPUTS:
    input_tpr_file_="$1"
    ndx_file_="$2"    
    atom_selection_="$3"
    
    ## PRINTING
    echo -e "\n----- index_make_ndx_residue_ -----"
    echo "  Creating index file for residue: ${atom_selection}"
    echo "  Using TPR file: ${input_tpr_file_} // index file: ${ndx_file_}"
    
### CREATING SPECIFIC INDEX FILE
gmx make_ndx -f "${input_tpr_file_}" -o "${ndx_file_}" &> /dev/null << INPUTS
keep 0
r ${atom_selection_}
q
INPUTS
}



## FUNCTION TO CHECK IF THE FILE EXISTS
# This function simply checks if a file exists and outputs true/false
# $1: Full path to file
# USAGE: check_file_exist $FILE
function check_file_exist () {
    if [ -e "$1" ]
    then
        echo "True"
    else
        echo "False"
    fi
}

## FUNCTION TO CHECK IF FILE EXISTS AND STOP IF IT DOES NOT
# The purpose of this function is to check if a file exists. If it does not, then stop the script!
#   $1: full path toe file
# USAGE: stop_if_does_not_exist $FILE
function stop_if_does_not_exist () {
    ## PRINTING
    echo "Checking if file exists: $1 ..."
    if [ -e "$1" ]
    then
        echo "--> $1 exists, continuing!"
    else
        echo "XXX $1 does not exist!"
        echo "Check for errors! You may have a non-existant, vital file!"
        echo "Pausing here for 5 seconds so you can see this error! ..."
        sleep 5
        exit 1
    fi
}

## FUNCTION TO CHECK IF DIRECTORY EXISTS
# The purpose of this function is to check if a directory exists. If not, then stop the script!
## INPUTS:
# $1: full path to the directory
## OUTPUTS:
#   This function will output a print function. It will cancel the script if the directory does in fact not exist
function check_dir_exist () {
    ## DEFINING INPUT
    input_dir_path_="$1"
    echo -e "\n~~~ check_dir_exists ~~~"
    ## CHECKING IF DIRECTORY EXISTS
    if [ -d "${input_dir_path_}" ]; then
        echo "${input_dir_path_} exists ... continuing script!"; sleep 2
    else
        echo "Error!!! ${input_dir_path_} does not exist!"
        echo "Pausing here so you can see this error!"
        echo "Stopping the script"; sleep 5
        exit 1
    fi
}

### FUNCTION TO CHECK IF DIRECTORY EXISTS
# This function checks if the directory exists. If so, prompt removal
# $1: Directory you are interested in
# $2: -f <-- forced, create create directory without prompt
# USAGE: create_dir $DIRECTORY -f
function create_dir () {
    # DEFINING VARIABLES
    directory="$1"
    # CHECK IF DIRECTORY EXISTS
    dir_exists=$(check_file_exist ${directory})
    if [ "${dir_exists}" == "True" ]; then
        if [ "$2" != "-f" ]; then
            echo "${directory} already exists! Do you want to delete and recreate? (y/n)"
            read deletion_criteria
            while [ "${deletion_criteria}" != "y" ] && [ "${deletion_criteria}" != "n" ]; do
                echo "Error! Incorrect prompt! Deletion criteria can only be \"y\" or \"n\", try again:"
                read deletion_criteria
            done
        else # Forcing is on!
            deletion_criteria="y"
        fi
        
        if [ "${deletion_criteria}" == "y" ]; then
            echo "Deleting and recreating ${directory}, pausing 3 seconds..."
            sleep 3
            rm -rv "${directory}"
            mkdir -p "${directory}"
        elif [ "${deletion_criteria}" == "n" ]; then
            echo "Stopping here -- deletion prevented"
            echo "Check if you need these files. This error message is to prevent overwriting of data"
            exit
        fi
    else
        echo "Creating ${directory}"
        mkdir -p "${directory}"
    fi
}

### FUNCTION TO ADD TO TOPOLOGY AFTER LAST INCLUDE
# $1: TOPOLOGY FILE
# $2: STRING TO INCLUDE
# USAGE: add_include_topology gold_ligand.top butanethiol.itp
function add_include_topology () {
    ## DEFINING VARIABLES
    topology_file="$1"
    string_to_include="$2"
    ## GETTING LINE NUMBER OF LAST INCLUDE
    line_num_last=$(grep -nE "\#include" ${topology_file} | sed 's/\([0-9]*\).*/\1/' |tail -n1)
    ## ADDING TO THE TOPOLOGY FILE
    sed -i "$(( line_num_last+1 ))i#include \"$2\"" "$1"
}

### FUNCTION TO ADD TO TOPOLOGY AFTER SPECIFIC INCLUDE
# $1: TOPOLOGY FILE
# $2: STRING TO SEARCH
# $3: STRING TO INCLUDE
# USAGE: topology_add_include_specific sam.top ROT_NS.itp "#include \"postre.itp\""
function topology_add_include_specific () {
    ## DEFINING VARIABLES
    topology_file="$1"
    string_to_search="$2"
    string_to_include="$3"
    ## GETTING LINE NUMBER OF LAST INCLUDE
    line_num_last=$(grep -nE "$2" ${topology_file} | sed 's/\([0-9]*\).*/\1/' |tail -n1)
    ## ADDING TO THE TOPOLOGY FILE
    sed -i "$(( line_num_last+1 ))i$3" "$1" # #include \"$3\""
}

### FUNCTION TO ADD TO TOPOLOGY BEFORE SPECIFIC LOCATION
# $1: TOPOLOGY FILE
# $2: STRING TO SEARCH
# $3: STRING TO INCLUDE
# USAGE: topology_add_b4_include_specific sam.top ROT_NS.itp "#include \"postre.itp\""
function topology_add_b4_include_specific () {
    ## DEFINING VARIABLES
    topology_file="$1"
    string_to_search="$2"
    string_to_include="$3"
    ## GETTING LINE NUMBER OF LAST INCLUDE
    line_num_last=$(grep -nE "$2" ${topology_file} | sed 's/\([0-9]*\).*/\1/' |tail -n1)
    ## ADDING TO THE TOPOLOGY FILE
    sed -i "$(( line_num_last ))i$3" "$1" # #include \"$3\""
}

### FUNCTION TO CORRECT TOPOLOGY RESIDUE NUMBERS
# The purpose of this function is to correct for the total number of residues there are in a topology file.
# NOTE: This function only works for the last residue index found. You may have issues later when you have multiple residues defined! (how confusing)
## INPUTS:
#       $1: topology file
#       $2: residue_name
#       $3: correct number for residue name
## OUTPUTS:
#       Updated topology with the correct residue number
## USAGE: top_update_residue_num sam.top SOL 1
function top_update_residue_num () {
    ## DEFINING INPUTS
    top_file_="$1"
    residue_name_="$2"
    total_residues_="$3"
    
    ## CREATING TEMPORARY FILE
    temp_file="temp_file_.temp"
    
    ## REPLACING THE LAST RESIDUE NAME (done by using tac commands)
    tac ${top_file_} | sed "0,/${residue_name_}/{s/${residue_name_}.*/${residue_name_}   ${total_residues_}/}" | tac > ${temp_file}
    ## NEED TO FIX LATER THE LAST RESIDUE SHOULD BE DONE ONLY ONCE!
    
    ## REPLACING CURRENT TOPOLOGY, AND REMOVING TEMP
    cp -r "${temp_file}" "${top_file_}"
    rm "${temp_file}"

}

### FUNCTION TO CORRECT NUMBER OF RESIDUE NUMBERS FOR A TOPOLOGY
# The purpose of this function is to correct a specific topology residue number. We are assuming you know which molecule you want and how many atoms there are in a residue.
# NOTE: You may have issues later if your residue names are too similar. Be careful in selecting residue names that are not the same as water (e.g. SOL) for example.
## INPUTS:
#       $1: topology file
#       $2: gro file
#       $3: residue_name
#       $4: conversion betwen atom and residue numbers (e.g. for water: 3 atoms per one residue)
function top_correct_residue_name () {
    ## DEFINING INPUTS
    top_file_="$1"
    gro_file_="$2"
    residue_name_="$3"
    conversion_atom_to_res="$4"
    
    ## PRINTING
    echo -e "\n------ top_correct_residue_name ------"
    echo "Correcting topology for residue name: ${residue_name_}"
    echo "Conversion atom to residue: ${conversion_atom_to_res}"
    
    ## USING GREP TO FIND THE TOTAL NUMBER OF INSTANCES
    total_res_instances=$(grep "${residue_name_}" "${gro_file_}" | wc -l)
    
    ## CONVERTING TO TOTAL NUMBER OF RESIDUES
    total_residues=$(awk -v total_res=${total_res_instances} -v conv=${conversion_atom_to_res} 'BEGIN{ printf "%d", total_res/conv }')
    
    ## CORRECTING FOR TOPOLOGY
    top_update_residue_num "${top_file_}" "${residue_name_}" "${total_residues}"

}


### FUNCTION TO EXTRACT TOPOLOGY RESIDUE NUMBER
# The purpose of this function is to extract the total number of residues from a topology file. 
# INPUTS:
#       $1: topology file
#       $2: residue name
# OUTPUTS:
#   total residue number    
## USAGE:
#   top_extract_total_residue mixed_solv.top SOL
function top_extract_total_residue () {
    ## DEFINING INPUTS
    top_file_="$1"
    residue_name_="$2"
    
    ## DEFINING TEMPORARY FILE
    temp_file_="extract_total_res_temp.txt"
    
    ## FINDING MOLECULE LINE
    molecule_line=$(grep -nE '\[ molecules \]' "${top_file_}"  | sed 's/\([0-9]*\).*/\1/')
    
    ## TAILING AND GENERATING A TEMPORARY FILE
    tail -n+"${molecule_line}" "${top_file_}" > "${temp_file_}"
    
    ## FINDING LINE OF RESIDUE
    residue_num=$(grep -w "${residue_name_}" "${temp_file_}" | awk '{print $2}')
    
    ## PRINTING
    echo "${residue_num}"
    
    ## REMOVING TEMPORARY FILE
    if [ -e "${temp_file_}" ]; then
        rm "${temp_file_}"
    fi
}

### FUNCTION TO EXTARCT MULTIPLE RESIDUES
function top_extract_multiple_residues () {
    ## DEFINING INPUTS
    top_file_="$1"
    output_file_="${2:-output_multiple_residue.csv}" # output file
    
    ## FINDING PATH OF TOPOLOGY FILE
    topology_path="$(dirname $(dirname ${top_file_}))"
    
    ## DEFINING ARRAY RESIDUES
    declare -a residue_array=("SOL" "DIO" "dmso" "TET" "THF" "NMP" "GVLL")
    
    ## CREATING OUTPUT FILE AND ADDING HEADER
    if [ ! -e "${output_file_}" ]; then
        ## CREATING OUTPUT
        touch "${output_file_}"
        
        ## CREATING HEADER ARRAY
        header="File_name, "

        ## ADDING EACH RESIDUE TO HEADER
        for each_residue in ${residue_array[@]}; do
            header="${header}${each_residue}, "
        done
        
        ## ADDING FIRST HEADER
        echo "${header}" > "${output_file_}"
    fi
    
    data="${topology_path}, "
    ## LOOPING AND CHECKING EACH RESIDUE NAME
    for each_residue in ${residue_array[@]}; do
        ## FINDING TOTAL RESIDUE
        total_residue=$(top_extract_total_residue "${top_file_}" "${each_residue}")
        ## ADDING TO DATA
        data="${data}${total_residue}, "
    done
    ## AT THE END, SEND TO FIRST HEADER
    echo "${data}" >> "${output_file_}"

}



### FUNCTION TO DESIGNATE A SPLIT 
# The purpose of this function is to split a number of set windows, for instance, into chunks.
# For example, suppose you have a free energy calculation and you want to split it across
# multiple submission scripts. You will need to divide, say 17 windows, into 6 separate jobs.
# This function is designed to tell you which windows are required for each submission script.
# INPUTS:
#   $1: Starting window, e.g. 0
#   $2: Ending window, e.g. 17
#   $3: Number of splits, e.g. 3
# OUTPUTS:
#   _initial_window_array_: initial window array
#   _final_window_array_: final window array
## USAGE: 
#   > generate_split_array 0 17 6
#   > echo "${_initial_window_array_[@]}"
#   > echo "${_final_window_array_[@]}"
#   This function will generate output variables: _initial_window_array_, and _final_window_array_
function generate_split_array () {
    ## DEFINING INPUTS
    starting_value_="$1"
    ending_value_="$2"
    num_split="$3"
    
    ## FINDING TOTAL NUMBER OF WINDOWS
    n_total="$(( ${ending_value_}-${starting_value_}+1 ))"

    ## DIVIDING TOTAL WINDOWS WHILE ROUNDING DOWN
    divided_value="$((${n_total} / ${num_split}))"
    
    ## DEFINING INITIAL AND FINAL LAMBDAS
    _initial_window_array_=($(seq ${starting_value_} ${divided_value} ${ending_value_}))
    _final_window_array_=()
    
    ## LOOPING THROUGH INITIAL LAMBDA ARRAY TO GET FINAL LAMBDA
    for each_lambda in $(seq 0 $(( ${#_initial_window_array_[@]} - 1 ))); do
        ## SEEING IF LAMBDA IS THE FINAL VALUE
        if [[  "${each_lambda}" !=  "$(( ${#_initial_window_array_[@]} - 1 ))" ]]; then
            # SUBTRACT ONE FROM LAMBDA AND ADD IT TO FINAL
            final_lambda=$(( ${_initial_window_array_[each_lambda+1]}-1 ))
        ## FINAL LAMBDA, USE LAST VALUE
        else
            final_lambda="${ending_value_}"
        fi
        ## STORING FINAL LAMBDA
        _final_window_array_+=(${final_lambda})
    
    done

}

### FUNCTION TO FIND THE INTERSECTION BETWEEN TWO ARRAYS
# The purpose of this function is to find the intersection between two array variables in bash
## INPUTS:
#   $1: Number of elements in array 1
#   $2: Array 1
#   $3: Number of elements in array 2
#   #4: Array 2
## OUTPUTS:
#   Intersection as an array
## USAGE EXAMPLE: read -a residue_name_intersection <<< $(intersection_two_arrays "${#index_list[@]}" "${index_list[@]}" "${#l_res_names[@]}" "${l_res_names[@]}")
function intersection_two_arrays ()
{
    ## READ FOR MULTIPLE ARRAYS: https://stackoverflow.com/questions/10953833/passing-multiple-distinct-arrays-to-a-shell-function
    
    ##########################
    ### RE-CREATING ARRAYS ###
    ##########################
      declare -i num_args array_num;
      declare -a curr_args;
      while (( $# )) ; do
        curr_args=( )
        num_args=$1; shift
        while (( num_args-- > 0 )) ; do
          curr_args+=( "$1" ); shift
        done
        ## DEFINING ARRAY NUMBER
        array_num="$((++array_num))"
        # printf "$((++array_num)): %q\n" "${curr_args[@]}"
        ## STORING ARRAY 1 AND ARRAY 2
        if [ "${array_num}" == "1" ]; then
            ## COPYING ARRAY
            array_1=("${curr_args[@]}")
        elif [ "${array_num}" == "2" ]; then
            ## COPYING ARRAY
            array_2=("${curr_args[@]}")
        fi
        done
    ## PRINTING
    # echo "FINDING INTERSECTION BETWEEN TWO ARRAYS:"
    # echo "Array 1: ${array_1[@]}"
    # echo "Array 2: ${array_2[@]}"
    ###############################################
    ### FINDING INTERSECTION BETWEEN TWO ARRAYS ###
    ###############################################
    ## ADDING FRAME TO BLANKS
    l2=" ${array_1[*]} "
    ## LOOPING THROUGH EACH ITEM
    for item in ${array_2[@]}; do
      if [[ $l2 =~ " $item " ]] ; then    # use $item as regexp
        ## STORING THE ITEM
        result+=($item)
      fi
    done
    ## PRINTING THE ITEMS
    echo  ${result[@]}
}

### FUNCTION TO COMPUTE NUMBER OF CORES REQUIRED
# The purpose of this function is to compute the number of cores required. The assumption here is that we have 500 atoms / core and in SWARM, the maximum number of cores is 28. Therefore, we can compute the most efficient number of cores possible.
# INPUTS:
#   $1: total number of atoms
#   $2: num cores available (default: 28)
#   $3: num_atom per core (default: 500 / core)
# OUTPUTS:
#   variable: number of core to use
# USAGE: num_cores=$(compute_num_cores_required ${total_atoms})
function compute_num_cores_required () {
    ## DEFINING INPUTS
    num_atoms_="$1"
    num_cores_avail_="${2:-28}"
    num_atom_per_core_="${3:-500}"
    
    ## COMPUTING NUMBER OF CORES
    num_cores=$(awk -v n_atoms=${num_atoms_} -v n_atoms_per_core=${num_atom_per_core_} 'BEGIN{ printf "%d", n_atoms/n_atoms_per_core }')
    
    ## SEEING WHETHER THE CORES IS TOO MANY
    if [ "${num_cores}" -gt "${num_cores_avail_}" ]; then
        ## PRINTING ONLY AVAILABLE C OREES
        echo "${num_cores_avail_}"
    else
        echo "${num_cores}"
    fi
}


##################################
### GROMACS-SPECIFIC FUNCTIONS ###
##################################

### FUNCTION TO CHECK XTC FILE FOR THE TOTAL FRAMES
# The purpose of this function is to output the total frames in a given trajectory
function checkxtc() {
    # The purpose of this script is to take an xtc file and find the very last frame in time (ps).
    ## INPUTS: 
    #   $1 -- xtc file name
    #   $2 -- Name of your variable to store it in for the total frames in ps
    #   $3 -- name of your variable to store the final frame number (e.g. 100000.000 ps)
    ## OUTPUTS:
    #   current_xtc_equil_time: current xtc time
    ## USAGE: checkxtc ${equil_xtc} current_xtc_equil_time
                    
    # Defining local variable to keep my result
    local __lastFrameTime=$2
    local __final_frame_ps=$3
    
    # Defining files
    xtcfile="$1"
    tempxtcFile="checkxtc.check"
    
    # Creating temporary file using gmx check (standard error output is 2)
    gmx check -f "$1" 2> $tempxtcFile
    
    # Now, let's find the value
    currentFrame=$(grep "Step" ${tempxtcFile} | awk {'print $2'}) # total frames
    timestep=$(grep "Step" ${tempxtcFile} | awk {'print $NF'}) # dt in ps
    
    # Calculating last frame based on the current frame and time steps
    last_frame_ps=$(awk -v nFrames=$currentFrame -v dt=$timestep 'BEGIN {printf "%.3f",(nFrames-1)*dt; exit}' )
    ## FINDING FINAL FRAME
    final_frame_ps=$(grep "Last frame" ${tempxtcFile} | awk {'print $NF'})
    
    eval $__lastFrameTime=${last_frame_ps}
    eval $__final_frame_ps=${final_frame_ps}
    
    # Deleting the temporary file
    if [ -e $tempxtcFile ]; then
        rm $tempxtcFile
    fi
    
    # Printing
    echo "Looking at the XTC file: $1; which has ${last_frame_ps} ps worth of data"
    echo "Last frame is: ${final_frame_ps}"
}

### FUNCTION TO MEASURE BOX SIZE BASED ON GRO FILE
## The purpose of this function is to measure the box size of a gro file.
## INPUTS:
#   $1: gro file
## OUTPUTS:
#   array: array of length 3 with x, y, z components
## USAGE: read -a box_size <<< $(gro_measure_box_size "${gro_file}")
#   ACCESS X, Y, Z BOX LENGTH BY:
#   box_x=${array[0]}
#   box_y=${array[1]}
#   box_z=${array[2]}
function gro_measure_box_size ()
{
    ## DEFINING INPUTS
    input_gro_file_="$1"
    
    ## GET BOX VECTORS AND STORE THEM INTO COMPONENTS
    output=$(tail -n 1 ${input_gro_file_})
    ## CREATING ARRAY OF THE OUTPUTS
    array=($output)
    
    ## OUTPUTTING
    echo "${array[@]}"
}

### FUNCTION TO COMPUTE VOLUME FROM GRO FILE
# INPUTS:
#       $1: input gro file
# OUTPUTS:
#   volume in nm3
# USAGE:
#   volume=$(gro_compute_volume gro_file)
function gro_compute_volume () {
    ## INPUTS
    input_gro_file="$1"
    
    ## READING BOX VOLUME
    read -a box_size <<< $(gro_measure_box_size "${input_gro_file}")
    ## FINDING VOLUME
    volume="$(compute_volume_from_box_length "${box_size[@]}")"
    
    ## PRINTING VOLUME
    echo "${volume}"

}

### FUNCTION TO COMPUTE NUMBER DENSITY BASED ON A GRO FILE
# INPUTS:
#   input_num_solvents: number of solvents that are within the gro file
#   input_gro_file: input gro file
# OUTPUTS:
#   num_density: number density = num / volume
# USAGE:
#   gro_compute_num_density 216 gro_file
function gro_compute_num_density () {
    ## INPUTS
    input_num_solvents="$1" # number of solvents
    input_gro_file="$2" # gro file
    
    ## COMPUTING VOLUME
    volume=$(gro_compute_volume "${input_gro_file}")
    
    ## COMPUTING NUMBER DENSITY
    num_density=$(awk -v num=${input_num_solvents} -v vol=${volume} 'BEGIN{ printf  "%.4f",num/vol}')
    
    ## PRINTING
    echo "${num_density}"

}

### FUNCTION TO FIND THE TOTAL NUMBER OF GRO FILES
## The purpose of this function is to find the total number of atoms for a gro file
## INPUTS:
#   $1: gro files
## OUTPUTS:
#   total number of atoms for the gro file
## USAGE: total_atoms=$(gro_count_total_atoms gro_file)
function gro_count_total_atoms () {
    ## DEFINING INPUTS
    input_gro_file_="$1"
    
    ## FINDING TOTAL NUMBER OF ATOMS
    num_atoms=$(head -n 2 "${input_gro_file_}" | tail -n 1) # Gets total number of water molecules
    ## TRIMMING ANY EXCESS
    trim="${num_atoms%"${num_atoms##*[![:space:]]}"}"
    #remove leading whitespace too
    trim="${trim#"${trim%%[![:space:]]*}"}"
    ## PRINTING
    echo "${trim}"
}


### FUNCTION TO COMBINE TWO GRO FILES
## The purpose of this function is to combine two gro files. Note that this uses gmx commands, so gromacs must be installed. Note that you may need to shift each of the gro files to ensure that it is placed correctly. 
## INPUTS:
#   $1: gro file 1
#   $2: gro file 2
#   $3: output gro file name
## USAGE:  gro_combine_two_files beta_lactoglobulin_enlarge.gro np_no_water_center_enlarge.gro combined.gro
function gro_combine_two_files ()
{   
    ## DEFINING INPUTS
    input_gro_file_1_="$1"
    input_gro_file_2_="$2"
    
    ## DEFINING OUTPUTS
    output_gro_file_="$3"
    
    ## DEFINING TEMPORARY FILES
    temp_file_1_="${input_gro_file_1_%.gro}_temp.gro"
    temp_file_2_="${input_gro_file_2_%.gro}_temp.gro"
    
    ## CREATING TEMPORARY FILES
    cp "${input_gro_file_1_}" "${temp_file_1_}"
    cp "${input_gro_file_2_}" "${temp_file_2_}"
    
    ## FINDING TOTAL NUMBER OF ATOMS FOR EACH GRO FILE
    gro_1_total_atoms_="$(gro_count_total_atoms "${input_gro_file_1_}")"
    gro_2_total_atoms_="$(gro_count_total_atoms "${input_gro_file_2_}")"
    
    ## FINDING TOTAL NUMBER OF WATER MOLECULES
    total_num_atoms_=$((${gro_1_total_atoms_} + ${gro_2_total_atoms_}))
    
    ## ADJUSTING GRO FILES
    # REMOVING FIRST TWO LINES OF FILE 2
    sed -i '1d' "${temp_file_2_}"
    sed -i '1d' "${temp_file_2_}"
    # REMOVING LAST LINE OF FILE 1
    sed -i '$d' "${temp_file_1_}"
    
    ## COMBINING
    cat "${temp_file_1_}" "${temp_file_2_}" > "${output_gro_file_}"

    ## REPLACE NUMBER OF ATOMS WITH COMBINED NUMBER
    sed -i -e 0,/${gro_1_total_atoms_}/s/${gro_1_total_atoms_}/${total_num_atoms_}/ "${output_gro_file_}"
    
    ## RENUMBER USING GENCONF -- NO CENTERING
    gmx genconf -f "${output_gro_file_}" -o "${output_gro_file_}" -renumber &>/dev/null
    
    ## REMOVING TEMPORARY FILES
    rm "${temp_file_1_}" "${temp_file_2_}"
}

### FUNCTION TO COMBINE TWO FILES BY THE ZTH DIMENSION
# This function combines gro files based on z dimensions. 
# This is useful for lots of systems, such as NP-lipid bilayers systems
# NOTES:
#   - First gro file will be the top
#   - Second gro file will be the bottom
#   - Buffer zone is used to prevent overlap of the two gro file systems
# INPUTS:
#   $1: input gro file 1 (top)
#   $2: input gro file 2 (bottom)
#   $3: output gro file
#   $4: buffer z-value
# OUTPUTS:
#   new gro file with combined gro 1 and 2
# USAGE:
### COMBINING GRO FILES
#gro_combine_two_files_by_z "${np_gro_5}" \
#                           "${output_lipid_membrane_prefix}.gro" \
#                           "${combine_file_1}" \
#                           "${buffer_zone}"
function gro_combine_two_files_by_z () {
    ## DEFINING INPUTS
    input_gro_file_1_="$1"
    input_gro_file_2_="$2"
    output_gro_file="${3-combined.gro}"
    buffer_z_="${4:-0}" # zero by default
    
    ## DEFINING TEMP FILE
    temp_file_1="${input_gro_file_1_%.gro}_temp_resize.gro"
    temp_file_2="${input_gro_file_2_%.gro}_temp_resize.gro"
    
    ## GETTING DIMENSION
    read -a gro_file_1_dims_ <<< $(read_gro_box_lengths "${input_gro_file_1_}")
    read -a gro_file_2_dims_ <<< $(read_gro_box_lengths "${input_gro_file_2_}")
    
    ## GETTING TOTAL BOX SIZE
    total_box_size=$(awk -v file_1_z=${gro_file_1_dims_[2]} \
                         -v file_2_z=${gro_file_2_dims_[2]} \
                         -v buffer_z=${buffer_z_} \
                             'BEGIN{ printf "%.5f", \
                             file_1_z + file_2_z + buffer_z}')
                             
    ## GETTING TOTAL TRANSLATION
    translate_z=$(awk -v z_1_end=${gro_file_2_dims_[2]} \
                      -v buffer_z=${buffer_z_} \
                       'BEGIN{ printf "%.5f", \
                         z_1_end + buffer_z}')
                             
    ## GMX EDITCONF FOR FILE 1
    gmx editconf -f "${input_gro_file_1_}" \
                 -o "${temp_file_1}" \
                 -noc \
                 -box "${gro_file_2_dims_[0]}" \
                      "${gro_file_2_dims_[1]}" \
                       "${total_box_size}" \
                 -translate 0 0 "${translate_z}" &>/dev/null
    
    ## GMX EDIT CONF FOR FILE 2
    gmx editconf -f "${input_gro_file_2_}" \
                 -o "${temp_file_2}" \
                 -noc \
                 -box "${gro_file_2_dims_[0]}" \
                      "${gro_file_2_dims_[1]}" \
                       "${total_box_size}" &>/dev/null
                       
    ## COMBINING GRO FILES
    gro_combine_two_files "${temp_file_1}" \
                          "${temp_file_2}" \
                          "${output_gro_file}"
                          
    ## REMOVING TEMP FILES
    rm -f "${temp_file_1}" "${temp_file_2}"
}


################################################
################## MDP_FILES  ##################
################################################

### FUNCTION TO FIND MDP OPTIONS
# The purpose of this function is to find mdp options. We are assuming here that you have some 'options = value', where value represents the third column via awk
# INPUTS:
#   $1: mdp file
#   $2: variable that you are looking for (e.g. nsteps)
# OUTPUTS:
#   Values that corresonds to options.
# USAGE:
#   mdp_nsteps=$(mdp_find_option_values ${mdp_file} "nsteps")
function mdp_find_option_values () {
    ## DEFINING INPUTS
    input_mdp_file_="$1"
    input_mdp_options_="$2"
    
    ## USING GREP TO FIND THE FINAL VALUE
    echo "$(grep "$2" "$1"  | awk '{print $3}')"
    
}

### FUNCTION TO ADD ENERGY GROUPS TO MDP FILE
# The purpose of this function is to add energy groups to your MDP file assuming that no energy groups have been presented. As a result, we can simply tag on your mdp file. Note that we will copy over the mdp file from a currently existing one.
# INPUTS:
#   $1: input_mdp_file
#   $2: output_mdp_file
#   $3: enery groups, e.g. 'SOLUTE SOLVENT'
# OUTPUTS:
#   output mdp file with energy groups, e.g.
#       energygrps                       = SOLUTE SOLVENT
# USAGE: mdp_add_energy_groups "${input_mdp_file}" "${output_mdp_file}" "${energy_group_name}"
function mdp_add_energy_groups () {
    ## DEFINING INPUTS:
    input_mdp_file_="$1"
    output_mdp_file_="$2"
    energy_groups_="$3"
    
    ## COPYING OVER MDP FILE
    cp -r "${input_mdp_file_}" "${output_mdp_file_}"
    
    ## ADDING TO MDP FILE
    echo "" >> "${output_mdp_file_}"
    echo "; ENERGY GROUPS " >> "${output_mdp_file_}"
    echo "energygrps          = ${energy_groups_}" >> "${output_mdp_file_}"
    
    ## PRINTING
    echo "----- mdp_add_energy_groups -----"
    echo "--> Created ${output_mdp_file_} from ${input_mdp_file_}"
    echo "--> Added enery groups: ${energy_groups_}"
}


### FUNCTION TO COMMENT OUT INFORMATION IN MDP
# The purpose of this function is to comment out information in the MDP file.
# This is useful when you are trying things like commenting out freeze groups. 
# INPUTS:
#   $1: input_mdp_file
#   $2: comment out item
# USAGE:
#   mdp_comment_out_specific_flags "nvt_double_prod_gmx5_charmm36_frozen_gold.mdp" "freezegrps"
function mdp_comment_out_specific_flags () {
    ## DEFINING INPUTS
    input_mdp_file_="$1"
    comment_out_name_="${2:-freezegrps}"
    
    ## LOOKING FOR ITEM
    read -a locations <<< $(grep -nE "^${comment_out_name_}" "${input_mdp_file_}" | sed 's/\([0-9]*\).*/\1/')
    
    ## LOOPING THROUGH EACH LINE NUMBER
    for line_num in ${locations[@]}; do
        echo "Editing line ${line_num} in ${input_mdp_file_}"
        ## GETTING LINE TEXT
        line_text="$(sed -n ${line_num}p ${input_mdp_file_})"
        ## ADDING COMMENT
        line_text_with_comment="; ${line_text}"
        
        ## REPLACING THE LINE
        sed -i "s#${line_text}#${line_text_with_comment}#g" "${input_mdp_file_}"
        
    done

}

########################################################
################## ITP_FILE FUNCTIONS ##################
########################################################

### FUNCTION TO TURN OFF ITP CHARGES AND OUTPUT A NEW ITP
# The purpose of this function is to turn off charges for a specific itp file. The idea here is that we want a itp file with no charges to calculate energies, ignoring electrostatic interactions.
# INPUTS:
#   $1: input_itp_file
#   $2: output_itp_file
# OUTPUTS:
#   output itp file with no charges
# USAGE:
#   convert_itp_nocharge acetone.itp acetone_nocharge.itp
function convert_itp_nocharge () {
    ## DEFINING INPUTS
    input_itp_file_="$1"
    output_itp_file_="$2"
    
    ## DEFINING COLUMN OF CHARGE
    charge_column=7
    charge_value="0.000"
    
    ## FINDING ATOM LINE
    line_atom=$(grep -nE '\[ atoms \]' "${input_itp_file_}"  | sed 's/\([0-9]*\).*/\1/')
    line_atom_plus_one=$(( ${line_atom} + 1 )) # Used to ignore the atoms directive
    
    ## CREATING COPY OF INPUT TO OUTPUT
    cp "${input_itp_file_}" "${output_itp_file_}"
    
    # INDEX NUMBER
    idx="${line_atom_plus_one}"
    first_line_num=""
    re='^[0-9]+$' # Regular expression for numbers
    
    ## PRINTING
    echo "----- convert_itp_nocharge -----"
    echo -e "\n----- FINDING FIRST AND LAST LINE OF THE ATOM DIRECTIVE -----"
    ## READING LINE
    while read line; do
        ## PRINTING
        echo "${line}"
        
        ## FINDING FIRST LINE NUMBER THAT HAS "1" IN FIRST COLUMN
        if [ $(echo "${line}" | awk {'print $1'} ) == "1" ]; then
            first_line_num="${idx}"
            echo "First line found: ${first_line_num}"
        fi
        
        ## NOW THAT FIRST LINE IS FOUND, CONTINUING UNTIL THE PRINT FUNCTION BREAKS
        if [[ ! -z "${first_line_num}" ]]; then
        
            if ! [[ $(echo "${line}" | awk {'print $1'} ) =~ $re ]] ; then
                last_line_num="$(( ${idx} - 1 ))" # CORRECTING FOR LAST NUMBER
                echo "Last line found: ${last_line_num}"
                echo "Stopping while loop"
                break
            fi
        
        fi
        
        ## ADDING INDEX NUMBER
        idx=$(( ${idx}+1 ))
        
    done < <(tail -n +${line_atom_plus_one} ${input_itp_file_})
    
    echo -e "--------------------------------------------------------------------------\n"
    echo "Finding complete! First line number: ${first_line_num}, last line number: ${last_line_num}"
    echo "First line: $(sed -n -e ${first_line_num}p ${input_itp_file_})"
    echo "Last line: $(sed -n -e ${last_line_num}p ${input_itp_file_})"    
    echo "Check if this is correct!!!"
    # sleep 2
    
    ## CREATING FILE
    echo -e "\n--> Creating ${output_itp_file_}..."
    
    ## USING HEAD TO GET THE TOP OF THE FILE
    head -n $(( ${first_line_num}-1 )) "${input_itp_file_}" > "${output_itp_file_}"
    
    ## LOOPING THROUGH EACH LINE, EDITING, THEN ADDING
    while read line; do
        # echo "Adding no charge line: $(echo "${line}" | awk '$7="0.000"')"
        echo "Adding no charge line: $(echo "${line}" | awk '{print gensub (/[^[:blank:]]+/, v, 7)}' v="0.000")"
        echo "${line}" | awk '{print gensub (/[^[:blank:]]+/, v, 7)}' v="0.000" >> "${output_itp_file_}"
        
    done < <(sed -n "${first_line_num},${last_line_num}p" "${input_itp_file_}")
    
    ## USING TAIL TO GET THE REST OF THE FILE
    tail -n +$(( ${last_line_num}+1 )) "${input_itp_file_}" >> "${output_itp_file_}"
}

### FUNCTION TO CORRECT GENRESTR INDEX VALUES
# The purpose of this script is to correct genrestr's indexing. Normally, you 
# get index values relative to the gro file. However, this is non-ideal when 
# you want to restraint a particular molecule, since you will need the 
# [ position_restraints ] portion to be directly after the itp file.
# INPUTS:
#   $1: posre.itp file -- input posre.itp
# OUTPUTS:
#   Updated posre itp file
# USAGE: itp_fix_genrestr_numbering solute_posre.itp
function itp_fix_genrestr_numbering () {
    ## DEFINING INPUTS
    input_posre_file_="$1"
    
    ## DEFINING TEMP FILE
    temp_posre_file="temp_posre.temp"
    
    ## FINDING POSITION RESTRAINT LINE
    pos_rest_line=$(grep -nE "\[ position_restraints \]" "${input_posre_file_}" | sed 's/\([0-9]*\).*/\1/')
    echo "${pos_rest_line}"
    
    ## CREATING TEMP FILE
    sed -n "1,${pos_rest_line}p" "${input_posre_file_}" > ${temp_posre_file}
    
    ## ADDING COMMENT LINE
    echo ";  i funct       fcx        fcy        fcz" >> ${temp_posre_file}
    
    ## DEFINING COUNTER
    atom_counter="1"
    
    ## USING WHILE LOOP
    while read line; 
    do 
    echo ${line} | awk '{$1='${atom_counter}'; print ;}' >> ${temp_posre_file}
    # echo ${line} | sed '1s/^/1/'
    
    ## ADDING TO COUNTER
    atom_counter=$(( ${atom_counter} + 1 ))
    
    done <<< "$(sed '/^;/ d' ${input_posre_file_} | sed -e '1,/\[ position_restraints \]/ d' )" # <-- Reads itp file, removes all comments, looks for everything following position restraints
    
    ## COPYING TEMP AND REMOVING IT
    cp -r "${temp_posre_file}" "${input_posre_file_}"
    if [ -e "${temp_posre_file}" ]; then
        rm "${temp_posre_file}"
    fi

}

### FUNCTION TO EXTRACT NAME OF RESIDUE GIVEN ITP FILE
# This function looks into your itp file and extracts the residue name from it
# $1: ITP file
# Usage: itp_get_resname itp_file_name
function itp_get_resname () {
    itp_file="${1%.itp}"
    # Defining itp file
    itp_file="${itp_file}.itp"

    # Finding residue name
    resname=$(sed '/^;/ d' "${itp_file}" | grep -E -A 1 '\[ moleculetype \]' | tail -n1 | awk '{print $1}') # Removes comments, finds molecule type, looks at second line, then prints first column
    
    # Printing to variable
    echo "${resname}"
}

### FUNCTION TO CREATE A PRM FILE FROM ITP FILE
# The purpose of this function is to create a parameter 
# file using the itp file. We assume you may have messed up and 
# had atomtypes, etc. inside an itp file, which should have been 
# outside in a parameter file.
# INPUTS:
#   $1: itp_file
# OUTPUTS:
#   Updated itp file with parameters as a .prm file
## USAGE:
#   extract_prm_from_itp sam.itp
function itp_extract_prm () {
    ## DEFINING INPUTS
    input_itp_file_="$1"
    
    ## GETTING PARMATER FILE
    parameter_file="${input_itp_file_%.itp}".prm
    
    ## DEFINING COPY OF ITP FILE
    itp_file_copy="__temp_itp_file.itp"
    
    ## COPYING ITP FILE
    cp "${input_itp_file_}" "${itp_file_copy}"
    
    ## FINDING INDEX OF MOLECULE TYPE
    line_num_molecule_type=$(grep -nE '\[ moleculetype \]' "${input_itp_file_}" | sed 's/\([0-9]*\).*/\1/')
    
    ## GENERATING TWO FILES
    head -n $((line_num_molecule_type-1)) "${itp_file_copy}" > "${parameter_file}"
    
    ## ITP FILE
    tail "${itp_file_copy}" -n+"${line_num_molecule_type}" > "${input_itp_file_}"
    
    ## REMOVING TEMP FILE
    if [ -e "${itp_file_copy}" ]; then
        rm "${itp_file_copy}"
    fi
    
}



################################################
################## TRAJECTORY ##################
################################################

### FUNCTION TO CHECK XTC TOTAL TIME IN PS
# The purpose of this script is to take an xtc file and find the very last frame in time (ps).
# INPUTS:
#   $1 -- xtc file name
#   $2 -- Name of your variable to store it in
## OUTPUTS:
#   Variable last frame time
## USAGE:
#   xtc_check_total_time_ps ${prod_xtc} current_prod_xtcTime
function xtc_check_total_time_ps_updated() {
    # Defining local variable to keep my result
    local __lastFrameTime=$2

    # Defining files
    xtcfile="$1"
    tempxtcFile="checkxtc.check"

    # Creating temporary file using gmx check (standard error output is 2)
    gmx check -f "$1" 2> ${tempxtcFile}
    
    ## WAITING FOR GMX CHECK TO COMPLETE
    wait
    
    ## CHECKING TO SEE IF THERE IS A TEMP FILE
    if [ -e ${tempxtcFile} ]; then
    
        # Now, let's find the value
        currentFrame=$(grep "Step" ${tempxtcFile} | awk {'print $2'}) # total frames
        timestep=$(grep "Step" ${tempxtcFile} | awk {'print $NF'}) # dt in ps

        # Calculating last frame
        last_frame_ps=$(awk -v nFrames=${currentFrame} -v dt=${timestep} 'BEGIN {printf "%.3f",(nFrames-1)*dt; exit}' )
        
        ## WAITING FOR LAST FRAME TO BE FOUND
        wait
        
        eval $__lastFrameTime=${last_frame_ps}
    else
        echo "Error! ${tempxtcFile} does not seem to exist!"
        echo "Pausing here for 3 seconds so you can see this error message"
        echo "Check xtc_check_total_time_ps_updated command in server_general_research_functions.sh"
        sleep 3
    fi
    
    ## Printing
    echo "Looking at the XTC file: $1; which has ${last_frame_ps} ps worth of data"
    
    ## Deleting the temporary file
    if [ -e ${tempxtcFile} ]; then
        rm ${tempxtcFile}
    fi
}

### FUNCTION TO CHECK EDR FILE
# The purpose of this function is to check the time frames for an .edr file
# INPUTS:
#   $1: edr file
#   $2: time frame variable
# OUTPUTS:
#   Last  frame time
## USAGE:
#   edr_check_total_time_ps ${edr_file} current_prod_time
function edr_check_total_time_ps() {
    ## DEFINING LOCAL VARIABLE
    local __lastFrameTime=$2
    ## DEFINING VALUES
    file_="$1"
    tempFile="checkedr.check"
    
    # Creating temporary file using gmx check (standard error output is 2)
    gmx check -e "${file_}" 2> ${tempFile}
    
    ## EXTRACTING LAST FRAME
    last_frame_ps=$(grep "Last energy frame" "${tempFile}" | awk '{print $NF'})
    
    ## EVALUATING LAST  FRAME
    eval $__lastFrameTime=${last_frame_ps}

    # Deleting the temporary file
    if [ -e ${tempFile} ]; then
        rm ${tempFile}
    fi

    # Printing
    # echo "Looking at file: $1; which has ${last_frame_ps} ps worth of data"
}

### FUNCTION TO CHECK IF EDR FILE HAS SAME SIMULATION TIME AS SUGGESTED
# The purpose of this function is to compute the total time frame available for a *.edr file. This function will output either "True" or "False".
## INPUTS:
#   $1: edr file
#   $2: time to check in ps (e.g. 300000.000) -- note the decimal places
## OUTPUTS:
#   variable: true or false
## USAGE:
#   edr_log=$(edr_logical_same_time_ps $edr_file 300000.000)
#   will output true if your edr file matches the time, otherwise it is false.
function edr_logical_same_time_ps () {
    ## DEFINING INPUTS
    edr_file_="$1"
    time_check_ps_="$2"
    
    ## RUNNING ONLY IF EDR FILE EXISTS
    if [ -e "${edr_file_}" ]; then
        ## COMPUTING TIMES
        edr_check_total_time_ps ${edr_file_} edr_file_time_ps

        ## SEEING EQUIVALENCE
        check_time_equality=$(awk 'BEGIN{ print "'${time_check_ps_}'"!="'${edr_file_time_ps}'" }')

        ## USING IF LOGICALS
        if [ "${check_time_equality}" -eq 1 ];then 
            # IF 1, THEN NOT EQUAL
            echo "false"
        else
            # Otherwise, they are equal
            echo "true"
        fi
    else
        echo "false" # No EDR file
    fi
}

### FUNCTION TO CHECK IF GENEREALIZD FILE HAS SAME SIMULATION TIME AS SUGGESTED
# The purpose of this function is to compute the total time frame available for a *.edr file. This function will output either "True" or "False".
## INPUTS:
#   $1: any type of file
#   $2: time to check in ps (e.g. 300000.000) -- note the decimal places
#   $3: variable to store the true/false logical in
## OUTPUTS:
#   variable: true or false
## USAGE:
#   check_logical_same_time_ps edr_file 11000.000 variable
#   edr_log=$(edr_logical_same_time_ps $edr_file 300000.000) <-- depreciated
#   will output true if your edr file matches the time, otherwise it is false.
function check_logical_same_time_ps () {
    ## DEFINING INPUTS
    file_="$1"
    time_check_ps_="$2"
    
    ## DEFINING LOCAL VARIABLE
    local __logical_same_time="$3"
    
    ## CHECKING EXTENSION
    extension_="${file_##*.}"
    
    ## RUNNING ONLY IF EDR FILE EXISTS
    if [ -e "${file_}" ]; then
        ## COMPUTING TIMES
        if [[ "${extension_}" == "edr"  ]]; then
            edr_check_total_time_ps ${file_} file_time_ps
        elif [[ "${extension_}" == "xtc"  ]]; then
            xtc_check_total_time_ps_updated ${file_} file_time_ps
        else
            echo "Error! Extension (${extension_}) does not have a reading time value!"
            echo "Please check to see if you can include the extension within the check_logical_same_time_ps code"
            echo "Stopping here to prevent subsequent errors..."
            sleep 5
            exit 1
        fi
        
        ## WAITING
        wait
        ## SEEING EQUIVALENCE
        check_time_equality=$(awk 'BEGIN{ print "'${time_check_ps_}'"!="'${file_time_ps}'" }')

        ## USING IF LOGICALS
        if [ "${check_time_equality}" -eq "1" ];then 
            # IF 1, THEN NOT EQUAL
            eval $__logical_same_time="False"
        else
            eval $__logical_same_time="True"
        fi
    else
        eval $__logical_same_time="False"
    fi
}

##############################################
################## TOPOLOGY ##################
##############################################

### FUNCTION TO CLEAR TOPOLOGY MOLECULES SECTION
## The purpose of this function is to clear the [ molecules ] section for any topology file
## INPUTS:
#   $1: topology file
## OUTPUTS:
#   Updated topology file with [ molecules ] cleared
## USAGE EXAMPLE: top_clear_molecules sam.top
function top_clear_molecules () {
    ## DEFINING INPUTS
    input_top_file_="$1"
    
    ## DEFINING TEMPORARY TOPOLOGY FILE
    temp_top_file_="${input_top_file_%.top}_temp.top"
    
    ## FINDING LINE WHERE MOLECULES SHOW UP
    line_num_molecules="$(grep -nE '\[ molecules \]' "${input_top_file_}" | sed 's/\([0-9]*\).*/\1/')"
    
    ## REMOVING ALL LINES
    head -n"${line_num_molecules}" "${input_top_file_}" > "${temp_top_file_}"
    
    ## COPYING OVER TOP FILE AND REMOVING THE TEMPORARY
    cp -r "${temp_top_file_}" "${input_top_file_}"
    rm "${temp_top_file_}"
    
    ## ADDING COMMENTS FOR THE TOPOLOGY FILE
    echo "; compound   nr_mol" >> "${input_top_file_}"
}

### FUNCTION TO READ TOPOLOGY MOLECULES (NO COMMENTS)
# The purpose of this function is to read topology file molecules section
## INPUTS:
#   $1: topology file
function top_read_molecules () {
    ## DEFINING INPUTS
    input_top_file_="$1"
    
    ## FINDING LINE WHERE MOLECULES SHOW UP
    line_num_molecules="$(grep -nE '\[ molecules \]' "${input_top_file_}" | sed 's/\([0-9]*\).*/\1/')"
    
    ## ADDING ONE TO LINE NUMBER
    line_num_molecules=$((${line_num_molecules}+1))
    
    ## PRINTING AND IGNORING COMMENTS (;)
    sed -n ''${line_num_molecules}',$p' "${input_top_file_}" | grep -vE "^;"

}

####################################
##### MOLECULAR INFO FUNCTIONS #####
####################################

## FUNCTION TO CHECK MOLECULAR VOLUME TO SEE IF ENTRY EXISTS
# INPUTS:
#   $1: name of the molecule
#   $2: path of the file
# OUTPUTS:
#   True/False
# USAGE: read logical <<< $(check_molecule_info $MOLECULE_NAME $FILE_PATH)
function check_molecule_info () {
    # DEFINING VARIABLES
    molecule_name="$1"
    file_path="$2"
    
    # GREPPING TO SEE EXISTANCE IN FIRST LINE
    line_location=$(grep -nE "^${molecule_name}" "${file_path}" | sed 's/\([0-9]*\).*/\1/')
    
    # SEEING IF MOLECULAR INFORMATION IS THERE!
    if [ -z "${line_location}" ]; then
        echo "False" # LINE DOES NOT EXIST
    else
        echo "True" # LINE DOES EXIST
    fi
}


### THE PURPOSE OF THIS FUNCTION IS TO WRITE INTO MOLECULAR INFO INTO A FILE
# INPUTS:
#   $1: molecular_file: File of interest
#   $2: name: Name of molecule you are interested in
#   $3: value: Value accompanied with the molecule
# USAGE: write_molecular_info $MOLECULAR_FILE $SOLUTE_NAME $VALUE
function write_molecular_info () {
    # DEFINING VARIABLES
    molecular_file="$1"
    name="$2"
    value="$3"

    # READING FILE TO SEE IF MOLECULE EXISTS
    read within_molecule_vol <<< $(check_molecule_info ${name} ${molecular_file})
    if [ "${within_molecule_vol}" == "False" ]; then
        echo "Adding line to ${molecular_file}"
        printf "%s, %0.2f\n" "${name}" "${value}" >> "${molecular_file}"
    else
        echo "${molecular_file} already has ${value} -- doing nothing"
    fi
}

### FUNCTION TO CALCULATE MOLECULAR WEIGHT
# INPUTS:
#   $1: PDB file
#   $2: reference file which has all atom names and molecular weights
# OUTPUTS:
#   Value for molecular weight
## USAGE: calc_pdb_molecular_weight urea_ini.pdb
function calc_pdb_molecular_weight () {
    ## DEFINING INPUT VARIABLE
    input_pdb="$1" # INPUT PDB
    ref_file=${2-"/home/akchew/scratch/SideProjectHuber/prep_mixed_solvent/input_files/solvents/molecular_weights_ref.txt"} # REFERENCE FILE

    ## GETTING THE ATOM LIST
    atom_list=$(head -n -1 "${input_pdb}" | awk '{print $3}' | sed 's/[0-9]*//g')
    
    ## GETTING REFERENCE FILE DATA WITHOUT COMMENTS
    ref_file_data="$(grep "^[^;]" ${ref_file})"
    
    ## DEFINING MOLECULAR WEIGHT VARIABLE
    molecular_weight=0
    
    ## LOOPING THROUGH THE LIST AND GETTING MOLECULAR WEIGHT, THEN ADDING IT
    while read -r current_atom; do
        ## FINDING THE ATOMIC WEIGHT
        # Removes all comments | finds atom name | prints mass | cuts all white space
        atom_weight=$(grep "^[^;]" "${ref_file}" | grep "${current_atom}" | awk -F',' '{print $2}' | tr -d ' ')        
        
        ## ADDING TO MOLECULAR WEIGHT
        # Adds the weights, removes all issues with windows, then makes math happen
        molecular_weight=$(echo "${molecular_weight} + ${atom_weight}" | tr -d $'\r' | bc )
    done <<< "${atom_list}"
    ## PRINTING MOLECULAR WEIGHT
    echo "${molecular_weight}"
}


### FUNCTION TO ROUND TO NEAREST INTEGER
# The purpose of this function is to round to the nearest integer. For instances, 
# 7.3 would get you 8 and so on.
# INPUTS:
#   $1: value you want to round
# OUTPUTS:
#   rounded value
# USAGE:
#   round_nearest_integer 5
#   rounded_value=$(round_nearest_integer 5)
function round_nearest_integer () {
    ## DEFINING INPUTS
    value_="$1"
    ## PRINTING
    echo "a=${value_}; b=1; if ( a%b ) a/b+1 else a/b" | bc
    # Note, b is 1 for nearest integer. In theory, this should work for any division.
}

### FUNCTION TO GET FIRST RESIDUE NAME IN GRO FILE
# The purpose of this function is to get the first residue name in a gro file. 
# INPUTS:
#   $1: gro file path
# OUTPUTS:
#   $2: first residue name from gro file
# USAGE:
#   gro_get_first_residue_name gro_file_name.gro
#   res_name=$(gro_get_first_residue_name gro_file_name.gro)
function gro_get_first_residue_name  () {
    ## DEFINING INPUTS
    gro_file_="$1"
    ## GETTING RESIDUE NAME
    first_residue_name_=$(head -n3 "${gro_file_}"  | tail -n1 | cut -c 6-10 | tr -d '[:space:]')
    ## PRINTING RESIDUE NAME
    echo "${first_residue_name_}"
}

### FUNCTION TO GET RESIDUE NAME WITHIN ATOMS OF ITP FILE
# The purpose of this function is to get the residue names of a
# itp file. 
# INPUTS:
#   $1: itp file, which you're looking for "LIG" as shown below
#   [ atoms ]
# ;   nr  type  resi  res  atom  cgnr     charge      mass       ; qtot   bond_type
#      1   c3     1   LIG     C    1    -0.095100     12.01000 ; qtot -0.095
# OUTPUTS:
#   residue name
# USAGE:
#   itp_get_residue_name_within_atom 10_mer_PE.itp
function itp_get_residue_name_within_atom () {
    ## DEFINING INPUTS
    input_itp_file_="$1"
    
    ## DEFINING TEMPORARY FILE
    temp_file_="${input_itp_file_%.itp}_temp_file.itp"
    ## GETTING ITP WITHOUT COMMENTS
    grep -o '^[^;]*' ${input_itp_file_} > "${temp_file_}"
    
    ## FINDING ATOM
    atom_linenum=$(grep -nE '\[ atoms \]' "${temp_file_}" | sed 's/\([0-9]*\).*/\1/')
    
    ## GETTING RESIDUE NAME
    res_name_=$(sed -n "$((${atom_linenum}+1))p" "${temp_file_}"  | awk '{print $4}')
    
    ## REMOVING TEMP FILE
    if [ -e "${temp_file_}" ]; then
        rm "${temp_file_}"
    fi  
    
    ## PRINTING RESIDUE NAME
    echo "${res_name_}"
}

### FUNCTION TO GET TOTAL NUMBER OF ATOMS
# The purpose of this function is to find the total number of atoms 
# in an ITP file. This code assumes you have a [ bonds ] in your itp file 
# which could be used as an end point
# INPUTS:
#   $1: itp file
#   [ atoms ]
#;   nr    type   resnr  residue    atom    cgnr   charge    mass
#   1      SG311  1        DOD       S1    1  -0.082   32.0600
#   2      CG321  1        DOD       C2    1  -0.099   12.0110
# USAGE:
#   itp_find_total_atoms "charmm36-jul2017.ff/dodecanethiol.itp"
function itp_find_total_atoms () {

    ## DEFINING INPUTS
    input_itp_file_="$1"

    ## DEFINING TEMPORARY FILE
    temp_file_="${input_itp_file_%.itp}_temp_file.itp"
    ## GETTING ITP WITHOUT COMMENTS
    grep -o '^[^;]*' ${input_itp_file_} > "${temp_file_}"

    ## FINDING ATOM
    atom_linenum=$(grep -nE '\[ atoms \]' "${temp_file_}" | sed 's/\([0-9]*\).*/\1/')
    ## ADDING 1 TO REMOVE ATOM LINE
    atom_linenum=$((${atom_linenum} + 1))
    ## DEFINING COUNT
    count=0
    
    ## LOOPING THROUGH EACH LINE
    while read line; do
        ## SEEING IF THE LINE STARTS WITH [ ]
        if [[ ${line} == \[* ]]; then
            ## PRINTING NUMBER OF ATOMS
            echo "${count}"
            break
            
        fi
        count=$((${count}+1))
        
    done <<< "$(tail -n +${atom_linenum} "${temp_file_}")"
    ## REMOVING TEMP FILE
    if [ -e "${temp_file_}" ]; then
        rm "${temp_file_}"
    fi  
    
}

### FUNCTION TO FIX POSITION RESTRAINTS
# The purpose of this function is to fix the position restraints for a single residue. The assumption is that your residue is on top of the list of indices (i.e. it starts at 1). This is useful for fixing restraints for a single residue but avoiding the issue where genrestr creates an index for all possible residues. 
# INPUTS:
#   $1: posre index file
#   $2: itp_file
# OUTPUTS:
#   updated position restraint file
# USAGE:
#   fix_posre_for_single_residue lig_posre.itp "charmm36-jul2017.ff/dodecanethiol.itp"
function fix_posre_for_single_residue () {
    ## DEFINING INPUTS
    posre_index_="${1}"
    itp_file_="${2}"
    
    ## DEFINING TEMP FILE
    temp_file_="${posre_index_%.itp}_temp_file.itp"
    ## DEFINING DEFAULT COMMENT LINE
    comment_line=";  i funct       fcx        fcy        fcz"
    
    ## ADDING HEADER
    ## GETTING POSITION RESTRAINT LINE
    posre_line=$(grep -nE '\[ position_restraints \]' "${posre_index_}" | sed 's/\([0-9]*\).*/\1/')
    ## ADDING HEADER
    head -n ${posre_line} "${posre_index_}" > "${temp_file_}"
    echo "${comment_line}" >> "${temp_file_}"
    
    ## ADDING ONE TO POSRE
    posre_line=$((${posre_line}+1))
    
    ## GETTING WITHOUT COMMENTS
    # grep -o '^[^;]*' ${posre_index_} > "${temp_file_}"
    
    ## GETTING TOTAL NUMBER OF ATOMS
    num_atoms_=$(itp_find_total_atoms ${itp_file_})
    
    ## LOOPING THROUGH EACH LINE
    while read line; do
        ## GETTING ATOM INDEX
        atom_index=$(awk '{print $1}' <<< ${line})
    
        ## SEEING IF THE LINE STARTS WITH [ ]
        if [[ ${atom_index} -le "${num_atoms_}" ]]; then
            ## PRINTING NUMBER OF ATOMS
            echo "${line}" >> "${temp_file_}"
        else
            ## STOP THE SCRIPT
            break
        fi
        
    done <<< "$(tail -n +${posre_line} "${posre_index_}" | grep -o '^[^;]*')"
    ## REMOVING TEMP FILE
    if [ -e "${temp_file_}" ]; then
        ## COPYING OVER THE TEMP FILE
        cp -r "${temp_file_}" "${posre_index_}"
        rm "${temp_file_}"
    fi 
}

### FUNCTION TO EXPAND THE BOX
# INPUTS:
    # $1: GRO File
    # $2: Amount to expand in nm's (assuming cubic/rectangular)
# OUTPUTS:
    # Expanded gro file
function gro_expand_box () {
    # The purpose of this function is that it takes the gro file and expands it by a given value
    # Input files
    gro_file="$1"
    expansionValue="$2" # Amount expanded in nms
    ## SEEING IF THERE IS ANY EXPANSION
    if (( $(awk 'BEGIN {print ("'${expansionValue}'" != "'0'")}') )); then
        # Get box vector and splitting into components
        output=$(tail -n 1 $gro_file)

        # create as array, can now access each component separately:
        array=($output)
        box_x=${array[0]}
        box_y=${array[1]}
        box_z=${array[2]}

        # Now, adding to your x, y, and z:
        new_x=$(awk "BEGIN {print $box_x + $expansionValue; exit}")
        new_y=$(awk "BEGIN {print $box_y + $expansionValue; exit}")
        new_z=$(awk "BEGIN {print $box_z + $expansionValue; exit}")

        # Using editconf to change the file (no center, box)
        gmx editconf -f "$gro_file" -o "$gro_file" -box $new_x $new_y $new_z -noc &>/dev/null

        # Printing what you did
        echo "--- Running gro_expand_box function ---"
        echo "Changing the following gro file: ${gro_file}"
        echo "Changing x box size from ${box_x} to ${new_x}"
        echo "Changing y box size from ${box_y} to ${new_y}"
        echo "Changing z box size from ${box_z} to ${new_z}"
    else
        echo "Since the expansion was set to: ${expansionValue}, we are not doing anything!"
        echo "No expansion of the gro file was made!"
    fi
}

### FUNCTION TO JOIN ARRAY TO A STRING
# USAGE 1: join_array_to_string , "${data[@]}"
# USAGE 2: rdf_xvg_input=$(join_array_to_string , "${output_file_array[@]}")
function join_array_to_string () {
  local IFS="$1"
  shift
  echo "$*"
}


### FUNCTION TO EXTRACT VARIABLE NAME
# The purpose of this function is to extract the variable name for a given script.
# INPUTS:
#   $1: file name
#   $2: string that you need to find the variable for. Note that you should have a "=" sign at the end
# OUTPUTS:
#   variable_output: [str]
#       variable output from the bash script
# USAGE:
#   extract_variable_name_from_file "extract_deep_cnn_separated_instances.sh" 'instance_pickle_name=' 
function extract_variable_name_from_file () {
    ## DEFINING INPUTS
    file_input_="$1"
    str_to_look_for_="$2"

    ## PRINTING THE OUTPUT, REMOVES ALL DOUBLE QUOTES
    variable_output=$(grep "${str_to_look_for_}" "${file_input_}" | awk -F'[= ]' '{print $2}' | sed 's/"//g' )
    ## PRINTING
    echo "${variable_output}"
}

### FUNCTION TO CHECK XTC TIME
# This function checks the trajectory time using the *.chk file.
# INPUTS:
#		$1: trajectory time variable
#		$2: prefix for the files, e.g. sam_prod for sam_prod.xtc
# OUTPUTS:
#	xtc time 
# USAGE:
#	check_xtc_time xtctime sam_prod
function check_xtc_time () {
	## INPUTS
	local __xtctime="${1:xtc_time}"
	prefix_="${2-sam_prod}"

	## DEFINING LOCAL FILE
	temp_file_="${prefix_}_temp.txt"

	## USING GMX CHECK
	gmx check -f "${prefix_}".cpt 2> "${temp_file_}"

	## GETTING LAST FRAME
	last_frame=$(grep "Last frame" "${temp_file_}" | awk '{print $NF}')

	## SAVING
	eval $__xtctime="${last_frame}"

	## REMOVING TEMP FILE
	if [ -e "${temp_file_}" ]; then
		rm "${temp_file_}"
	fi
}




#!/bin/bash

# nanoparticle_functions.sh
# This contains all global functions for ligand builder
# Created by: Alex K. Chew (12/15/2017)

### FUNCTIONS AVAILABLE TO YOU
## print_script_name: prints name of the script
## extract_lig_name_from_dir: extracts ligand name from directory name
## check_file_exist: checks if a file exists (True/False output)
## solvate_with_correct_total_atoms: Solvate a gro file that is not water
## check_output_dir: Check if output directory exists -- if so, delete
## check_ff_suffix: Check force field suffix to output something you can use in a directory name
## itp_turn_pos_off: Turns off itp file position restraints given itp file
## topology_turn_on_off: Turns on / off itp files within topology
## checkNMoveFile: Checks if file exists, then moves it accordingly
## extract_itp_resname: Extracts residue name given ITP file
## extract_itp_multiple_resname: Extracts multiple resnames
## itp_fix_genrestr: Fixes GMX GENRESTR position restraints by removing extraneous atoms
### MDP FILES
## mdp_remove_freezegrps: removes freeze groups from mdp files
## mdp_check_freeze_grps: Check if the mdp file freeze groups is correct

### OTHER FUNCTIONS
# find_residue_name_from_ligand_txt: finds residue name from ligand
#   index_read_list_: reads index list as an array
#   intersection_two_arrays: finds intersection between two arrays
#   find_outputname: finds output name
### GRO FUNCTIONS

### MISCELLANEOUS FUNCTIONS
##  str_find_col_num: finds the column number of a string list
#   check_setup_file: checks input setup files
#   get_residue_num_from_sum_file: gets residue numbers from summary file

## get_hostname_details: Gets hostname details, walltimes, etc. 

## NOMENCLATURE
#   get_output_name_gnp_water: output name for GNP in pure water
#   get_output_name_for_mixed_solvents: gets outputname for mixed-solvent systems
#   make_ndx_gold_with_ligand: makes index file with gold and ligands only
#   get_output_name_for_mixed_solvents_multiple: output name for mixed_solvent_environments
#   extract_output_name_for_mixed_solvents_multiple: extraction for mixed solvent multiple

## TRAJECTORY FUNCTIONS
#   trunc_: truncates the trajectory
#   create_xtc_centered_: creates a truncated trajectory with the gold centered


### USEFUL GREP / SED / AWK COMMANDS
# REMOVES EVERYTHING FOR GREP EXCEPT LINE NUMBER: sed 's/\([0-9]*\).*/\1/'
# GREP PRINT LINE AND LINE NUMBER: grep -nE
# AWK PRINT EACH COLUMN: awk '{print $1}' <-- prints first column

## PRINTING
echo "*** LOADING NANOPARTICLE FUNCTIONS (nanoparticle_functions.sh) ***"

#### FUNCTIONS ####

## FUNCTION TO PRINT SCRIPT NAME (No inputs)
# USAGE: print_script_name
function print_script_name () {
    script_name=`basename "$0"`
    echo "*** RUNNING SCRIPT: ${script_name} ***"
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

### FUNCTION TO COUNT TOTAL NUMBER OF ATOMS AND NUMBER OF RESIDUES (ASSUME HOMOGENOUS)
# DESIGNED TO FIX GMX SOLVATION ERROR
# INPUTS:
#   $1: GRO FILE (Assumed single molecule / homogenous)
#   $2: TOP FILE (Will edit this file)
#   $3: OUTPUT GRO FILE
# USAGE: read_total_atoms butanethiol.gro x butanethiol_solv.gro
function solvate_with_correct_total_atoms () {
    echo "---- solvate_with_correct_total_atoms ----"
    echo "This script was designed to fix gmx solvate error with topology inputs"
    ## DEFINING FILES
    gro_file="$1"
    top_file="$2"
    output_gro="$3"

    ## FINDING RESIDUE NAME
    residue_name=$(sed "\$d" ${gro_file} | sed "1,2d" | head -n1 | awk '{print $1}' | sed 's/[0-9]//g')
    echo "RESIDUE NAME: ${residue_name}"
    
    ## FINDING TOTAL ATOMS IN THE NEW SOLVATED GRO FILE
    solv_total_atoms=$(sed "\$d" ${output_gro} | sed "1,2d" | grep "${residue_name}" | wc -l )
    echo "TOTAL ATOMS IN SOLVATED: ${solv_total_atoms}"
    
    ## FINDING TOTAL ATOMS IN ONE SINGLE RESIDUE
    # PRE-FOR LOOP
    current_count=0; first_try="false"
    
    # LOOPING THROUGH SOLVATED GRO FILE TO FIND THE NUMBER OF ATOMS IN 1 RESIDUES
    while read line; 
    do 
        # Look through each line, find the second column, and remove all non-numerical values
        current_atom_value=$(echo "$line" | awk '{print $2}' | sed 's/[^0-9]//g')
        echo "Current atom number: ${current_atom_value}"
        # Assuming you have more than 1 atom
        if [[ ${current_atom_value} -eq "1" ]]; then 
            if [[ ${first_try} == "true" ]]; then
                echo "Found the total atoms as: ${current_count}"
                solv_total_atoms_in_res=$(echo ${current_count})
                break
            else
                # Stating that you already tried the first one
                first_try="true"
                echo "Turning off first try"
            fi
        fi
        # Adding 1 to current count
        current_count=$(( $current_count+1 ))
    
    done <<< "$(sed "\$d" "${output_gro}" | sed "1,2d" )" # <-- Reads gro and delete 1st two lines, and last line
    
    ## FINDING TOTAL RESIDUES TO ADD TO TOPOLOGY
    total_res=$(awk -v total_atoms=$solv_total_atoms -v atoms_in_res=$solv_total_atoms_in_res 'BEGIN{ printf "%d", total_atoms/atoms_in_res }')
    echo "TOTAL RESIDUES IN SOLVATION: ${total_res}"
    
    ## CHANGING TOP FILE
    # REMOVING LAST LINE WITH THE RESIDUE NAME
    Line_of_residue=$(grep -nE "${residue_name}" "${top_file}" | sed 's/\([0-9]*\).*/\1/')
    current_res_in_top=$(grep -E "${residue_name}" "${top_file}" | awk '{print $2}')
    # REPLACING THE NUMBERS
    sed -i "${Line_of_residue}s/${current_res_in_top}/${total_res}/" "${top_file}"
    echo "REPLACING TOPOLOGY RESIDUE OF: ${current_res_in_top}"
    echo "NEW RESIDUE NUMBER: ${total_res}"
}

### FUNCTION TO COUNT TOTAL NUMBER OF ATOMS IN AN ITP FILE
# This function counts total atoms from the [ atoms ] directories
# $1: itp file (no itp extension)
# USAGE: itp_count_atoms ROT_NS.itp
function itp_count_atoms () {
    # DEFINING VARIABLES
    itp_file="$1.itp"
    
    # Using a placeholder for atom numbers
    atom_num="0"
    
    # LOOPING THROUGH AND COUNTING
    while read line;
    do
        current_atom_num=$(echo "$line" | awk '{print $1}' | sed 's/[^0-9]//g')
        
        ## CHECKING IF ATOM NUMBER IS VALID
        if [[ "${current_atom_num}" =~ ^-?[0-9]+$ ]]; then
            atom_num=$(( $atom_num+1 ))
        fi

    done <<< "$(sed '/^;/ d' ${itp_file} | sed -e '1,/\[ atoms \]/d' | sed '/\[ bonds \]/q')" # Removes comments, looks for atoms, then looks for bonds (i.e. truncate file into atoms to bonds)
    
    echo "${atom_num}"
}


### FUNCTION TO CHECK IF OUTPUT DIRECTORY EXISTS -- IF SO, DELETE
# $1: OUTPUT DIRECTORY
function check_output_dir () {
    output_dir="$1"
if [ -e $output_dir ];
    then
        echo "Removing duplicate output directory: $output_dir"
        echo "Deleting..... pause 5 seconds in case you want to cancel"
        sleep 5
        rm -rfv "$output_dir"
fi

mkdir -p ${output_dir}
}

### FUNCTION TO CHECK THE FORCE FIELD SUFFIX
# INPUTS:
#   $1: forcefield: Name of force field
# OUTPUTS:
#   forcefield_suffix: Suffix of the force field
# USAGE: forcefield_suffice=$(check_ff_suffix opls-aa)
function check_ff_suffix () {
    forcefield="$1"
    # Creating names for force fields
    if [[ "${forcefield}" == "opls-aa" ]]; then
        forcefield_suffix="OPLS"
    elif [[ "${forcefield}" == "charmm36-nov2016.ff" ]]; then
        forcefield_suffix="CHARMM36"
    elif [[ "${forcefield}" == "charmm36-jul2017.ff" ]]; then
        forcefield_suffix="CHARMM36jul2017"
    elif [[ "${forcefield}" == "pool.ff" ]]; then
        forcefield_suffix="POOL"
    else
        echo "Error! Check forcefield_suffix in nanoparticle_functions.sh"
        exit
    fi
        echo "${forcefield_suffix}"
}

### FUNCTION TO TURN OFF POSITION RESTRAINTS GIVEN ITP FILE
# INPUTS
#   $1: ITP FILE
# USAGE: itp_turn_pos_off butanethiol.itp
function itp_turn_pos_off () {
    # DEFINING ITP FILE
    itp_file="$1"

    # FINDING LINE NUMBER OF POSITION RESTRAINTS
    line_of_pos_rest=$(grep -nE "position_restraints" ${itp_file} | sed 's/\([0-9]*\).*/\1/')
    echo $line_of_pos_rest
    
    # LOOPING TO FIND ALL THE LINES OF POSITION RESTRAINTS
    if [ -z "${line_of_pos_rest}" ]; then
        echo "No position restraints found"
    else
        echo "LINE OF POSITION RESTRAINT: ${line_of_pos_rest}";
        
        # DECLARING ARRAY
        declare -a pos_res_array=("${line_of_pos_rest}")
        counter=0; previous_value="${line_of_pos_rest[@]}";
        while read line;
        do
            if [ -z "$line" ]; then
                echo "Found empty line, stopping turning off itp file"
                break
            fi
            declare -a pos_res_array=("${pos_res_array[@]}" "$(( ${previous_value} + 1 ))")
            counter=$(( ${counter}+1 ))
        done <<< $(tail -n +$(( ${line_of_pos_rest}+1 )) ${itp_file})

        # PRINTING NUMBERS
        echo "${pos_res_array[@]}"
    fi
    
    # LOOPING AND TURNING OFF POSITION RESTRAINTS
    for current_array_value in "${pos_res_array[@]}"; do
        echo "Turning off position restraints in line ${current_array_value} in $1"
        ## CHECKING IF THE POSITION RESTRAINTS ARE ALREADY COMMENTED OUT
        if [[ $(sed -n "${current_array_value} p" $1) != \;* ]]; then
            sed -i "${current_array_value}s/^/; /" "$1"
        fi
    done
}

### FUNCTION TO TURN ITP FILES ON OR OFF IN TOPOLOGY FILE
# $1: ON for on, OFF for off
# $2: TOPOLOGY FILE
# $3: Within "" to turn off
# USAGE:  topology_turn_on_off ON sam.top ROT_SN.itp
function topology_turn_on_off () { 
    ## DEFINING VARIABLES
    turn_on_off_option="$1"
    topology_file="$2"
    string_to_turn_off="$3"
    
    ## GETTING LINE NUMBER INSTANCE
    line_num=$(grep -nE "${string_to_turn_off}" ${topology_file} | sed 's/\([0-9]*\).*/\1/')
    
    # IF TURNING OFF
    if [[ ${turn_on_off_option} == "OFF" ]]; then
        ## COMMENTING OUT
        # CHECKING IF ALREADY COMMENTED OFF
        if [[ $(sed -n "${line_num} p" ${topology_file}) != \;* ]]; then
            sed -i "${line_num}s/^/; /" ${topology_file}
            echo "Turning OFF ${string_to_turn_off} in ${topology_file}"
        else
            echo "${string_to_turn_off} is already OFF in ${topology_file} -- not doing anything"
        fi
    # TURNING ON
    else
        if [[ $(sed -n "${line_num} p" ${topology_file}) == \;* ]]; then
            sed -i "${line_num}s/^; *//" ${topology_file}
            echo "Turning ON ${string_to_turn_off} in ${topology_file}"
        else
            echo "${string_to_turn_off} is already ON in ${topology_file} -- not doing anything"
        fi
    fi
}


### FUNCTION TO CHECK IF FILE EXISTS -- THEN MOVES IT
# $1: File of interest
# $2: Location you want to move it to if it exists
function checkNMoveFile () {
if [ -e $1 ]; then
    echo "$1 exists! Moving to $2"
    mv $1 $2
else
    echo "$1 does not exist! Continuing without moving."
fi
}

### FUNCTION TO EXTRACT NAME OF RESIDUE GIVEN ITP FILE
# This function looks into your itp file and extracts the residue name from it
# $1: ITP file
# Usage: extract_itp_resname itp_file_name
function extract_itp_resname () {
    itp_file="${1%.itp}"
    # Defining itp file
    itp_file="${itp_file}.itp"

    # Finding residue name
    resname=$(sed '/^;/ d' "${itp_file}" | grep -E -A 1 '\[ moleculetype \]' | tail -n1 | awk '{print $1}') # Removes comments, finds molecule type, looks at second line, then prints first column
    
    # Printing to variable
    echo "${resname}"
}

### FUNCTION TO EXTRACT RESIDUE NAMES FOR MULTIPLE LIGAND NAMES
# This function looks into ligand_names and finds the residue names for each ligand name
## INPUTS:
#   ARGS: ligand names separated by spaces, e.g. C11NH3p C11NH2
## OUTPUTS:
#   array, separated by spaces of residue name, e.g. NH3 NH2
function extract_itp_multiple_resname () {
    ## DEFINING ALL ITP FILE NAMES
    ARGS=("$@")
    
    ## DEFINING RESIDUE NAME LIST
    residue_name_list=()
    
    ## LOOPING THROUGH ARGUMENTS AND FINDING RESIDUE NAME
    for each_arg in ${ARGS[@]}; do
        ## FINDING RESIDUE NAME
        res_name="$(extract_itp_resname ${each_arg})"
        
        ## STORING TO AN ARRAY
        residue_name_list+=(${res_name})
    done
    ## PRINTING
    echo "${residue_name_list[@]}"
}



### FUNCTION TO EDIT THE ITP FILE TO REMOVE ALL EXTRANEOUS ATOMS FROM GENRESTR
# Given total number of atoms for a particular ligand, this script will look into your position restraint file and remove all extra atoms. 
# $1: ITP file
# $2: total number of atoms in a single ligand
# USAGE: itp_fix_genrestr posre.itp 84
function itp_fix_genrestr () {
    ## DEFINING VARIABLES
    itp_file="$1"
    total_atoms="$2"
    
    ## PRINTING
    echo "--- itp_fix_genrestr ---"
    echo "Fixing ITP file for ${itp_file}, given total atoms: ${total_atoms}"
    
    ## USING FOR LOOP TO FIND THE LINE THAT MATCHES
    while read line; 
    do 
        # Look through each line, find the second column, and remove all non-numerical values
        current_atom_value=$(echo "$line" | awk '{print $1}' | sed 's/[^0-9]//g')
        # Assuming you have more than 1 atom
        if [[ ${current_atom_value} -eq "${total_atoms}" ]]; then 
            lineOfInterest="$line"
            echo "Line of interest: ${lineOfInterest}"
            break
        fi
    
    done <<< "$(sed '/^;/ d' ${itp_file} | sed -e '1,/\[ position_restraints \]/ d' )" # <-- Reads itp file, removes all comments, looks for everything following position restraints
    
    ## FINDING LINE NUMBER TO REMOVE
    line_num_ending=$(grep -nE "${lineOfInterest}" ${itp_file}  | head -n1 | sed 's/\([0-9]*\).*/\1/')
    echo -e "Removing lines ${line_num_ending} onward in ${itp_file}\n"
    
    ## USING SED TO REMOVE ALL LINES FOLLOWING
    sed -i "$(( ${line_num_ending} + 1 ))"',$d' ${itp_file}
    # head -n${line_num_ending} ${itp_file} > ${itp_file}
    echo $'' >> ${itp_file} # ADDING BLANK LINE
}


### FUNCTION TO LOOK AT THE SERVER AND OUTPUT THE DETAILS WOULD LIKE
## INPUTS: None
## OUTPUTS:
#   ServerKeyWord: Keyword of the server
#   walltime: Maximum walltime in the server
#   numberOfCPUs: Number of CPUs desired

# USAGE: read ServerKeyWord walltime numberOfCPUs <<< $(get_hostname_details)
function get_hostname_details () {

    # IFS="|" # Separation for divider
    # Changing work folder for different servers
    if [[ $(hostname) == *"stampede"* ]]; then # Check if we are in the stampede server
        ServerKeyWord="STAMPEDE" # Server name (Used to find the folder with submission scripts)
        walltime="48:00:00" # Change WALLTIME in days-hh:mm:ss
        numberOfCPUs=32 # Change NUMBEROFCPUS

    elif [[ $(hostname) == *"comet"* ]]; then
        ServerKeyWord="COMET" # Server name (Used to find the folder with submission scripts)
        walltime="48:00:00" # Change WALLTIME in days-hh:mm:ss
        numberOfCPUs=24 # Change NUMBEROFCPUS

    elif [[ $(hostname) == *"swarm"* ]]; then
        ServerKeyWord="SWARM" # Server name (Used to find the folder with submission scripts)
        walltime="7-00:00:00" # Change WALLTIME in days-hh:mm:ss
        numberOfCPUs=28 # Change NUMBEROFCPUS

    elif [[ $(echo $HOSTNAME) == *"chtc"* ]]; then # Check if we are at the CHTC Supercomputers
        ServerKeyWord="CHTC_HPC" # Server name (Used to find the folder with submission scripts)
        walltime="5-00:00:00" # Change WALLTIME in days-hh:mm:ss
        numberOfCPUs=20 # Change NUMBEROFCPUS
    fi
    
    ## PRINTING TO VARIABLE
    echo "${ServerKeyWord} ${walltime} ${numberOfCPUs}"
}

### FUNCTION TO FIND THE COLUMN NUMBER GIVEN A LIST OF VALUES
# The purpose of this function is that given a list of values (e.g. BUT Ion, etc.), I would like to extract the column number of the matched string
## INPUTS:
#   $1: string with a list of values
#   $2: matching pattern
## OUTPUTS:
#   Column number
function str_find_col_num () {
    # PRINTS $1, SEPARATES BY SPACE, CREATES LIST, LOOK FOR ITEM IN LIST, THEN PRINT COLUMN
    echo "$1" | sed 's/ /\n/g' | nl | grep "$2" | awk '{print $1}'
}


### FUNCTION TO REMOVE FREEZE GROUPS OF A SPECIFIC TYPE
# The purpose of this function is to remove freeze groups of specific types, e.g. Ion. This will simply comment them off.
# Note: we assume that your mdp file has "freezegrps" and "freezedim"
## INPUTS:
#   $1: mdp file name 
#   $2: Freeze group name (e.g. Ion)
## OUTPUTS: updated mdp file with the group name removed
## USAGE EXAMPLE: mdp_remove_freezegrps nvt_double_equil_gmx5_freeze_lig_charmm36.mdp Ion
function mdp_remove_freezegrps () {
    ## PRINTING
    echo "----- MDP REMOVE FREEZE GROUPS ------"
    ## LOCATING ALL LINES WITH FREEZE GROUPS
    freeze_grps_line_num=$(grep -nE "freezegrps" "$1" | grep "$2" | sed 's/\([0-9]*\).*/\1/')
    freeze_dim_line_num=$(grep -nE "freezedim" "$1" | sed 's/\([0-9]*\).*/\1/')
    ## FINDING WITHIN THE LINE IF YOUR NAME IS THERE
    if [[ ! -z "${freeze_grps_line_num}" ]]; then
        ## PRINTING
        echo "Found $2 within $1 in line ${freeze_grps_line_num}, removing..."
        echo "Fixing dimensions in line ${freeze_dim_line_num}..."
        ######## FREEZE GROUPS
        ## FINDING ALL RESIDUES (Using sed to print everything after equal signs)
        residues="$(sed -n "${freeze_grps_line_num}p" $1 | sed 's/^.*=//g')"
        # RETURNS: BUT Ion
        ## FINDING COLUMN INDEX
        res_col_index=$(str_find_col_num "${residues}" $2)
        ## FINDING FINAL COLUMN INDEX
        init_col_index=$(awk -v col_index=$res_col_index 'BEGIN{ printf "%d",(col_index-1)*3+1}')
        last_col_index=$(awk -v col_index=$init_col_index 'BEGIN{ printf "%d",(col_index+2)}')
        ## FINDING RESULTING STRING WITHOUT THE RESIDUES
        output_residues=$(echo "${residues}" | sed "s/$2//g")
        ######### FREEZE DIMENSIONS
        ## NOW FINDING STRING AFTER EQUALS FREEZZE GROUPS
        freeze_grps_dim=$(sed -n "${freeze_dim_line_num}p" $1 | sed 's/^.*=//g' )
        ## FINDING FINAL FREEZE GROUP DIMENSIONS (Using awk to print every other column)
        output_freeze_grp_dims=$(echo "${freeze_grps_dim}" | awk -v f=${init_col_index} -v t=${last_col_index} '{for(i=1;i<=NF;i++)if(i>=f&&i<=t)continue;else printf("%s%s",$i,(i!=NF)?OFS:ORS)}')
        
        ### EDITING THE MAIN FILE
        echo "CHANGING FREEZE GROUP LINE TO: freezegrps=${output_residues}"
        echo "CHANGING FREEZE GROUP DIM. TO: freezedim=${output_freeze_grp_dims}"
        sed -i "${freeze_grps_line_num}s/.*/freezegrps=${output_residues}/" $1
        sed -i "${freeze_dim_line_num}s/.*/freezedim=${output_freeze_grp_dims}/" $1
    else
        echo "$2 residue is not found as a freeze group in $1..."
        echo "No changes have been made to $1"
    fi
}

### FUNCTION TO FIND RESIDUE NAME
# The purpose of this function is to find the residue name from the ligand name text
# INPUTS:
#   $1: ligand name text
#   $2: path to ligand names
# OUTPUTS:
#   ligand residue name
# USAGE:
#   find_residue_name_from_ligand_txt butanethiol
# As a variable: var=$(find_residue_name_from_ligand_txt butanethiol)
function find_residue_name_from_ligand_txt () {
    ## DEFINING INPUTS
    input_ligand_name_="$1"
    input_ligand_text_path_="${2:-${LIGAND_NAMES_TXT}}"

    ## LOCATING RESIDUE NAME
    residue_name_=$(grep "\b${input_ligand_name_}\b" "${input_ligand_text_path_}" | awk '{print $2}')

    ## PRINTING
    echo "${residue_name_}"

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

### FUNCTION TO GET NOMENCLATURE OF THE NAME
# The purpose of this function is to find the outputname given multiple inputs
# INPUTS:
#   $1: shape type
#   $2: temperature that you want it in
#   $3: diameter in nm
#   $4: ligand names
#   $5: force field name
#   $6: trial number
# OUTPUTS:
#   output directory name
# USAGE:
#   output_dirname="$(find_outputname ${shape_type} ${temp} ${diameter} ${ligand_names} ${forcefield} ${trial_num})"
function find_outputname () {
    ## DEFINING INPUTS
    shape_type_="$1"
    temp_="$2"
    diameter_="$3"
    ligand_names_="$4"
    force_field_="$5"
    trial_num_="$6"
    
    ## FINDING FORCE FIELD SUFFIX
    ### CHECKING FORCEFIELD SUFFIX
    forcefield_suffix=$(check_ff_suffix "${force_field_}")
    
    ## DEFINING OUTPUT DIRECTORY NAME
    output_dirname="${shape_type_}_${temp_}_K_${diameter_}_nmDIAM_${ligand_names_}_${forcefield_suffix}_Trial_${trial_num_}"
    
    ## PRINTING
    echo "${output_dirname}"
}


## FUNCTION TO CHECK INTEGRITY OF SETUP FILE
# $1: Full path to setup file
function check_setup_file () {
    # DEFINING VARIABLES
    setup_file_path="$1" # file path to setup file
    molecule_name="$2"
    # PRINTING
    echo "Checking integrity of setup file........"
    echo "Setup file directory: $1"
    
    # CHECKING FILES
    check_itp="$(check_file_exist ${setup_file_path}/${molecule_name}.itp)"
    check_gro="$(check_file_exist ${setup_file_path}/${molecule_name}.gro)"
    check_prm="$(check_file_exist ${setup_file_path}/${molecule_name}.prm)"
    echo "ITP file exists: ${check_itp}"
    echo "GRO file exists: ${check_gro}"
    echo "PRM file exists: ${check_prm}"
    
    if [ "${check_itp}" == "True" ] && [ "${check_gro}" == "True" ] && [ "${check_prm}" == "True" ]; then
        echo "Integrity intact, itp, gro, and prm file is correctly located"
    else
        echo "ERROR! Check setup files. You should have an itp, gro, and prm file... stopping here."
        exit
    fi
}

####################
### NOMENCLATURE ###
####################

### FUNCTION TO GET NAME OF OUTPUT
# The purpose of this function is to get the output name for mixed solvents
# INPUTS:
#   $1: cosolvent: name of the cosolvent
#   $2: cosolvent_perc: percentage of the cosolvent
#   $3: system_temp: temperature of the mixed-solvent system
#   $4: initial_box_length: initial box length of the system
#   $5: perctype: type of percentage it is (e.g. massfrac)  
# OUTPUTS:
#   output_name: output name for mixed solvent directory
function get_output_name_for_mixed_solvents () {
    ## DEFINING INPUTS
    cosolvent_="$1"
    cosolvent_perc_="$2"
    system_temp_="$3"
    initial_box_length_="$4"
    perctype_="${5-massperc}"   
    
    ## GETTING OUTPUT NAME
    output_name="${initial_box_length_}_${cosolvent_}_${cosolvent_perc_}_${perctype_}_${system_temp_}"
    ## PRINTING
    echo "${output_name}"
}

### FUNCTION TO GET OUTPUTNAME FOR GNPS IN WATER
# The purpose of this function is to get the output name for GNP in pure water
# INPUTS:
#   $1: shape type
#   $2: temp
#   $3: diameter
#   $4: ligand names
#   $5: force field suffix
#   $6: trial nubmer
# OUTPUTS:
#   output_name: outputname for gnp in water
# USAGE EXAMPLE:
# output_dirname=$(get_output_name_gnp_water "${shape_type}" "${temp}" "${diameter}" "${ligand_names}" "${forcefield}" "${trial_num}")
function get_output_name_gnp_water (){
    ## DEFINING INPUTS
    shape_type_="$1"
    temp_="$2"
    diameter_="$3"
    ligand_names_="$4"
    forcefield_="$5"
    trial_num_="$6"
    
    ### CHECKING FORCEFIELD SUFFIX
    forcefield_suffix=$(check_ff_suffix "${forcefield_}")
    
    ## GETTING DIRECTORY NAME
    output_dirname="${shape_type_}_${temp_}_K_${diameter_}_nmDIAM_${ligand_names_}_${forcefield_suffix}_Trial_${trial_num_}"
    
    ## PRINTING
    echo "${output_dirname}"

}

### FUNCTION TO EXTRACT OUTPUT NAME FOR GNP IN WATER
# The purpose of this function is to extract the output name for mixed solvents
# INPUTS:
#   $1: input name, e.g. 
#       EAM_300.00_K_2_nmDIAM_C11CONH2_CHARMM36jul2017_Trial_1 
# OUTPUTS:
#   Returns divide details:
#       -0: temperature
#       -1: box size
#       -2: ligand names
#       -3: trial
# USAGE:
#   read -a extract_array <<< $(extract_output_name_gnp_water EAM_300.00_K_2_nmDIAM_C11CONH2_CHARMM36jul2017_Trial_1 )
function extract_output_name_gnp_water () {
    ## DEFINING INPUTS
    input_name_="$1"
    
    ## SPLITTING ARRAY
    my_array=($(echo ${input_name_} | tr "_" "\n")) 
    #"_"
    # RETURNS: 8_nm 300_K 1_mf aceticacid_formate_methylammonium_propane
    
    ## SPLITTING ARRAY FURTHER TO GET DETAILS
    # BOX SIZE
    box_size_=${my_array[3]}
    
    # TEMPERATURE
    temp_=${my_array[1]}
    
    # TEMPERATURE
    ligname_=${my_array[5]}

    # TEMPERATURE
    trial_=${my_array[-1]}

    ## DECLARING OUTPUT ARRAY
    declare -a output_array=("${temp_}" \
                             "${box_size_}" \
                             "${ligname_}" \
                             "${trial_}")
    
    ## PRINTING
    echo "${output_array[@]}"

}



### FUNCTION TO CREATE INDEX FILE
# The purpose of this function is to generate an index file of the gold with the ligand only. 
# INPUTS:
#   $1: TPR file
#   $2: index file (which will be outputted)
#   $3: ligand residue name
#   $4: gold residue name
# OUTPUTS:
#   Index file with gold, system, and together gold + ligands
function make_ndx_gold_with_ligand () {
    ## DEFINING INPUTS
    input_tpr_="$1"
    index_file_="$2"
    ligand_residue_name_="$3"
    ## DEFINING GOLD RESIDUE NAME
    gold_residue_name_="${4-AUNP}"
    
    ## SEEING IF EXISTING AND REMOVING IF TRUE
    if [ -e "${index_file_}" ]; then
        rm "${index_file_}"
    fi
    
## GMX INDEX
gmx make_ndx -f "${input_tpr_}" -o "${index_file_}" >/dev/null 2>&1 << INPUTS
keep 0
r ${gold_residue_name_}
r ${gold_residue_name_} | r ${ligand_residue_name_}
q
INPUTS
}

### FUNCTION TO GET RESIDUE NAME FROM A SUMMARY FILE
# The purpose of this function is to get the residue name from a summary file.
# INPUTS:
#   $1: residue name
#   $2: summary file
# OUTPUTS:
#   total residue number, assuming your second column has the residue number
## USAGE:
#   num_res=$(get_residue_num_from_sum_file res_num summary_file)
function get_residue_num_from_sum_file () {
    ## DEFINING INPUTS
    residue_name_="$1"
    summary_file_="$2"
    ## GETTING RESIDUE NUMBER
    num_residues_=$(grep "${residue_name_}" "${summary_file_}"| awk '{print $2}')
    ## PRINTING
    echo "${num_residues_}"
}


### FUNCTION TO LOOK UP MOST LIKELY CONFIGURATIONS 
# <-- depreciated from BASH_MOST_LIKELY_SCRIPT
# The purpose of this function is to extract a gro / pdb file
# from the most likely configuration.
## INPUTS:
#   $1: path_output_summary - output summary file
#   $2: simulation directory location
#   $3: simulation output directory location
#   $4: tpr file within sim directory
#   $5: extracted simulation trajectory
#   $6: ligand residue name
#   $7: new gro file name
#   $8: most likely index
#   $9: gold-ligand index file
#   $10: distance to the edge of the box
## OUTPUTS:
#   gro file that is extracted and gro file that is aligned
## USAGE:
#   extract_most_likely_configuration np_most_likely.summary ./ ./ sam_prod.tpr sam_prod_10_ns_whole_center.xtc DOD np_most_likely 1 gold_with_ligand.ndx 0.8
function extract_most_likely_configuration () {
    ## DEFINING INPUTS
    path_output_summary="$1" # output summary file

    ## DEFINING SIM INPUT DIRECTORY LOCATION
    path_to_sim_output_directory="$2"
    
    ## DEFINING SIM OUTPUT DIRECTORY
    path_to_new_sim_output_directory="$3"
    
    ## TPR file
    sim_tpr_file="${4:-sam_prod.tpr}"
    
    ## XTC FILE
    extract_sim_xtc_file="${5:-sam_prod_10_ns_whole_center.xtc}"
    
    ## DEFINING LIGAND RESIDUE NAME
    lig_residue_name="$6"
    
    ## DEFINING NEW GRO FILE
    new_gro_file="$7"
    
    ## DEFINING MOST LIKELY INDEX
    most_likely_index="${8:-1}" # most likely index desired
    ## DEFINING INDEX FILE
    index_file="${9:-gold_with_ligand.ndx}" # index file
    ## DEFINING DISTANCE TO EDGE
    dist_to_edge="${10:-0.8}"
    
    ### DEFAULTS
    ## DEFINING GOLD RESIDUE NAME
    gold_residue_name="AUNP"

    
    ######################################################
    ### LOCATING MOST LIKELY CONFIGURATION AND DUMPING ###
    ######################################################
    ## DEFINING LINES TO LOOK UP TO
    lines_to_look_up="$((${most_likely_index}+1))"

    ## FINDING NEAREST TRAJECTORY TIME
    most_likely_time=$(head -n "${lines_to_look_up}" "${path_output_summary}" | tail -n1 | awk '{print $4}')
    echo "The most likely time found: ${most_likely_time} ps"; sleep 1

    ## CREATING INDEX FILE
    make_ndx_gold_with_ligand "${path_to_sim_output_directory}/${sim_tpr_file}" \
                              "${path_to_sim_output_directory}/${index_file}"   \
                              "${lig_residue_name}"   \
                              "${gold_residue_name}"  

    ## DEFINING COMBINED NAME
    combined_name="${gold_residue_name}_${lig_residue_name}"

    ## DUMPING NEAREST TRAJECTORY
echo "Creating ${new_gro_file}.gro from ${path_to_sim_output_directory}..."
gmx trjconv -f "${path_to_sim_output_directory}/${extract_sim_xtc_file}" -s "${path_to_sim_output_directory}/${sim_tpr_file}" -o "${path_to_new_sim_output_directory}/${new_gro_file}.gro" -n "${path_to_sim_output_directory}/${index_file}" -dump "${most_likely_time}" -pbc mol -center >/dev/null 2>&1  << INPUTS
${gold_residue_name}
${combined_name}
INPUTS
# >/dev/null 2>&1 

echo "Creating aligned gro file: ${new_gro_file%.gro}_align.gro"
    ## EDITING SAM FILE TO ALIGN AND MAKE DISTANCE OF NP (WITH CUBIC BOX)
gmx editconf -f "${new_gro_file}.gro" -o "${new_gro_file%.gro}_align.gro" -d "${dist_to_edge}" -princ -bt cubic >/dev/null 2>&1  << INPUTS
System
INPUTS

echo "Creating aligned pdb file: ${new_gro_file%.gro}_align.pdb"
## CREATING PDB FILE
gmx editconf -f ${new_gro_file%.gro}_align.gro -o ${new_gro_file%.gro}_align.pdb >/dev/null 2>&1 

}


### FUNCTION TO GET NAME OF OUTPUT
# The purpose of this function is to get the output name for mixed solvents
# INPUTS:
#   $1: initial_box_length: initial box length of the system
#   $2: system_temp: system temperature
#   $3: cosolvent_perc: percentage of the cosolvent
#   $4: cosolvent_name_: name of the cosolvent system
# OUTPUTS:
#   output_name: output name for mixed solvent directory
function get_output_name_for_mixed_solvents_multiple () {
    ## DEFINING INPUTS
    initial_box_length_="$1"
    system_temp="$2"
    cosolvent_perc_="$3"
    cosolvent_name_="$4"
    
    ## GETTING OUTPUT NAME
    output_name="${initial_box_length_}_nm-${system_temp}_K-${cosolvent_perc_}_mf-${cosolvent_name_}"
    ## PRINTING
    echo "${output_name}"
}

### FUNCTION TO EXTRACT OUTPUT NAME FOR MIXED SOLVENTS
# The purpose of this function is to extract the output name for mixed solvents
# INPUTS:
#   $1: input name, e.g. 
#       8_nm-300_K-1_mf-aceticacid_formate_methylammonium_propane
# OUTPUTS:
#   Returns divide details:
#       - box size
#       - temperature
#       - mole fraction
#       - solvents
# USAGE:
#   read -a extract_array <<< $(extract_output_name_for_mixed_solvents_multiple ${name})
function extract_output_name_for_mixed_solvents_multiple (){
    ## DEFINING INPUTS
    input_name_="$1"
    
    ## SPLITTING ARRAY
    my_array=($(echo ${input_name_} | tr "-" "\n")) #"_"
    # RETURNS: 8_nm 300_K 1_mf aceticacid_formate_methylammonium_propane
    
    ## SPLITTING ARRAY FURTHER TO GET DETAILS
    # BOX SIZE
    box_size_=($(echo ${my_array[0]} | tr "_" "\n"))
    
    # TEMPERATURE
    temp_=($(echo ${my_array[1]} | tr "_" "\n"))
    
    # MOLE FRAC
    molfrac_=($(echo ${my_array[2]} | tr "_" "\n"))
    
    # SOLVENTS
    solvents_=$(echo ${my_array[3]} | sed "s/_/,/g" )
    
    ## DECLARING OUTPUT ARRAY
    declare -a output_array=("${box_size_[0]}" \
                             "${temp_[0]}" \
                             "${molfrac_[0]}" \
                             "${solvents_}")
    
    ## PRINTING
    echo "${output_array[@]}"

}




### FUNCTION TO UPDATE MOLECULAR SOLVENTS
# The purpose of this function is to update the molecular details of solvents
# INPUTS:
#   $1: solvent name
# OUTPUTS:
#   - molecular details should be updated on a list
function update_molecular_solvents () {
    ## SOURCING FILES
    source "${BASH_EXTRACT_SOLVENT_NP_DEFAULTS}"
    ## DEFINING INPUTS
    input_solvent_name_="$1"
    
    ## RUNNING SCRIPT
    bash "${BASH_EXTRACT_SOLVENT_DETAILS}" "${path_to_dir}"  \
                                           "${input_solvent_name_}" \
                                           "${equil_suffix}" \
                                           "${equil_folder}" \
                                           "${initial_frame}" \
                                           "${output_list_name}" \
                                           "${rewrite}"
    
}

### FUNCTION TO GET MOLECULAR VOLUME OF A COSOLVENT NAME
# The purpose of this function is to get molecular volume of a solvent
# INPUTS:
#   $1: solvent name
# OUTPUTS:
#   output molecular volume
function get_molecular_volume_solvent () {
    ## SOURCING FILES
    source "${BASH_EXTRACT_SOLVENT_NP_DEFAULTS}"
    ## DEFINING INPUTS
    input_solvent_name_="$1"
    
    ## GETTING MOLECULAR VOLUME
    molecular_volume="$(grep "^${input_solvent_name_}" ${path_free_volume_list} | awk '{print $3}' )"
    
    ## CHECKING IF MOLECULAR VOLUME IS THERE
    if [ -z "${molecular_volume}" ]; then
        update_molecular_solvents "${input_solvent_name_}"
        
    ## GETTING MOLECULAR VOLUME
    molecular_volume="$(grep "^${input_solvent_name_}" ${path_free_volume_list} | awk '{print $3}' )"
    
    fi
    
    ## EXTRACTING MOLECULAR VOLUME
    echo "${molecular_volume}"

}

### FUNCTION TO FIND THE NEAREST SOLVENT WITH THE BOX LENGTH
# The purpose of this function is to find the best mixed solvent for 
# your desired system. 
## INPUTS:
#   $1: multiple solvent name, separated by commas
#   $2: mole fraction of each solvent
#   $3: box length desired for the solvent
#   $4: output text file
## OUTPUTS:
#   output file with details of solvents for best directory
function get_multiple_solvents_best_dir () {
    ## DEFINING INPUTS
    solvent_name_="$1"
    #"aceticacid,formate,methylammonium,propane"
    mole_frac_="$2"
    #"1"
    ## DEFINING DESIRED BOX LENGTH
    box_length_desired_="$3"
    #"7"
    
    ## DEFINING OUTPUT TEXTBOX
    output_details="$4"
    # "output.txt"
    
    ## READING TO ARRAY
    declare -a box_lengths_avail_=($( grep "${solvent_name_} ${mole_frac_}" "${MULTIPLE_SOLVENT_INFO_TXT}" | awk '{print $NF}'))
    declare -a dir_names=($( grep "${solvent_name_} ${mole_frac_}" "${MULTIPLE_SOLVENT_INFO_TXT}" | awk '{print $5}'))
    
    ## PRINTING
    echo "; Extraction of details for ${solvent_name_}, ${mole_frac_} mole frac" > "${output_details}"
    echo -e "; Desired box length: ${box_length_desired_} \n" >> "${output_details}"
    
    ## DEFINING BOX SIZE AVAILABLE
    if [ ! -z "${box_lengths_avail_}" ]; then
        ## PRINTING
        echo "; Solvent name, box size, difference value" >> "${output_details}"
        ## COMPUTING DIFFERENCE IN THE ARRAY
        len=${#box_lengths_avail_[*]}
        difference_array=()
        for (( i=0; i<=$(( $len -1 )); i++ ))
        do
            ## COMPUTING DIFFERENCEVALUE
            diff_value=$(echo "scale=4;${box_lengths_avail_[$i]}-${box_length_desired_}" | bc)
            ## COMPUTING DIFFERENCE ARRAY
            difference_array+=(${diff_value})
            ## ADDING TO LIST
            echo ""${dir_names[$i]}" "${box_lengths_avail_[i]}" "${diff_value}"" >> "${output_details}"
            
        done
        
        echo "" >> "${output_details}"
        echo "; Finding best solvent" >> "${output_details}"
        
        ## DEFINING SMALLEST VALUE
        smallest_value=10000000
        
        ## LOOPING THROUGH DIFFERENCE ARRAY AND FINDING THE SMALLEST POSITIVE VALUE
        for (( i=0; i < ${#difference_array[@]}; ++i )); do
            ## DIFFERENCE VALUE
            diff_value=${difference_array[i]}
            ## UPDATING SMALLEST VALUE
            # if [[ "${diff_value}" < "${smallest_value}" ]] && [[ "${diff_value}" > 0 ]]; then
            if (( $(echo "${diff_value} < ${smallest_value}" | bc -l) && $(echo "${diff_value} > 0" | bc -l) )); then
                smallest_value=${diff_value}
                index=$i
            fi
        done
        
        ## CHECKING IF IT WORKS
        if [[ "${smallest_value}" == 10000000 ]]; then
            echo "ERROR: No smallest value found" >> "${output_details}"
        else
            ## GETTING THE SMALLEST ARRAY VALUE
            solvent_name=$( grep "${solvent_name_} ${mole_frac_}" "${MULTIPLE_SOLVENT_INFO_TXT}" | head -$((${index}+1)) | tail -1 | awk '{print $5}')
            ## PRINTING
            echo "BEST_SOLVENT: ${solvent_name}" >> "${output_details}"
            echo "BEST_DISTANCE: ${smallest_value}" >> "${output_details}"
        
        fi  
    
    else
        echo "ERROR: No solvent system found for ${solvent_name_} and ${mole_frac_} mol frac" >> "${output_details}"
    fi

}

### FUNCTION TO USE TRJCONV TO SHORTEN THE TIME
# INPUTS:
#   $1: input tpr file
#   $2: input xtc file
#   $3: system name
#   $4: output xtc file
#   $5: trunction time in ps (lower bound)
function trunc_ {
    ## DEFINING VARIABLES
    input_tpr_file="$1"
    input_xtc_file="$2"
    system_name="$3"
    output_xtc_file="$4"
    truncation_time="$5"
    
gmx trjconv -s "${input_tpr_file}" -f "${input_xtc_file}" -o "${output_xtc_file}" -b "${truncation_time}" -pbc whole << INPUT
${system_name}
INPUT
}


### FUNCTION TO CREATE A TRAJECTORY THAT HAS A MOLECULE OF INTEREST CENTERED
# The purpose of this function is to get trajectories with the molecules centered
# INPUTS:
#   $1: input tpr file
#   $2: input xtc file
#   $3: centering group
#   $4: output_xtc file
#   $5: output ndx file
# OUTPUTS:
#   xtc file with the molecule of interest centered
# USAGE: create_xtc_centered_ tpr_file xtc_file AUNP output_xtc_file
function create_xtc_centered_ {
    ## DEFINING INPUTS
    input_tpr_file_="$1"
    input_xtc_file_="$2"    
    atom_selection_="$3"
    system_selection="System"
    ## DEFINING OUTPUTS
    output_xtc_file_="$4"
    output_ndx_file_="$5"

### CREATING SPECIFIC INDEX FILE
gmx make_ndx -f "${input_tpr_file_}" -o "${output_ndx_file_}" << INPUTS
keep 0
keep 1
r ${atom_selection_}
q
INPUTS


## RUNNING TRJCONV
gmx trjconv -f "${input_xtc_file_}" -s "${input_tpr_file_}" -o "${output_xtc_file_}" -n "${output_ndx_file_}" -center -pbc mol << INPUTS
${atom_selection_}
${system_selection}
INPUTS

}

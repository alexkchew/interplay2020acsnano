#!/bin/bash

# gridrc.sh
# The purpose of this script is to load all variables for running 
# grid functions. 

# Written by: Alex K. Chew (03/29/2020)

## DEFINING IMPORT PYTHON FUNCTIONS

## GETTING CURRENT LOCATION
dir_loc="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

## DEFINING PATH TO MODULES
PATH_TO_MODULES="$(cd ${dir_loc}/../modules; pwd)"

## DEFINING PYTHON PATH
export PYTHONPATH="${PYTHONPATH}:${PATH_TO_MODULES}"

## DEFINING MDDESCRIPTORS
MDDESCRIPTORS="${PATH_TO_MODULES}/MDDescriptors"

## DEFINING PATH TO GRID
PATH_TO_GRID_BASH="${dir_loc}/generate_grid.sh"

## DEFINING GRIDDING FUNCTION
PYTHON_SCRIPTS_GRIDDING="${MDDESCRIPTORS}/surface/generate_wc_grid.py"
PYTHON_SCRIPTS_MU="${MDDESCRIPTORS}/surface/combine_neighbors_array.py"
PYTHON_SCRIPTS_NUM_NEIGHBORS="${MDDESCRIPTORS}/surface/generate_hydration_maps_parallel.py"

PYTHON_SCRIPTS_REMOVE_GRID_FOR_PLANAR="${MDDESCRIPTORS}/application/np_hydrophobicity/remove_grids_for_planar_SAMs.py"

#################
### FUNCTIONS ###
#################

function gro_measure_box_size () 
{ 
    input_gro_file_="$1";
    output=$(tail -n 1 ${input_gro_file_});
    array=($output);
    echo "${array[@]}"
}

### FUNCTION TO JOIN ARRAY TO A STRING
# USAGE 1: join_array_to_string , "${data[@]}"
# USAGE 2: rdf_xvg_input=$(join_array_to_string , "${output_file_array[@]}")
function join_array_to_string () {
  local IFS="$1"
  shift
  echo "$*"
}

### FUNCTION TO CREATE DIRECTORIES
create_dir () 
{ 
    directory="$1";
    dir_exists=$(check_file_exist ${directory});
    if [ "${dir_exists}" == "True" ]; then
        if [ "$2" != "-f" ]; then
            echo "${directory} already exists! Do you want to delete and recreate? (y/n)";
            read deletion_criteria;
            while [ "${deletion_criteria}" != "y" ] && [ "${deletion_criteria}" != "n" ]; do
                echo "Error! Incorrect prompt! Deletion criteria can only be \"y\" or \"n\", try again:";
                read deletion_criteria;
            done;
        else
            deletion_criteria="y";
        fi;
        if [ "${deletion_criteria}" == "y" ]; then
            echo "Deleting and recreating ${directory}, pausing 3 seconds...";
            sleep 3;
            rm -rv "${directory}";
            mkdir -p "${directory}";
        else
            if [ "${deletion_criteria}" == "n" ]; then
                echo "Stopping here -- deletion prevented";
                echo "Check if you need these files. This error message is to prevent overwriting of data";
                exit;
            fi;
        fi;
    else
        echo "Creating ${directory}";
        mkdir -p "${directory}";
    fi
}

check_file_exist () 
{ 
    if [ -e "$1" ]; then
        echo "True";
    else
        echo "False";
    fi
}

stop_if_does_not_exist () 
{ 
    echo "Checking if file exists: $1 ...";
    if [ -e "$1" ]; then
        echo "--> $1 exists, continuing!";
    else
        echo "XXX $1 does not exist!";
        echo "Check for errors! You may have a non-existant, vital file!";
        echo "Pausing here for 5 seconds so you can see this error! ...";
        sleep 5;
        exit 1;
    fi
}


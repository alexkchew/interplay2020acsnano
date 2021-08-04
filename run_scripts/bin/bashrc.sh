#!/bin/bash

# bashrc.sh
# This contains general codes and functions for np_hydrophobicity_project

## DEFINING MAIN DIRECTORY
dir_loc="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

## DEFINING MAIN DIRECTORY
main_dir="${dir_loc}/../../"

## DEFINING CURRENT WORKING DIRECTORY
export CURRENTWD=$(cd ${main_dir}; pwd)

## DEFINING SCRIPTS DIRECTORY
export PATH2SCRIPTS="${CURRENTWD}/run_scripts"

## DEFINING BIN DIRECTORY
export PATH2BIN="${PATH2SCRIPTS}/bin"

## DEFINING MAJOR FILES TO LOAD
### LOADING GLOBAL VARIABLES
source "${PATH2BIN}/global_vars.sh"
### LOADING GLOBAL GENERAL FUNCTIONS
source "${PATH2BIN}/functions.sh"
## LOADING GENERAL FUNCTIONS
source "${PATH2BIN}/server_general_research_functions.sh"

## NANOPARTICLE FUNCTIONS
source "${PATH2BIN}/nanoparticle_functions.sh"

## NANOPARTICLE FUNCTIONS
source "${PATH2BIN}/nanoparticle_functions_shared.sh"

### LOADING GLOBAL FUNCTIONS ACROSS ALL PROJECTS <-- will decide if necessary
# source "${HOME}/bin/bashfiles/server_general_research_functions.sh"

## EXPORTING PYTHON PATHS
export PYTHONPATH="${HOME}/bin/pythonfiles/modules:$PYTHONPATH"
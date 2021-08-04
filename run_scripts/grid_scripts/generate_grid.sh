#!/bin/bash

# generate_grid.sh
# The purpose of this script is to generate a grid for a NP system

# Written by: Alex K. Chew (10/24/2019)

## VARIABLES:
#   $1: path_sims -- path to simulations
#   $2: grid size
#   $3: output prefix
#   $4: input TPR file
#   $5: input XTC file
#   $6: n_procs -- number of processors
#   $7: analysis directory
#   $8: rewrite: true/false
#   $9: input_prefix for gro traj, etc.
#   $10: last gridding frame to compute, default is 1000 frames

################################
### LOADING GLOBAL VARIABLES ###
################################
dir_loc="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
source "${dir_loc}/gridrc.sh"

#########################
### INITIAL VARIABLES ###
#########################

## DEFINING WORKING DIRECTORY
path_sims="$1"

## DEFINING GRID SIZE
mesh="${2:-0.1,0.1,0.1}"
# "0.1"

## OUTPUT PREFIX
output_prefix="${3:-out}"
# "out"

## GRO AND XTC
input_tpr_file="${4:-sam_prod.tpr}"
input_xtc_file="${5:-sam_prod.xtc}"

## DEFINING NUMBER OF CORES
n_procs="${6:-20}"

## DEFINING ANALYSIS DIRECTORY
analysis_dir="${7:-hyd_analysis}"

## DEFINING REWRITE
rewrite=${8:-false}

## DEFINING INPUT PREFIX
input_prefix="${9:-sam_prod}"

## DEFINING FIRST AND LAST FRAME
begin_traj="${14-0}" 
end_traj=${10:-"1000"} # 1 ns required for gridding

## DEFINNG ALPHA AND CONTOUR LEVEL
alpha="${11:-None}"
contour="${12-None}"

## DEFINING PYTHON
python_func="${13-python3.6}"

## DEFINING IF NORMALIZATION IS DESIRED
want_normalize_c="${13-false}"

## DEFINING OUTPUT PREFIX
output_traj_prefix="${input_prefix}-${begin_traj}_${end_traj}-watO_grid"

## DEFINING OUTPUT GRO, XTC, AND TPR FILES
output_xtc_file="${output_traj_prefix}.xtc"
output_gro_file="${output_traj_prefix}.gro"
output_tpr_file="${output_traj_prefix}.tpr"

## DEFINING OUTPUT NAME
output_dir="${analysis_dir}/grid-${begin_traj}_${end_traj}" # -${mesh}

#########################
### DEFAULT VARAIBLES ###
#########################

## DEFINING WATER AND NOT HYDROGEN
water_oxygen_only_input="r SOL & a O*"

## DEFINING OUTPUT INDEX NAME
water_oxygen_only_output="SOL_&_O*"

## DEFINING INDEX FILE
index_file="${output_traj_prefix}.ndx"

#################
### MAIN CODE ###
#################

## CHECKING IF EXISTS
stop_if_does_not_exist "${path_sims}"

## GOING TO PATH
cd "${path_sims}"

## MAKING INDEX
if [[ ! -e "${index_file}" ]] || [[ "${rewrite}" == true ]]; then
gmx make_ndx -f "${input_tpr_file}" -o "${index_file}" << INPUTS
keep 0
keep 1
${water_oxygen_only_input}
q
INPUTS

fi

## USING TRJCONV TO GET TRUNCATED SYSTEM
if [[ ! -e "${output_xtc_file}" ]] || [[ "${rewrite}" == true ]]; then
## XTC FILE
gmx trjconv -f "${input_xtc_file}" \
            -s "${input_tpr_file}" \
            -o "${output_xtc_file}" \
            -pbc mol \
            -b "${begin_traj}" \
            -e "${end_traj}" \
            -n "${index_file}" << INPUTS
${water_oxygen_only_input}
INPUTS

## TPR FILE
gmx convert-tpr -s "${input_tpr_file}" \
                -o "${output_tpr_file}" \
                -n "${index_file}" << INPUTS
${water_oxygen_only_input}
INPUTS

## DUMPING GRO FILE
gmx trjconv -f "${output_xtc_file}" \
            -s "${output_tpr_file}" \
            -o "${output_gro_file}" \
            -dump "${begin_traj}" << INPUTS
${water_oxygen_only_input}
INPUTS

fi

## REMOVING ANY EXTRAS HASHTAGS
rm -f \#*


## CREATING DIRECTORY
if [ ! -e "${output_dir}" ]; then
    if [[ "${rewrite}" == true ]]; then
        create_dir "${output_dir}" -f
    else
        mkdir -p "${output_dir}"
    fi
fi

## CHECKING IF EXISTING
if [[ ! -e "${path_sims}/${output_dir}/${output_prefix}_willard_chandler.dat" ]] || [[ "${rewrite}" == true ]]; then

## RUNNING PYTHON CODE
${python_func} "${PYTHON_SCRIPTS_GRIDDING}" --path "${path_sims}" \
                                       --gro "${output_gro_file}" \
                                       --xtc "${output_xtc_file}" \
                                       --output_prefix "${output_prefix}" \
                                       --mesh "${mesh}" \
                                       --n_procs "${n_procs}" \
                                       --output_file "${output_dir}" \
                                       --alpha "${alpha}" \
                                       --contour "${contour}" \
                                       --want_normalize_c "${want_normalize_c}" \
                                       --debug

else
    echo "Since ${output_prefix}_willard_chandler.dat exists, continuing!"
fi

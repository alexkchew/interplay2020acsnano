#!/bin/bash

# extend_hydrophobicity_sims.sh
# The purpose of this script is to simply extend the simulations for simulations that are already complete. We will extend based on your desired extension period.

# Written by: Alex K. Chew (03/05/2020)


################################
### LOADING GLOBAL VARIABLES ###
################################
dir_loc="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
source "${dir_loc}/bin/bashrc.sh"

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

#########################
### INITIAL VARIABLES ###
#########################

## TOTAL TIME
time_desired="150000.000"
# "100000.000" # ps for production

## DEFINING SIMULATION LOCATION
sim_folder_name="20200224-planar_SAM_spr50"
# "20200224-GNP_spr_50"
# 
# "20200224-GNP_spr_50"
# "$1"

## DEFINING PATH TO SIM
path_sims="${PATH_SIMULATIONS}/${sim_folder_name}"

## checking path
stop_if_does_not_exist "${path_sims}"

## DEFINING NUMBER OF CORES
num_cores="28"

## DEFINING PATH TO SCRIPTS
input_run_script="run_extend_job.sh"
input_submit_script="submit_run_extend_job.sh"

## PATHS TO IT
path_input_run_script="${PATH_BASH_RUN_SCRIPTS}/${input_run_script}"
path_input_submit_script="${PATH_BASH_RUN_SCRIPTS}/${input_submit_script}"

## DEFINING OUTPUT PREFIX
output_prefix="sam_prod"

## DEFINING LIST OF DIRECTORIES
read -a sim_list <<< $(ls ${path_sims}/* -d)

## LOOPING THROUGH SIMS
for each_sim in ${sim_list[@]}; do
	## GO TO DIRECTORY
	cd "${each_sim}"
	## FINDING BASENAME
	sim_name="$(basename ${each_sim})"
	## PRINTING
	echo "Working on ${sim_name}"
	## CHECKING THE TIME
	check_xtc_time current_xtc_time "${output_prefix}"
	if [[ "${time_desired}" != "${current_xtc_time}" ]]; then
		## DEFINING JOB NAME
		job_name="${sim_name}_extended_${time_desired}"
		## COPYING SUBMIT FILE
		cp -r "${path_input_run_script}" ./"${input_run_script}"
		cp -r "${path_input_submit_script}" ./"${input_submit_script}"
		## EDITING SUBMISSION FILE
		sed -i "s/_NUMCORES_/${num_cores}/g" "${input_submit_script}"
		sed -i "s/_OUTPUTPREFIX_/${output_prefix}/g" "${input_submit_script}"
		sed -i "s/_TOTALPRODTIME_/${time_desired}/g" "${input_submit_script}"
		sed -i "s/_RUNCODE_/${input_run_script}/g" "${input_submit_script}"
		sed -i "s/_JOB_NAME_/${job_name}/g" "${input_submit_script}"
		## ADDING SUBMISSION INTO JOB LIST
		echo "${each_sim}/${input_submit_script}" >> "${JOB_LIST}"

	fi

done

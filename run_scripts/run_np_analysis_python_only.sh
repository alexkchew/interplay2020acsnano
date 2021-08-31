#!/bin/bash

# run_np_analysis_python_only.sh
# The purpose of this script is to run the nanoparticle analysis.
# This is an updated code with no gmx select
# This code is written purely in python. 
# 
# Written by: Alex K. Chew (01/03/2020)

# DEFINING VARIABLES
#   $1: simulation folder name
#   $2: simulation folder within $1
#   $3: contour level
#   $4: selection_type:
#           water_and_counterions: water + counterions
#           water_only: water only
#           all_heavy: all heavy atoms
#   $10: grid proc

# MAKE INDEX FILE
# r SOL & ! a H*
# OUTPUT: SOL_&_!H*
# resname SOL CL NA and not name \"H.*\
# r SOL & a O*
#   r SOL & a H*
#   !6 & !2
################################
### LOADING GLOBAL VARIABLES ###
################################
dir_loc="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
source "${dir_loc}/bin/bashrc.sh"

### FUNCTION TO GET THE NAME OF OUTPUT FILE
# The purpose of this function is to get the output file 
# of the hydration map name
# INPUTS:
#   $1: contour
#   $2: sigma
#   $3: mesh
#   $4: cutoff radius
#   $5: residue list that you care about
# OUTPUTS:
#   output name in the form of contour-sigma-mesh
# USAGE:
#   output_file=$(get_hydration_map_analysis_output_name ${contour} ${sigma} ${mesh})
function get_hydration_map_analysis_output_name () {
    ## DEFINING INPUTS
    contour_="$1"
    sigma_="$2"
    mesh_="$3"
    cutoff_radius_="$4"
    residue_list_="$5"
    begin_traj="$6"
    end_traj="$7"
    want_contour_c="${8:-false}"
    ## DEFINING WC INTERFACE TIME
    wc_grid_initial="$9"
    wc_grid_final="${10}"

    
    ## CONVERTING MESH TO STRING
    mesh_underscores=$(sed "s/,/_/g" <<< ${mesh_} )
    
    ## CONVERTING RESIDUE LIST
    residue_list_underscores=$(sed "s/,/_/g" <<< ${residue_list_} )
    
    if [[ "${want_contour_c}" == true ]]; then
        prefix="norm-"
    else
        prefix=""
    fi

    ## GETTING OUTPUT NAME
    output_name_="${prefix}${contour_}-${sigma_}-${mesh_}-${cutoff_radius_}-${residue_list_underscores}-${begin_traj}-${end_traj}-wc_${wc_grid_initial}_${wc_grid_final}"
    
    ## PRINTING OUTPUT NAME
    echo "${output_name_}"
}

#########################
### INITIAL VARIABLES ###
#########################

## DEFINING FORCEFIELD
forcefield="${FORCEFIELD}"


## DEFINING LOCATION
sim_folder_name="$1"
# "PLANAR"
# "$1"
# "NP_SPHERICAL"
# "191010-Most_likely_np_sims_NVT_charged_groups"

## DEFINING SPECIFIC FILE
sim_file_name="$2"
# "FrozenGoldPlanar_300.00_K_dodecanethiol_10x10_CHARMM36jul2017_intffGold_Trial_1-50000_ps"
# "$2"
# "MostlikelynpNVT_EAM_300.00_K_2_nmDIAM_C11OH_CHARMM36jul2017_Trial_1_likelyindex_1"
# "MostlikelynpNVT_EAM_300.00_K_2_nmDIAM_C11NH2_CHARMM36jul2017_Trial_1_likelyindex_1"

##############################
### USER DEFINED VARIABLES ###
##############################



## OUTPUT PREFIX FOR WILLARD CHANDLER
output_prefix="out"

## DEFINING GRID SIZE
mesh="0.1,0.1,0.1"

## DEFINING REWRITE OPTIONS
rewrite=false
# true
# true
# false
# true
# false
# true if rewrite desired
# false if no rewriting

## DEFINING REWRITE WC INTERFACE
rewrite_wc_interface=false

## DEFINING CONTOUR
contour=${3-"25.6"}  # 80% of the bulk (bulk ~ 32)
# contour="16"    # 50% of the bulk

## DEFINING SIGMA
sigma="0.24"

## DEFINING TOTAL TIME TO TRUNCATE FROM BEGINNING
truncate_time="2000" # ps

## DEFINING IF YOU WANT THE COUNTERIONS
selection_type="${4:-all_heavy}"

## DEFINING BEGINNING AND END TRAJ
begin_traj="2000"
# "0"   # 0 ns
end_traj="${5-50000}"
# "50000" # 50 ns

## DEFINING IF YOU WANT WC INTERFACE ONLY
want_wc_only="${6-false}"

## DEFINING CUTOFF RADIUS
cutoff_radius="${7-0.33}"

## DEFINING WANT CONTOUR LEVEL
want_contour_c="${8-false}"

## DEFINING PURE WATER SIMULATIONS
pure_water_sim="${9-false}"

## DEFINING DETAILS FOR PROCESSING
grid_n_procs="${10-20}"

## DEFINING PROCS FOR PYTHON CODE
n_procs="${11-1}"

## FINDING LIGAND NAME
if [[ "${pure_water_sim}" == false ]]; then
    lig_name=$(extract_lig_name_from_dir "${sim_file_name}")

    ## DEFINING TYPE
    initial_text=$(cut -d'_' -f"1" <<< "${sim_file_name}")

    ## GETTING TEXT
    if [[ "${initial_text}" == "FrozenPlanar"* ]] || \
        [[ "${initial_text}" == "FrozenGoldPlanar"* ]] || \
        [[ "${sim_file_name}" == *"Planar"* ]]; then
        ## DEFINING TYPE
        job_type="planar"
        ## DEFINING PLANAR SAM
        planar_sam="True"
    else
        job_type="spherical"
        planar_sam="False"
    fi
else
    job_type="water_only"
fi

## DEFINING RESIDUE LIST BASED ON selection type
if [[ "${selection_type}" == "water_and_counterions" ]]; then
    residue_list="HOH,CL,NA"
elif [[ "${selection_type}" == "water_only" ]]; then
    residue_list='HOH'
elif [[ "${selection_type}" == "all_heavy" ]]; then
    residue_list='all_heavy'
fi

## PRINTING
echo "Residue list: ${residue_list}"

## DEFINING FRAME RATE
frame_rate="28"
# "5000"
# "10000"
# "None"
# "5000"
# "None"

#########################
### DEFAULT VARAIBLES ###
#########################
## PYTHON SCRIPT
python_script_mu="${MDDESCRIPTORS}/surface/combine_neighbors_array.py"
python_script_num_neighbors="${MDDESCRIPTORS}/surface/generate_hydration_maps_parallel.py"

## SCRIPT LOCATION
scripts_folder="scripts"
module_folder=${scripts_folder}/modules

## DEFINING GRID BASHRC
gridrc="gridrc.sh"

## DEFINING RELATIVE GRIC RC
relative_gridrc_path="${scripts_folder}/$(basename ${PATH_GRID_SCRIPTS})/${gridrc}"

## GRID LOCATION
grid_folder="wc_grid"

## COMBINING PICKLE ARRAYS
python_script_combine_pickles="${MDDESCRIPTORS}/surface/combine_neighbors_array.py"

## DEFINING IF YOU WANT FROM THE END
want_grid_from_end=true
# true

## GETTING END TRAJ GRID
begin_traj_grid="0"
end_traj_grid="5000"
# "10000"
# "1000"

## GRO AND XTC
if [[ "${job_type}" != "water_only" ]]; then
    input_prefix="sam_prod"
    ## DEFINING LIGAND ITP FILE
    ligand_itp_file="${lig_name}.itp"
    echo "Ligand itp file: ${ligand_itp_file}"

else
    input_prefix="water_prod"
fi

## DEFINING INPUT FILES
input_gro_file="${input_prefix}.gro"
input_tpr_file="${input_prefix}.tpr"
input_xtc_file="${input_prefix}.xtc"


## DEFINING INDEX FILE
index_file="no_hydrogens.ndx"
# "gold_ligand.ndx"

###############
### SCRIPTS ###
###############

## BASH SCRIPT
bash_hydration_map_code="${PATH2SCRIPTS}/extract_hydration_maps_with_python.sh"
## DEFINING GRID CODE
bash_grid_code="${PATH2SCRIPTS}/generate_grid.sh"

## INPUT SUBMISSION SCRIPT
input_submission_script="${PATH_SUBMISSIONS}/submit_hydration_maps.sh"
output_submission_script="submit.sh"

## DEFINING SLURM 
slurm_output="hydration.out"

######################
### DEFINING PATHS ###
######################
## PATH TO SIMS
path_sims="${PATH_SIMULATIONS}/${sim_folder_name}/${sim_file_name}"

## CHECKING IF EXISTS
stop_if_does_not_exist "${path_sims}"

## GOING INTO DIRECTORY
cd "${path_sims}"

### GETTING GRID INFO ###

## CHECKING IF XTC IS CORRECT
check_xtc_time current_xtc_time "${input_prefix}"

## CHECKING THE TIME
if (( $(echo "${end_traj} > ${current_xtc_time}" | bc -l) )); then
    echo "--------------------------------------------"
    echo "Warning! Trajectory specified (${end_traj} ps) is greater than available trajectory of ${current_xtc_time}"
    echo "If you do not have the correct trajectory, it may skew your analysis!"
    echo "Stopping here:"
    echo "${path_sims}"
    echo "Please check your trajectory time!"
    echo "Pausing for 3 seconds"
    echo "--------------------------------------------"
    sleep 3
    exit
fi

## GETTING TRAJECTORY TIMES
if [[ "${want_grid_from_end}" == true ]]; then
    echo "Computing WC interface from the end"
    ##  COMPUTING INITIAL TRAJ
    begin_traj_grid=$(awk -v current_time=${end_traj} \
                       -v end_traj=${end_traj_grid} 'BEGIN{ printf "%d", current_time - end_traj }')

    ## FINDING END TRAJ
    end_traj_grid=$(awk -v current_time=${end_traj} \
                        'BEGIN{ printf "%d", current_time }')
    echo "New traj times: ${begin_traj_grid} - ${end_traj_grid} ps"
fi


## GETTING OUTPUT FILE NAME
if [[ "${job_type}" != "water_only" ]]; then
    output_file=$(get_hydration_map_analysis_output_name ${contour} ${sigma} ${mesh} ${cutoff_radius} ${residue_list} ${begin_traj} ${end_traj} ${want_contour_c} ${begin_traj_grid} ${end_traj_grid})
else
    output_file="${sigma}-${cutoff_radius}-${begin_traj}-${end_traj}"
fi


## CREATING FILES
mkdir -p "${module_folder}"
# mkdir -p "${grid_folder}"

## COPYING OVER SCRIPTS
cp -r "${MDDESCRIPTORS}" "${module_folder}"
cp -r "${MDBUILDERS}" "${module_folder}"

## COPYING GRID SCRIPTS
cp -r "${PATH_GRID_SCRIPTS}" "${scripts_folder}"


## CREATING DIRECTORY
if [ ! -e "${output_file}" ]; then
    mkdir -p "${output_file}"

else
    if [[ "${rewrite}" == true ]]; then
        create_dir "${output_file}" -f
    fi
fi


## DEPRECIATED -- NO LONGER NEED LIGAND RESIDUE NAMES
# ## DEFINING GOLD RESIDUE NAME
# if [[ "${job_type}" == "spherical" ]]; then
#     gold_residue_name="AUNP"
#     lig_residue_name=$(itp_get_resname ${ligand_itp_file})
# elif [[ "${job_type}" == "planar" ]]; then
#     gold_residue_name="AUI"
#     lig_residue_name=$(itp_get_resname ${forcefield}/${ligand_itp_file})
# elif [[ "${job_type}" == "water_only" ]]; then
#     gold_residue_name=""
#     lig_residue_name=""
# fi

# echo "${job_type}"
# ## CHECKINGI RESIDUE NAME FOUND
# if [[ -z "${lig_residue_name}" ]] && [[ "${job_type}" != "water_only" ]]; then
#     echo "Error! Missing ligand residue name"
#     echo "Check itp file: ${lig_residue_name}"
#     sleep 5
#     exit
# fi

## CREATING INDEX FILES
echo "Creating index file: ${index_file}"
gmx make_ndx -f "${path_sims}/${input_tpr_file}" \
             -o "${path_sims}/${index_file}" >/dev/null 2>&1 << INPUTS
keep 0
! a H*
name 1 no_hydrogens
q
INPUTS

## MAKING SURE INDEX FILE EXISTS
stop_if_does_not_exist "${path_sims}/${index_file}"

# ## CHECKING IF INDEX FILE ALREADY EXISTS
# if [[ ! -e "${path_sims}/${index_file}" ]] || [[ "${rewrite}" == true ]]; then

# fi

## REMOVING EXTRAS
rm -f "\#*"

## DEFINING COMBINED NAME
combined_name="no_hydrogens"

## DEFINING OUPTUT SIMS
path_sims_output="${path_sims}/${output_file}"

## GOING INTO ANALYSIS
cd "${path_sims_output}"

## COPYING OVER FILE
cp -r "${bash_hydration_map_code}" "${path_sims_output}"

## GETTING BASENAME
bash_hydration_file_name="$(basename ${bash_hydration_map_code})"

## USING SED TO EDIT DETAILS FOR BASH SCRIPT
sed -i "s/_MESH_/${mesh}/g" "${bash_hydration_file_name}"
sed -i "s/_NPROCS_/${n_procs}/g" "${bash_hydration_file_name}"
sed -i "s/_NGRIDPROCS_/${grid_n_procs}/g" "${bash_hydration_file_name}"
sed -i "s/_WANTPLANAR_/${planar_sam}/g" "${bash_hydration_file_name}"
sed -i "s/_GROFILE_/${input_gro_file}/g" "${bash_hydration_file_name}"
sed -i "s/_XTCFILE_/${input_xtc_file}/g" "${bash_hydration_file_name}"
sed -i "s/_TPRFILE_/${input_tpr_file}/g" "${bash_hydration_file_name}"
sed -i "s#_PYTHONSCRIPTMU_#${python_script_mu}#g" "${bash_hydration_file_name}"
sed -i "s#_OUTPUTFILE_#${output_file}#g" "${bash_hydration_file_name}"
sed -i "s#_OUTPUTPREFIX_#${output_prefix}#g" "${bash_hydration_file_name}"
## EDITING BASH GRID
sed -i "s#_BASHGRID_#${bash_grid_code}#g" "${bash_hydration_file_name}"
sed -i "s#_BEGINTRAJGRID_#${begin_traj_grid}#g" "${bash_hydration_file_name}"
sed -i "s#_ENDGRIDTRAJ_#${end_traj_grid}#g" "${bash_hydration_file_name}"
sed -i "s#_PATHTOGRIDRC_#${relative_gridrc_path}#g" "${bash_hydration_file_name}"
sed -i "s#_GRIDDIR_#${grid_folder}#g" "${bash_hydration_file_name}"
sed -i "s#_REWRITEWCINTERFACE_#${rewrite_wc_interface}#g" "${bash_hydration_file_name}"

## DEFINING REWRITE
sed -i "s#_REWRITE_#${rewrite}#g" "${bash_hydration_file_name}"
sed -i "s#_INDEXFILE_#${index_file}#g" "${bash_hydration_file_name}"
sed -i "s#_COMBINEDNAME_#${combined_name}#g" "${bash_hydration_file_name}"

sed -i "s#_BEGINTRAJ_#${begin_traj}#g" "${bash_hydration_file_name}"
sed -i "s#_ENDTRAJ_#${end_traj}#g" "${bash_hydration_file_name}"

## DEFINING ANALYSIS DIRECTORY
sed -i "s#_ANALYSISDIR_#${output_file}#g" "${bash_hydration_file_name}"
## GETTING NEIGHBORS CODE
sed -i "s#_PYTHGETNEIGHBORS_#${python_script_num_neighbors}#g" "${bash_hydration_file_name}"
## RESIDUE LIST
sed -i "s#_RESIDUELIST_#${residue_list}#g" "${bash_hydration_file_name}"
sed -i "s#_CUTOFFRADIUS_#${cutoff_radius}#g" "${bash_hydration_file_name}"
sed -i "s#_FRAMERATE_#${frame_rate}#g" "${bash_hydration_file_name}"

## CONTOUR LEVEL AND ALPHA
sed -i "s#_CONTOUR_#${contour}#g" "${bash_hydration_file_name}"
sed -i "s#_SIGMA_#${sigma}#g" "${bash_hydration_file_name}"

## WC INTEFACE ONLY
sed -i "s#_WANTWCONLY_#${want_wc_only}#g" "${bash_hydration_file_name}"

## NORMC
sed -i "s#_WANTNORMC_#${want_contour_c}#g" "${bash_hydration_file_name}"

## PURE WATER
sed -i "s#_WANTPUREWATER_#${pure_water_sim}#g" "${bash_hydration_file_name}"

####################################
### GENERATING SUBMISSION SCRIPT ###
####################################
## DEFINING JOB NAME
job_name="${sim_file_name}_${output_file}"

## DEFINING OUTPUT SUBMISSION FILE
output_submit_path="${path_sims_output}/${output_submission_script}"

## COPYING SUBMISSION FILE
cp -r "${input_submission_script}" "${output_submit_path}"

## EDITING SUBMISSION FILE
sed -i "s#_USER_#${USER}#g" "${output_submit_path}"
sed -i "s#_JOBNAME_#${job_name}#g" "${output_submit_path}"
sed -i "s#_BASHSCRIPT_#${bash_hydration_file_name}#g" "${output_submit_path}"
sed -i "s#_SLURMOUT_#${slurm_output}#g" "${output_submit_path}"

## ADDING TO JOB LIST
echo "${output_submit_path}" >> "${JOB_LIST}"


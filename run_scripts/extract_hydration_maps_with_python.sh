#!/bin/bash

# extract_hydration_maps_with_python.sh
# This script runs the python code to extract hydration maps. 
# Here, we will use the path of the python script to correctly run the nanoparticle analysis.

## DEFINING VARIABLES
# SCRIPTS
#   _PYTHONSCRIPTMU_ <-- python script to compute MU
#   _PYTHGETNEIGHBORS_ <-- python script to get neighbors
#   _BASHGRID_ <-- grid code

# FILE DETAILS
#   _GROFILE_ <-- gro path
#   _XTCFILE_ <-- xtc path
#   _TPRFILE_ <-- tpr file

## PROCESSORS
#   _NGRIDPROCS_ <-- number of procs for gridding
#   _NPROCS_ <-- number of processors
#   _MESH_ <-- mesh size
#   _OUTPUTFILE_ <-- output file
#   _OUTPUTPREFIX_ <-- output prefix
#   _WANTPLANAR_ <-- True if you want planar

#   _INDEXFILE_ <-- index file
#   _REWRITE_ <-- rewrite true/false
#   _COMBINEDNAME_ <-- combined name
#   _BEGINTRAJ_ <-- begining trajectory
#   _ENDTRAJ_ <-- ending trajectory
#   _ANALYSISDIR_ <-- analysis directory
#   _ENDGRIDTRAJ_ <-- end grid trajectory
#   _RESIDUELIST_ <-- residue list
#   _CUTOFFRADIUS_ <-- cutoff radius for searching
#   _SIGMA_ <-- sigma value for vibrations
#   _CONTOUR_ <-- contour level desired

##############
### INPUTS ###
##############

## EXPORTING PYTHON PATH
export PYTHONPATH="${HOME}/bin/pythonfiles/modules:${PYTHONPATH}"

## DEFINING FILE DETAILS
path_sims="$( cd "$( dirname "${BASH_SOURCE[0]}" )/../" >/dev/null 2>&1 && pwd )"

## DEFINING PATH TO SCRIPTS
path_to_rc="${path_sims}/_PATHTOGRIDRC_"

## LOADING RC FILE
source "${path_to_rc}"

## DEFINING BASH SCRIPT LOCATION
bash_grid_code="${PATH_TO_GRID_BASH}"

## DEFINING PYTHON SCRIPT
python_script_compute_mu="${PYTHON_SCRIPTS_MU}"

## DEFINING PYTHON SCRIPT
python_script_get_neighbors="${PYTHON_SCRIPTS_NUM_NEIGHBORS}"

## DEFINING ANALYSIS DETAILS
mesh="_MESH_"

## DEFINING FRAMES
n_procs="_NPROCS_"

## DEFINING FRAMES
begin_traj="_BEGINTRAJ_"
end_traj="_ENDTRAJ_"

## DEFINING GRID END TRAJ
begin_traj_grid="_BEGINTRAJGRID_"
end_traj_grid="_ENDGRIDTRAJ_"

## PROCS FOR GRIDDING
grid_procs="_NGRIDPROCS_"

## DEFINING IF PLANAR
planar_sam="_WANTPLANAR_"

gro_file="_GROFILE_"
xtc_file="_XTCFILE_"
tpr_file="_TPRFILE_"

## DEFINING INDEXFILE
index_file="_INDEXFILE_"

## DEFINING REWRITE
rewrite="_REWRITE_"
rewrite_wc_interface="_REWRITEWCINTERFACE_"

## DEFINING LIGAND-GOLD RESNAME
combined_name="_COMBINEDNAME_"

## DEFINING OUTPUT FILE
output_file="_OUTPUTFILE_"

## DEFINING OUTPUT PREFIX
output_prefix="_OUTPUTPREFIX_"

## DEFINING GRID DIRECTORY
grid_dir="_GRIDDIR_"

## DEFINING ANLAYSIS DIRECTORY
analysis_dir="_ANALYSISDIR_"

## PICKLE NAMES
pdb_file_name="${output_prefix}_hydration.pdb"

## DEFINING INPUT PREFIX
input_prefix="${gro_file%.gro}"

## DEFINING RESIDUE LIST
residue_list="_RESIDUELIST_"

## DEFINING CUTOFF RADIUS
cutoff_radius="_CUTOFFRADIUS_"

## DEFINING NORM C
want_contour_c="_WANTNORMC_"

## DEFINING SIGMA AND CONTOUR VALUE
sigma="_SIGMA_"
contour="_CONTOUR_"

## DEFINING IF YOU ONLY WANT WC INTERFACE
want_wc_only="_WANTWCONLY_"

## DEFINING IF YOU WANT PURE WATER SIMS
want_pure_water_only="_WANTPUREWATER_"

## IF YOU WANT DEBUG ON
want_debug="false"

########################
### STEP 1: GRIDDING ###
########################

## DEFINING GRID OUTPUT
grid_output="grid-${begin_traj_grid}_${end_traj_grid}" # -${mesh}
wc_output="${output_prefix}_willard_chandler.dat"
path_grid_data="${path_sims}/${output_file}/${grid_output}/${wc_output}"

## DEFINING PRUE WATER AND DEBUG
if [[ "${want_pure_water_only}" == false ]] || [[ "${want_debug}" == true ]]; then
  ## RUNNING BASH CODE
  bash "${bash_grid_code}" "${path_sims}" \
                           "${mesh}" \
                           "${output_prefix}" \
                           "${tpr_file}" \
                           "${xtc_file}" \
                           "${grid_procs}" \
                           "${analysis_dir}" \
                           "${rewrite_wc_interface}" \
                           "${input_prefix}" \
                           "${end_traj_grid}" \
                           "${sigma}" \
                           "${contour}" \
                           "${want_contour_c}" \
                           "${begin_traj_grid}"
fi

if [[ "${want_pure_water_only}" == true ]]; then
  ## GETTING GRO FILE
  read -a box_size <<< $(gro_measure_box_size ${path_sims}/${gro_file})

  ## CREATING DIRECTORY
  mkdir -p "${path_sims}/${output_file}/${grid_output}"

  ## GETTING HALF BOX LENGTH
  declare -a grid_pt=()

  ## LOOPING AND COMPUTING HALF BOX LENGTH
  for dim_size in ${box_size[@]}; do
    half_dim=$(awk -v dim=${dim_size} 'BEGIN{ printf "%.3f", dim/2.0}')
    ## ADDING
    grid_pt=(${grid_pt[@]} ${half_dim})
  done
  echo "Pure water simulations are turned on:"
  echo "Grid point dimensions: ${grid_pt[@]}"

  ## ADDING TO GRID
  echo "# x y z" > "${path_grid_data}"
  echo "" >> "${path_grid_data}"
  str_grid_pt=$(join_array_to_string , "${grid_pt[@]}")
  echo "${str_grid_pt}" >> "${path_grid_data}"

fi

## REMOVING GRID POINTS FOR PLANAR SAMS
if [[ "${planar_sam}" == True ]]; then
  
  ## DEFINING PATH GRID OUTPUT
  path_grid_output="${path_sims}/${output_file}/${grid_output}"
  ## MAKING COPY
  copy_dat="${wc_output%.dat}_orig_copy.dat"
  cp -r "${path_grid_data}" "${path_grid_output}/${copy_dat}"

  ## REMOVING ORIGINAL .DAT
  if [ -e "${path_grid_output}/${copy_dat}" ]; then
    rm "${path_grid_data}"
  fi

  ## PATH TO PICKLE
  path_pickle="${path_sims}/${output_file}/${grid_output}/remove_grid.pickle"

  ## REMOVING GRID POINTS
  python3.6 "${PYTHON_SCRIPTS_REMOVE_GRID_FOR_PLANAR}" --path_gro "${path_sims}/${gro_file}" \
                                                       --path_to_grid "${path_grid_output}/${copy_dat}" \
                                                       --path_to_pickle "${path_pickle}" \
                                                       --path_output_grid "${path_grid_data}"



   echo "Since planar SAMs is turned on, we will remove grid points!"
   echo "Copying over grid points to ${grid_output} -> ${copy_dat}"
   echo "Overwritting grid point file"
else
    echo "NOT PLANAR SAMS"
    echo "Not removing any grid points"

fi

### SEEING IF YOU ONLY WANT WC INTERFACE
if [[ "${want_wc_only}" == false ]]; then

##########################################
### STEP 2: GMX TO CREATE TRAJECTORIES ###
##########################################

## GO TO PATH SIM
cd "${path_sims}"

## GETTING !AUNP NAME -- CHANGED NOW TO NO HYDROGENS
combined_name_not="${combined_name}"
# "!${combined_name}_&_!SOL_H"

## DEFINING OUTPUT PREFIX FOR TRAJECTORY
output_traj_prefix="${input_prefix}_${begin_traj}_${end_traj}-heavyatoms"

## DEFINING OUTPUT
output_xtc_file="${output_traj_prefix}.xtc"
# output_tpr_file="${output_traj_prefix}.tpr"
output_gro_file="${output_traj_prefix}.gro"

## USING TRJCONV TO GET TRUNCATED SYSTEM
if [[ ! -e "${output_xtc_file}" ]] || [[ "${rewrite}" == true ]]; then
## XTC FILE
gmx trjconv -f "${xtc_file}" \
            -s "${tpr_file}" \
            -o "${output_xtc_file}" \
            -pbc mol \
            -b "${begin_traj}" \
            -e "${end_traj}" \
            -n "${index_file}" << INPUTS
${combined_name_not}
INPUTS
fi

#if [[ ! -e "${output_tpr_file}" ]] || [[ "${rewrite}" == true ]]; then
### TPR FILE
#gmx convert-tpr -s "${tpr_file}" -o "${output_tpr_file}" -n "${index_file}" << INPUTS
#${combined_name_not}
#INPUTS
#fi

if [[ ! -e "${output_gro_file}" ]] || [[ "${rewrite}" == true ]]; then
## DUMPING GRO FILE
gmx trjconv -f "${xtc_file}" \
            -s "${tpr_file}" \
            -o "${output_gro_file}" \
            -n "${index_file}" \
            -dump "${begin_traj}" << INPUTS
${combined_name_not}
INPUTS
fi

## REDEFINING GRO AND XTC, SKIPPING THE TRJCONV

###################################
#### STEP 3: GET NUM OCCURANCES ###
###################################

## DEFINING FRAME RATE
frame_rate="_FRAMERATE_"

## DEFINING PICKLE LOG FILE
pickle_log="neighbors.log"

## DEFINING PATH PICKLE
path_pickle="${path_sims}/${output_file}/compute_neighbors"

## DEFINING PATH TO GRID
path_grid="${path_sims}/${analysis_dir}/${grid_output}/${wc_output}"

## DEFINING PICKLE THAT SHOULD BE OUTPUT
final_pickle_name="${begin_traj}-${end_traj}.pickle"

## SEEING IF EXISTS
if [[ ! -e "${path_pickle}/${final_pickle_name}" ]]; then

### RUNNING PYTHON SCRIPT TO GET NUM OCCURANCES
python3.6 "${python_script_get_neighbors}" --path "${path_sims}" \
                                           --gro "${output_gro_file}" \
                                           --xtc "${output_xtc_file}" \
                                           --residue_list "${residue_list}" \
                                           --path_pickle "${path_pickle}" \
                                           --path_grid "${path_grid}" \
                                           --cutoff_radius "${cutoff_radius}" \
                                           --frame_rate "${frame_rate}" \
                                           --n_procs "${n_procs}" \
                                           --pickle_log "${pickle_log}" 
                                           
                                           

fi

########################################
### STEP 4: CONGLOMERATE ALL DETAILS ###
########################################

## DEFINING INPUTS
path_gro="${path_sims}/${gro_file}"

## PATH PICKLES
path_pdb_file="${path_sims}/${analysis_dir}/${pdb_file_name}"

## DEFINING PICKLE LOG
path_pickle_log="${path_pickle}/${pickle_log}"

# RUNNING PYTHON SCRIPT
python3.6 "${python_script_compute_mu}" --path_gro "${path_gro}" \
                                        --path_pdb_file "${path_pdb_file}" \
                                        --path_pickle_log "${path_pickle_log}" \
                                        --path_grid "${path_grid}" \
                                        --path_pickle "${path_pickle}"


fi

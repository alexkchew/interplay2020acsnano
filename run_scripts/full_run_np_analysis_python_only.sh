#!/bin/bash

# full_run_np_analysis_python_only.sh
# The purpose of this script is to run full nanoparticle analysis
# This uses gmx select protocols
# Written by: Alex K. Chew & Brad C. Dallin (10/26/2019)

################################
### LOADING GLOBAL VARIABLES ###
################################
dir_loc="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
source "${dir_loc}/bin/bashrc.sh"

#########################
### DEFAULT VARIABLES ###
#########################

## BASH FILE
bash_file="${PATH2SCRIPTS}/run_np_analysis_python_only.sh"

## PREFIX
folder_prefix="MostlikelynpNVT_EAM_300.00_K_2_nmDIAM"
# folder_prefix="MostlikelynpNVT_EAM_300.00_K_6_nmDIAM"
folder_suffix="CHARMM36jul2017_Trial_1_likelyindex_1"
# folder_suffix="CHARMM36jul2017_Trial_1_likelyindex_401"
#folder_suffix="CHARMM36jul2017_Trial_1_likelyindex_200"
# folder_suffix="CHARMM36jul2017_Trial_1_likelyindex_300"

############################
### DEFINING INPUT ARRAY ###
############################

## MAIN SIMULATION FOLDER
# main_sim_folder="191121-last_index_most_likely"
# main_sim_folder="PLANAR"
#main_sim_folder="NP_SPHERICAL_REDO"
#main_sim_folder="NP_ROT"
# main_sim_folder="20200114-NP_HYDRO_FROZEN"
#main_sim_folder="20200117-Frozen_Rot_Particles_6nm"
# main_sim_folder="20200212-planar_SAM_frozen"
# main_sim_folder="20200215-GNP_spring_const"
main_sim_folder="20200224-GNP_spr_50"
main_sim_folder="20200224-planar_SAM_spr50"

## DEFINING MAIN SIM FOR DOUBLE BONDS
main_sim_folder="20200227-50_doub"

## FOR LARGER PARTICLES
main_sim_folder="20200325-4_6nm"
main_sim_folder="20200325-6nm_dod"

## FOR PURE WATER
main_sim_folder="20200515-planar_short_lig_frozen"
# "20200512-6nm_Least_likely_config"
# "20200427-2nm_Least_likely_config_unsat"
# "20200427-2nm_Least_likely_config_charge"
# "20200421-6nm_OH"

## DEFINING PROCESSORS FOR GRID
grid_n_procs="20"

## NUMBER OF PROCESSORS FOR RUNNING COUNTING ALGORITHM
n_procs="28"
# "1"
#"20"
# "1"
# "20200421-2nm_Least_likely_config"
# "20200421-branch_frozen"
# "20200419-unsaturated_frozen"
# "20200414-6nm_particles_frozen"
# "20200411-mixed_sam_frozen"
# "20200413-Planar_SAMs-5nmbox_vac_with_npt_equil_NORESTRAINT"
# 
# "20200401-Renewed_GNP_sims_with_equil_other_molecules_PEG_4"
# "PURE_WATER_SIMS"

# ## NEW PLANAR SIMS
# main_sim_folder="20200326-Planar_SAMs_with_larger_z_frozen"

# ## MORE PLANAR SIMS
# main_sim_folder="20200326-Planar_SAMs_with_larger_z_8nm_frozen"
# main_sim_folder="20200328-Planar_SAMs_new_protocol-shorterequil_spr50"
# main_sim_folder="20200401-Renewed_GNP_sims_with_equil"
# "20200403-Planar_SAMs-5nmbox_vac_with_npt_equil"
# "20200403-Planar_SAMs-5nmbox_vac_with_npt_equil_other"
# "20200401-Renewed_GNP_sims_with_equil"
# "20200403-Planar_SAMs-5nmbox_vac_with_npt_equil_dod"
# "20200403-Planar_SAMs-5nmbox_novac_nptequil_long_nospr"
# "20200403-Planar_SAMs-5nmbox_novac_nptequil_long"
# "20200224-planar_SAM_spr50_dod_only"
# "20200401-Planar_SAMs_no_temp_anneal_or_vacuum_with_npt_equil"
# "20200331-Planar_SAMs_no_temp_anneal_or_vacuum_withsprconstant"
# "20200325-4_6nm"
# "20200327-PE-most_likely"
# "20200227-50_doub"
# "20200224-GNP_spr_50"
# "20200326-Planar_SAMs_with_larger_z_frozen_with_vacuum"
# 20200328-Planar_SAMs_new_protocol-shorterequil_spr50
# "20200326-Planar_SAMs_with_larger_z_frozen_with_vacuum"

# "20200224-planar_SAM_spr50"
# "20200224-GNP_spr_50"
# "20200215-planar_SAM_frozen_600spring"

if [ "${main_sim_folder}" == "PLANAR" ]; then
    folder_prefix="FrozenGoldPlanar_300.00_K"  
elif [ "${main_sim_folder}" == "20200114-NP_HYDRO_FROZEN" ]; then
    folder_prefix="FrozenPlanar_300.00_K"
fi

if [ "${main_sim_folder}" == "20200215-GNP_spring_const" ]; then
    folder_prefix="MostlikelynpNVTspr_600-EAM_300.00_K_2_nmDIAM"
fi

## PLANAR SURFACES SUFFIX
if [ "${main_sim_folder}" == "PLANAR" ] || [ "${main_sim_folder}" == "20200114-NP_HYDRO_FROZEN" ]; then
    folder_suffix="10x10_CHARMM36jul2017_intffGold_Trial_1-50000_ps"
fi

## 6 NM SURFACES
if [ "${main_sim_folder}" == "20200117-Frozen_Rot_Particles_6nm" ]; then
    folder_prefix="MostlikelynpNVT_EAM_300.00_K_6_nmDIAM"
fi

## SPRING CONSTANTS\
if [ "${main_sim_folder}" == "20200212-planar_SAM_frozen" ]; then
    folder_prefix="NVTspr_1000_Planar_300.00_K"
#    folder_prefix="NVTspr_50_Planar_300.00_K"
    folder_suffix="10x10_CHARMM36jul2017_intffGold_Trial_1-50000_ps"
fi

if [ "${main_sim_folder}" == "20200215-planar_SAM_frozen_600spring" ]; then
    folder_prefix="NVTspr_600_Planar_300.00_K"
    folder_suffix="10x10_CHARMM36jul2017_intffGold_Trial_1-50000_ps"
fi
 
if [ "${main_sim_folder}" == "20200224-GNP_spr_50" ] || [ "${main_sim_folder}" == "20200227-50_doub" ] ; then
    folder_prefix="MostlikelynpNVTspr_50-EAM_300.00_K_2_nmDIAM"
fi

if [ "${main_sim_folder}" == "20200224-planar_SAM_spr50" ] || [ "${main_sim_folder}" == "20200403-Planar_SAMs-5nmbox_vac_with_npt_equil_other" ]; then
    folder_prefix="NVTspr_50_Planar_300.00_K"
    folder_suffix="10x10_CHARMM36jul2017_intffGold_Trial_1-5000_ps"
fi


## DEFINING LIGAND ARRAY
declare -a ligand_names=("C11CONH2")
# "C11double67OH" "dodecen-1-thiol"
# "dodecanethiol" "C11OH" "C11OH" "C11NH2" "C11CF3" "C11CONH2" "C11COOH") # "C11NH2"
#  
# "dodecanethiol" "C11OH"
# "C11NH2" "C11CF3" "C11CONH2" "C11COOH"
# "dodecanethiol" "C11OH" 
# "C11OH"
# "C11NH2"  "C11CF3" "C11CONH2" "C11COOH"
# "dodecanethiol" "C11CF3" "C11CONH2" "C11COOH" "C11OH"
# "dodecanethiol" "C11OH"
# "C11CF3" "C11CONH2" "C11NH2" "C11COOH"
# "dodecanethiol" "C11CF3" "C11CONH2" "C11NH2" "C11COOH" "C11OH
# 
# ("dodecanethiol")
# 
# "C11CF3" "C11CONH2" "C11NH2" "C11COOH"
# "C11OH"
# "dodecanethiol"
# C11OH
#  "dodecanethiol"
# "ROTPE1" "ROTPE2"
# "ROT001" "ROT002" "ROT003" "ROT005" "ROT006" "ROT007" "ROT008" "ROT009"
# "ROT004"
# 
# "C11CF3"
# "dodecanethiol" "C11OH" "C11NH3" "C11NH2" "C11NCH33" "C11COOH" "C11COO" "C11CONH2" "C11CF3" "C11CCH33"
# "C11OH" "C11NH3" "C11COO" "C11C OOH" "C11CF3" "C11CONH2" "C11NH2"
# "dodecanethiol"
#  
# "C11CCH33" "C11NCH33" 
#  "C11COOH" 
# "C11COOH"
# "C11CCH33" "C11CF3" "C11CONH2" "C11COOH" "C11NCH33" "C11NH2"
# #  "C11COOH"
# "C11COO" "C11NH3"
# 
# "dodecanethiol" "C11OH"
# "C11COO" "C11NH3"
# "dodecanethiol" "C11OH"
# "C11COO" "C11NH3"
# "dodecanethiol" "C11OH"
#  "C11COO" "C11NH3"
# "ROTPE1" "ROTPE2" "ROTPE3" "ROTPE4"

## DEFINING CONTOUR LEVELS
declare -a contour_array=("26") # "25.6" 
# "25.6"
# "0.70"
# "0.80"
# "0.70"
# "27"
# "30"
# "30.5"
# "30" "31" "32" "32.5"
#  "33"
# "27" 
# "27" "30" "32" "32.5" "32.5" 
# "27" "30" "32" "32.5" 
#  "34" "36" "40"
# "28.8"  "20.8" "16" "8"
# 28.8 -- 90% of the bulk
# 25.6 -- 80% of the bulk
# 20.8 -- 65% of the bulk
# 16 -- 50% of the bulk
# 8 -- 25% of the bulk

## DEFINING END TRAJ
if [[ "${main_sim_folder}" == "20200224-planar_SAM_spr50" ]]; then
	end_traj="150000"
elif [[ "${main_sim_folder}" == "20200224-GNP_spr_50" ]]; then
	end_traj="50000"
	# "100000"
else
	end_traj="50000"
fi
# "150000"
# "100000"
# "150000"
# "100000"
# "150000"
# "100000"
# "150000"
# "100000"
# "50000"

## DEFINING CUTOFF RADIUS
declare -a cutoff_radius=("0.33")
# "0.25"
# "0.25"
# "0.33"
# "0.25"
# "0.20" "0.16" "0.10"
# "0.25"
# "0.33"

## DEFINING IF YOU WANT ALL POSSIBLE FOLDERS
want_all_sim_folders=true
# true

## DEFINING IF YOU WANT ONLY WC INTERFACE
want_wc_only=false
# false
# false

## SEEING IF YOU WANT NORMED
want_contour_c=false
# true

## TRUE IF YOU WANT PURE WATER SIMS
if [[ ${main_sim_folder} == "PURE_WATER_SIMS" ]]; then
	pure_water_sim=true
else
	pure_water_sim=false
fi


## REDEFINING LIGAND NAMES IF ALL SIMS IS TRUE
if [[ "${want_all_sim_folders}" == true ]]; then
	## DEFINING PATH TO SIMULATION
	path_to_sim="${PATH_SIMULATIONS}/${main_sim_folder}"

	## READING ALL FOLDERS
	IFS=' ' read -r -a ligand_names <<< $(ls ${path_to_sim})

fi

## DEFINING IF COUNTER IONS DESIRED
declare -a selection_type=("all_heavy") #  "false"  "false"
#  "false"

## LOOPING
for each_ligand in "${ligand_names[@]}"; do
	## CREATING SIM FOLDER
	if [[ "${want_all_sim_folders}" == true ]]; then
		current_sim_folder=$(basename ${each_ligand})
	else
		current_sim_folder="${folder_prefix}_${each_ligand}_${folder_suffix}"
	fi

	## LOOPING THROUGH CONTOUR LEVELS
	for each_contour in "${contour_array[@]}"; do
		## LOOPING THROUGH EACH COUNTERION
		for current_selection in "${selection_type[@]}"; do
			## LOOPING THROUGH CUTOFF
			for current_cutoff in "${cutoff_radius[@]}"; do
## PRINTING
echo "Running: "${bash_file}" 
							  "${main_sim_folder}" 
							  "${current_sim_folder}" 
							  "${each_contour}" 
							  "${current_selection}"
                              "${end_traj}"
                              "${want_wc_only}"
                              "${current_cutoff}"
                              "${want_contour_c}"
                              "${pure_water_sim}"
                              "${grid_n_procs}"
                              "${n_procs}""

				sleep 1

				## RUNNING BASH SCRIPT
				bash "${bash_file}" "${main_sim_folder}" \
									"${current_sim_folder}" \
									"${each_contour}" \
									"${current_selection}" \
	                                "${end_traj}" \
	                                "${want_wc_only}" \
	                                "${current_cutoff}" \
	                                "${want_contour_c}" \
	                                "${pure_water_sim}" \
	                                "${grid_n_procs}" \
	                                "${n_procs}"
	            done
		done
	done


done




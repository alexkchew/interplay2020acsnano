#!/bin/bash

# extract_hydration_maps.sh
# This script runs the python code to extract hydration maps. 
# Here, we will use the path of the python script to correctly run the nanoparticle analysis.

## DEFINING VARIABLES
#   _PYTHONSCRIPT_ <-- python script
#   _SIMPATH_ <-- simulation path
#   _GROFILE_ <-- gro path
#   _XTCFILE_ <-- xtc path
#   _NFRAMES_ <-- number of frames in trajectory
#   _NPROCS_ <-- number of processors
#   _MESH_ <-- mesh size
#   _OUTPUTFILE_ <-- output file
#   _OUTPUTPREFIX_ <-- output prefix
#   _WANTPLANAR_ <-- True if you want planar

## DEFINING ANALYSIS DETAILS
mesh="_MESH_"

## DEFINING FRAMES
n_frames="_NFRAMES_"
n_procs="_NPROCS_"

## DEFINING IF PLANAR
planar_sam="_WANTPLANAR_"

## DEFINING FILE DETAILS
path_sims="_SIMPATH_"
gro_file="_GROFILE_"
xtc_file="_XTCFILE_"

## DEFINING PYTHON SCRIPT
python_script="_PYTHONSCRIPT_"

## DEFINING OUTPUT FILE
output_file="_OUTPUTFILE_"

## DEFINIGN OUTPUT PREFIX
output_prefix="_OUTPUTPREFIX_"

## RUNNING PYTHON SCRIPT
python3.6 "${python_script}" --path "${path_sims}" \
                           --gro "${gro_file}" \
                           --xtc "${xtc_file}" \
                           --output_prefix "${output_prefix}" \
                           --mesh "${mesh}" \
                           --n_frames "${n_frames}" \
                           --n_procs "${n_procs}" \
                           --planar "${planar_sam}" \
                           --output_file "${output_file}"
                           
                           


#!/bin/bash 
# submit_hydration_maps.sh
# This submission script is designed to submit hydration maps
## VARIABLES:
#   _USER_ <-- User to email details to
#   _JOBNAME_ <-- job name
#   _BASHSCRIPT_ <-- bash script
#   _SLURMOUT_ <-- slurm output

###SERVER_SPECIFIC_COMMANDS_START

#SBATCH -p compute
#SBATCH -t 108:00:00
#SBATCH -J _JOBNAME_
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1              # total number of mpi tasks requested
#SBATCH --mail-user=_USER_@wisc.edu
#SBATCH --mail-type=all  # email me when the job starts

## ADDING PATH VARIABLES
export PYTHONPATH="/usr/lib64/python3.6/site-packages:${PYTHONPATH}"

###SERVER_SPECIFIC_COMMANDS_END

## DEFINING THE BASH SCRIPT
bash_script="_BASHSCRIPT_"
slurm_out="_SLURMOUT_"

## RUNNING BASH SCRIPT
time bash "${bash_script}"
# > "${slurm_out}" 2>&1


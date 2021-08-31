#!/bin/bash 
# submit_run_extend_job.sh
# This script extends the job
# This submission script is to submit for a planar SAM
## VARIABLES:

###SERVER_SPECIFIC_COMMANDS_START

#SBATCH -p compute
#SBATCH -t 108:00:00
#SBATCH -J _JOB_NAME_
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1              # total number of mpi tasks requested
#SBATCH --mail-user=akchew@wisc.edu
#SBATCH --mail-type=all  # email me when the job starts

## DEFINING NUMBER OF CORES
num_cores="_NUMCORES_"

###SERVER_SPECIFIC_COMMANDS_END

## VARIABLES
output_prefix="_OUTPUTPREFIX_"
total_prod_time_ps="_TOTALPRODTIME_"
run_code="_RUNCODE_"

## RUNNING CODE
bash "${run_code}" \
            "${num_cores}" \
            "${output_prefix}" \
            "${total_prod_time_ps}"


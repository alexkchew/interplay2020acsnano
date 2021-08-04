#!/bin/bash
# run_extend_job.sh
# This script extends the job

# DEFINING NUMBER OF CORES TO RUN THIS ON
num_cores="${1-28}" # number of cores you want to run with e.g. 28

# SUBMIT THIS SCRIPT using the command sbatch thisscriptname

## DEFINING PREFIX
output_prefix="${2-sam_prod}" # output prefix
# "sam_prod"

## DEFINING TOTAL PRODUCTION TIME
total_prod_time="${3-100000.000}" # Total production time
# "100000.000"

## DEFINING SIM PREFIX
sim_prefix="${output_prefix}"

## CONVERTING TPR
gmx convert-tpr -s "${sim_prefix}.tpr" \
                -until "${total_prod_time}" \
                -o "${sim_prefix}.tpr"

## SEEING IF XTC FILE EXISTS
if [ -e "${output_prefix}.xtc" ]; then
    ## RESTARTING
    gmx mdrun -nt "${num_cores}" \
              -v \
              -s "${output_prefix}.tpr" \
              -cpi "${output_prefix}.cpt" \
              -append \
              -deffnm "${output_prefix}"
else
    ## STARTING
    gmx mdrun -nt "${num_cores}" -v -deffnm "${output_prefix}"
fi



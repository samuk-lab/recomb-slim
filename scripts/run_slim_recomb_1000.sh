#!/bin/bash
# SLiM array job: distributes 24,000 parameter combinations across 100 jobs via round-robin
# Update --array, N_JOBS, paths, mem, and time as needed

#SBATCH --job-name=slim_recomb1000
#SBATCH --array=1-100
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=672:00:00
#SBATCH --output=logs/slim_%A_%a.out
#SBATCH --error=logs/slim_%A_%a.err

PARAM_FILE=slim_recomb_1000.txt
SLIM_SCRIPT=slim_recomb_stickleback.slim

N_JOBS=100  # must match --array upper bound

mkdir -p logs results
module load slim

# Parameter file columns (tab-delimited, one header row):
#   recomb_multiplier   k_adapt   gamma_shape   replicate   seed
awk -v job=${SLURM_ARRAY_TASK_ID} -v n=$N_JOBS \
    'NR>1 && (NR-2)%n==job-1' $PARAM_FILE | \
while IFS=$'\t' read -r RECOMB_MULTIPLIER K_ADAPT GAMMA_SHAPE REPLICATE SEED; do

    slim \
        -s ${SEED} \
        -d recomb_multiplier=${RECOMB_MULTIPLIER} \
        -d adapt_pop_size=${K_ADAPT} \
        -d gamma_shape=${GAMMA_SHAPE} \
        -d replicate=${REPLICATE} \
        -d seed=${SEED} \
        ${SLIM_SCRIPT}

done
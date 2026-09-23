#!/bin/bash
# SLiM array job: 400 parameter combinations (gamma_shape=0.4, k_adapt=10000,
# 4 recomb_multipliers x 100 replicates), one row per job, max 50 concurrent.

#SBATCH --job-name=slim_LD_100
#SBATCH --array=1-400%50
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=672:00:00
#SBATCH --output=logs/slim_LD_%A_%a.out
#SBATCH --error=logs/slim_LD_%A_%a.err

PARAM_FILE=slim_recomb_LD_100.txt
SLIM_SCRIPT=slim_recomb_stickleback_LD_D.slim

mkdir -p logs LD_D/results
module load slim

# Parameter file columns (tab-delimited, one header row):
#   recomb_multiplier   k_adapt   gamma_shape   replicate   seed
awk -v job=${SLURM_ARRAY_TASK_ID} 'NR==job+1' $PARAM_FILE | \
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
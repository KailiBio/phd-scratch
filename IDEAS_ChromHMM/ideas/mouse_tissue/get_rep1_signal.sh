#!/bin/bash
#SBATCH -n 3
#SBATCH -N 1
#SBATCH --mem=5G
#SBATCH --time=12:00:00
#SBATCH -o get_rep1_signal.out
#SBATCH -e get_rep1_signal.err
#SBATCH --partition=12hours
#SBATCH --job-name=get_rep1_signal
#SBATCH --array=1-68%30

echo "SLURM_JOBID:"$SLURM_JOBID
time bash /data/zusers/fankaili/ideas/code/get_rep1_signal/ENCODE_rep1_signal_code_${SLURM_ARRAY_TASK_ID}.sh

echo "Done! Congrats~"


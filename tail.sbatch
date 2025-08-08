#!/bin/bash
#SBATCH --time=1:00:00                  
#SBATCH --account=chenyliu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=UCNtau
#SBATCH --partition=physics
#SBATCH --output=./out/job_output_%A_%a.out
#SBATCH --mem=10000
#SBATCH --mail-type=All
#SBATCH --mail-user="ciyouh2@illinois.edu"
#SBATCH --array=1-11

source activate /projects/illinois/eng/physics/chenyliu/Ryan_ciyouh2/root

START=$((26308 + 500 * (${SLURM_ARRAY_TASK_ID} - 1)))
END=$((START + 500))

./Pulse_Tail ../UCNtau_2022_raw_data/ ./output/ ./betteroutput2022.json $START $END
#!/bin/bash
#SBATCH --time=10:00
#SBATCH --mem=1G
#SBATCH --output=./slurm/logs/%x-%j.log

# Set up software environment
module load StdEnv/2020   # This is already loaded, but for future compatibility...
module load gcc/9.3.0
module load python/3.9
module load scipy-stack
export R_LIBS=~/.local/R/$EBVERSIONR/
source ~/env/bin/activate

# Run Snakemake
snakemake -s snakefile_mockcomm_step3summarize.py --configfile configs/mockcomm_Round2.yaml --cores 1

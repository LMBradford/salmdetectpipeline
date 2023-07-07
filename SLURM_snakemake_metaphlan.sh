#!/bin/bash

#SBATCH --time=16:00:00
#SBATCH --mem=80G
#SBATCH --output=./slurm/logs/%x-%j.log
#SBATCH --cpus-per-task=16


# Load the required modules
module load gcc blast samtools bedtools bowtie2 python/3.10

# Generate your virtual environment in $SLURM_TMPDIR
virtualenv --no-download ${SLURM_TMPDIR}/env
source ${SLURM_TMPDIR}/env/bin/activate

# Install metaphlan and snakemake and their dependencies
pip install --no-index --upgrade pip
pip install --no-index metaphlan==4.0.3
pip install --no-index snakemake
pip install --no-index tabulate==0.8.10

#Run snakemake
snakemake -s snakefile_mockcomm_metaphlan.py --configfile configs/configfile.yaml --cores $SLURM_CPUS_PER_TASK --keep-going 

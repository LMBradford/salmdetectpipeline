#!/bin/bash

#SBATCH --time=15:00
#SBATCH --mem=6G
#SBATCH --output=./slurm/logs/%x-%j.log
#SBATCH --cpus-per-task=1

#Load bbmap module
module load gcc
module load bbmap/38.86

#Run bbmap filterbyname
filterbyname.sh in1=/home/bradford/projects/def-alexwong/bradford/MockComm_masterfiles/Mockcomm_master_R1.fq.gz in2=/home/bradford/projects/def-alexwong/bradford/MockComm_masterfiles/Mockcomm_master_R2.fq.gz names=Library/Names/${N}_names.txt include=t out1=Library/${N}_R1.fq.gz out2=Library/${N}_R2.fq.gz -Xmx6g

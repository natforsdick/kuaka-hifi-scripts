#!/bin/bash -e 
#SBATCH -A ga03186 # CHANGE!! 
#SBATCH -J fastqc # job name (shows up in the queue) 
#SBATCH -c 4
#SBATCH --mem=800M 
#SBATCH --time=00:35:00 #Walltime (HH:MM:SS) 
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err     

################### 
# FastQC 
# Nat Forsdick, 2022-06-01
################### 

# F+R fastqs of ~60G take ~2.5 hrs to run.

# This script takes two arguments, $1 : the path to the input directory, \ 
# and $2 : the prefix of the fastq file.
# e.g. to execute:
# sbatch run_array_fastqc.sl /nesi/nobackup/ga03048/data/illumina/ H07456-L1_S1_L00

#### MODULES ####
module purge
module load FastQC/0.11.9  
## MultiQC/1.9-gimkl-2020a-Python-3.8.2
#################

#### ENVIRONMENT #### 
# Example sampname: H07456-L1_S1_L002_R1_001.fastq.gz  
# We pass the L00* in as the array numbers

IN_DIR=$1
IN_PREFIX=$2

fq=.fastq
#####################

cd ${IN_DIR}

fastqc -t ${SLURM_CPUS_PER_TASK} ${IN_DIR}${IN_PREFIX}${fq} 

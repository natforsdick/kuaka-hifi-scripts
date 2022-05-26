#!/bin/bash -e

#SBATCH --account ga03048
#SBATCH --job-name trimmomatic
#SBATCH --time 01:40:00
#SBATCH --mem=10G
#SBATCH --output=%x-%j.out
#SBATCH --cpus-per-task=12
#SBATCH --error=%x-%j.err
#SBATCH --profile=task

#############
# Variables #
#############
indir=/nesi/project/ga03186/data/HiSeq-kaki-DNA1565/
outdir=/nesi/nobackup/ga03186/Kaki-HiSeq/trimmomatic/
prefix=H01392-L1_S7_L005

###############
# Modules #
###############
ml purge
module load Trimmomatic/0.39-Java-1.8.0_144

#############
# run       #
############# 

cd ${outdir}
date

trimmomatic PE -threads $SLURM_CPUS_PER_TASK ${indir}${prefix}_R1_001.fastq.gz ${indir}${prefix}_R2_001.fastq.gz \
-baseout ${outdir}${prefix}_trimmed.fastq.gz \
ILLUMINACLIP:/opt/nesi/mahuika/Trimmomatic/0.39-Java-1.8.0_144/adapters/TruSeq3-PE-2.fa:2:30:10:2:keepBothReads \
LEADING:3 TRAILING:3 SLIDINGWINDOW:5:20 MINLEN:50

date

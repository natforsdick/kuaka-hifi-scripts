#!/bin/bash -e
#SBATCH -A ga03186
#SBATCH -J racon
#SBATCH --mem=
#SBATCH --time=
#SBATCH -c
#SBATCH --outout %x.%j.out
#SBATCH --error %x.%j.err
#SBATCH --profile=task

#############
# MODULES
ml purge
ml load Racon/1.4.21-GCC-9.2.0-CUDA-11.2.0-hybrid
############

###########
# PARAMS
INRAWDIR=/nesi/project/ga03186/data/HiSeq-kaki-DNA1565/
INMAPDIR=/nesi/nobackup/ga03186/kaki-hifi-asm/asm3-hic-hifiasm-p/02-bwa/
OUTDIR=/nesi/nobackup/ga03186/kaki-hifi-asm/asm3-hic-hifiasm-p/03-racon/
INRAW=H01392-L1_S7_L005
R1=_R1_001.fastq.gz
R2=_R2_001.fastq.gz
INMAP=subsamp_hifi_bwa
ASM=

racon 

#!/bin/bash -e

#SBATCH --account ga03186
#SBATCH --job-name masurca
#SBATCH --cpus-per-task=32
#SBATCH --mem 10G 
#SBATCH --time 1-00:00:00
#SBATCH --output %x.%j.out
#SBATCH --error %x.%j.err
#SBATCH --profile=task

#############
# run-masurca.sl
# Nat Forsdick, 2021-11-01
# Masurca hybrid genome assembly
#############

############
# MODULES
ml purge
ml load MaSuRCA/4.0.5-gimkl-2020a
############

############
# PARAMS
hiseqdir=/nesi/project/ga03186/data/HiSeq-kaki-DNA1565/
hiseqpre=H01392-L1_S7_L005
hiseq=${hiseqdir}${hiseqpre}
R1=_R1_001.fastq.gz
R2=_R2_001.fastq.gz
hifiin=/nesi/project/ga03186/data/JF_PacBio-kaki-Steeves-Order260/
hifi=m54349U_210221_005741.fastq.gz
OUTDIR=/nesi/nobackup/ga03186/kaki-hifi-asm/asm5-masurca/
###########

cd $OUTDIR

echo Beginnning MaSuRCA at 
date

masurca -i ${hiseq}${R1},${hiseq}${R2} -r ${hifiin}${hifi} -t 32
date



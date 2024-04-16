#!/bin/bash -e

#SBATCH --account ga03186
#SBATCH --job-name hifiasm-kuaka
#SBATCH --cpus-per-task=36
#SBATCH --mem 120G 
#SBATCH --time 05:00:00 # When starting this job, set to 04:00:00, then change to 00:30:00 for final mem-intensive step
#SBATCH --output hifiasm.%j.out
#SBATCH --error hifiasm.%j.err
#SBATCH --profile=task

##############
# HIFIASM
# 2022-06-01
# Nat Forsdick
##############

##############
# MODULES
module purge
module load hifiasm/0.15.5-GCC-9.2.0
##############

##############
INDIR=/nesi/nobackup/ga03186/kuaka-genome/01-preprocessing/01a-adapfilt/
DATA=kuaka-hifi.filt.fastq.gz
INHIC=/nesi/project/ga03186/data/kuaka-HiC/
DATAHIC=kuaka_omnic_clean_  # suffix is _R1.fastq.gz, _R2.fastq.gz
OUTDIR=/nesi/nobackup/ga03186/kuaka-genome/asm3-hifiasm-hic/
OUTPRE=kuaka-hifiasm3
##############

mkdir -p $OUTDIR
cd $OUTDIR

echo "Starting hifiasm assembly for ${DATA} at " 
date
echo "hifiasm --primary -t 24 -o ${OUTDIR}${OUTPRE}.asm ${INDIR}${DATA}"
hifiasm --primary -t 36 -o ${OUTDIR}${OUTPRE} ${INDIR}${DATA} --h1 ${INHIC}${DATAHIC}R1.fastq.gz --h2 ${INHIC}${DATAHIC}R2.fastq.gz 2> ${OUTPRE}.log
echo "Finished at " 
date

#echo "Outputting .fa files at"
#date
#awk '/^S/{print ">"$2;print $3}' ${OUTPRE}.asm.hic.p_ctg.gfa > ${OUTPRE}.hic.p_ctg.fa
#awk '/^S/{print ">"$2;print $3}' ${OUTPRE}.asm.hic.a_ctg.gfa > ${OUTPRE}.hic.a_ctg.fa
#echo "Completed hifiasm pipeline at"
#date

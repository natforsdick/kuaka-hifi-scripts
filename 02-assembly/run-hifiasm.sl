#!/bin/bash -e

#SBATCH --account ga03186
#SBATCH --job-name hifiasm-kuaka
#SBATCH --cpus-per-task=36
#SBATCH --mem 50G 
#SBATCH --time 06:00:00 # When starting this job, set to 04:00:00, then change to 00:30:00 for final mem-intensive step
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
INDIR=/nesi/nobackup/ga03186/kuaka-genome/01-preprocessing/
DATA=kuaka-hifi.fastq
OUTDIR=/nesi/nobackup/ga03186/kuaka-genome/asm1-hifiasm/
OUTPRE=kuaka-hifiasm
##############

mkdir -p $OUTDIR
cd $OUTDIR

echo "Starting hifiasm assembly for ${DATA} at " 
date
echo "hifiasm --primary -t 24 -o ${OUTDIR}${OUTPRE}.asm ${INDIR}${DATA}"
hifiasm --primary -t 24 -o ${OUTDIR}${OUTPRE}.asm ${INDIR}${DATA}
echo "Finished at " 
date

echo "Outputting .fa files at"
date
awk '/^S/{print ">"$2;print $3}' ${OUTPRE}.asm.p_ctg.gfa > ${OUTPRE}.p_ctg.fa
awk '/^S/{print ">"$2;print $3}' ${OUTPRE}.asm.a_ctg.gfa > ${OUTPRE}.a_ctg.fa
echo "Completed hifiasm pipeline at"
date

#!/bin/bash -e

#SBATCH --account ga03186
#SBATCH --job-name hifiasm_test
#SBATCH --cpus-per-task 32
#SBATCH --mem 36G
#SBATCH --time 03:00:00
#SBATCH --output hifiasm.%j.out
#SBATCH --error hifiasm.%j.err
#SBATCH --profile=task

##############
# HIFIASM - HiFi only
# 2021-08-16
# Nat Forsdick
##############

##############
# MODULES
module purge
module load hifiasm/0.15.5-GCC-9.2.0
##############

##############
INDIR=/nesi/project/ga03186/data/JF_PacBio-kaki-Steeves-Order260/processed/
DATA=m54349U_210221_005741.fastq
OUTDIR=/nesi/nobackup/ga03186/kaki-hifi-asm/asm2-hifiasm-p/
OUTPRE=asm2-hifiasm-p
date=$(date)
##############

mkdir -p $OUTDIR
cd $OUTDIR

echo "Starting hifiasm assembly for ${DATA} at " $date

# Used -f0 for testing - used for small genomes to disable inital Bloom filters requiring 16 GB mem.
hifiasm --primary -o ${OUTPRE} -t $SLURM_CPUS_PER_TASK ${INDIR}${DATA} 2> ${OUTPRE}.log

awk '/^S/{print ">"$2;print $3}' ${OUTPRE}.p_ctg.gfa > ${OUTPRE}.p_ctg.fa  # get primary contigs in FASTA format
awk '/^S/{print ">"$2;print $3}' ${OUTPRE}.a_ctg.gfa > ${OUTPRE}.a_ctg.fa  # get alternate contigs in FASTA format

echo "Completed at " $date


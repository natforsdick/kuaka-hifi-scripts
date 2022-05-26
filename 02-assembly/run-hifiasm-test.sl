#!/bin/bash -e

#SBATCH --account ga03186
#SBATCH --job-name hifiasm_test
#SBATCH --cpus-per-task 36
#SBATCH --mem 20G
#SBATCH --time 01:00:00
#SBATCH --output hifiasm.%j.out

##############
# HIFIASM TEST
# 2021-08-16
# Nat Forsdick
##############

##############
# MODULES
module purge
module load hifiasm/0.15.5-GCC-9.2.0
##############

##############
#INDIR=/nesi/project/ga03186/data/JF_PacBio-kaki-Steeves-Order260/processed/
#DATA=m54349U_210221_005741.fastq
#INHIC=/nesi/nobackup/ga03048/data/HiC/
#DATAHIC=Kaki_Hi_C_J42W3_GAGCGCCA_L001_  # suffix is _R1.fastq.gz, _R2.fastq.gz
OUTDIR=/nesi/nobackup/ga03186/kaki-hifi-asm/asm2-hifiasm/
#OUTPRE=kaki-hifiasm
##############

cd $OUTDIR

hifiasm -o test -t $SLURM_CPUS_PER_TASK -f0 /nesi/project/ga03186/scripts/HiFi-scripts/chr11-2M.fa.gz 2> test.log
awk '/^S/{print ">"$2;print $3}' test.bp.p_ctg.gfa > test.p_ctg.fa  # get primary contigs in FASTA format


# Hi-C phasing with paired-end short reads in two FASTQ files
#hifiasm -o HG002.asm --h1 read1.fq.gz --h2 read2.fq.gz HG002-HiFi.fq.gz

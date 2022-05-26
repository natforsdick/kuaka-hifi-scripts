#!/bin/bash -e

#SBATCH --account ga03186
#SBATCH --job-name bam2fastq
#SBATCH --cpus-per-task=2
#SBATCH --mem 500M # 
#SBATCH --time 00:03:00 
#SBATCH --output %x.%j.out
#SBATCH --error %x.%j.err
#SBATCH --profile=task

##############
# run_bam2fastq.sl
# Convert PacBio bam to fastq
# 2021-10-11
# Nat Forsdick
##############

##############
# MODULES
module purge
ml SMRT-Link/8.0.0.80529-cli-tools-only
##############

##############
INDIR=/nesi/project/ga03186/data/JF_PacBio-kaki-Steeves-Order260/
DATA=m54349U_210221_005741
OUTDIR=/nesi/nobackup/ga03186/kaki-hifi/
##############

cd $OUTDIR

echo Starting at
date

# Once this is completed you will need to run:
# ml SMRT-Link/8.0.0.80529-cli-tools-only
# pbmerge -o ${DATA}.ccs.bam ${DATA}.ccs.*.bam
# pbindex ${DATA}.ccs.bam

# Then we can extract HiFi vs 'lowfi':
module purge
ml BamTools/2.5.1-gimkl-2020a
bamtools filter -in ${INDIR}${DATA}.subreads.bam -out ${DATA}.hifi_reads.bam -tag "rq":">=0.99"
bamtools filter -in ${INDIR}${DATA}.subreads.bam -out ${DATA}.lowfi_reads.bam -tag "rq":"<0.99"

# Then convert to fastq:
ml SMRT-Link/8.0.0.80529-cli-tools-only
#bam2fastq ${DATA}.hifi_reads.bam -o ${DATA}.hifi_reads 
#bam2fastq ${DATA}.lowfi_reads.bam -o ${DATA}.lowfi_reads

echo Finished at
date

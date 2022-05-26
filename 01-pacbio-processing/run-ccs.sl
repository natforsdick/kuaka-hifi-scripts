#!/bin/bash -e

#SBATCH --account ga03186
#SBATCH --job-name ccs
#SBATCH --cpus-per-task=24
#SBATCH --mem 50G # 
#SBATCH --time 03:00:00 
#SBATCH --array=1-2 # need to run 20 total, but let's just try with 1 to see the resource use.
#SBATCH --output %x.%j.out
#SBATCH --error %x.%j.err
#SBATCH --profile=task

##############
# ccs.sl
# Convert PacBio to HiFi
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

mkdir -p $OUTDIR
cd $OUTDIR

echo Starting hifi conversion date

# Hi-C phasing with paired-end short reads in two FASTQ files
ccs ${INDIR}${DATA}.subreads.bam ${DATA}.ccs.${SLURM_ARRAY_TASK_ID}.bam --chunk ${SLURM_ARRAY_TASK_ID}/20 -j $SLURM_CPUS_PER_TASK --log-file out_${SLURM_ARRAY_TASK_ID}
echo Finished at date

#############
# Once this is completed you will need to run:
# ml SMRT-Link/8.0.0.80529-cli-tools-only
# pbmerge -o ${DATA}.ccs.bam ${DATA}.ccs.*.bam
# pbindex ${DATA}.ccs.bam

# Then we can extract HiFi vs 'lowfi':
# ml BamTools/2.5.1-gimkl-2020a
# bamtools filter -in ${DATA}.ccs.bam -out ${DATA}.hifi_reads.bam -tag "rq":">=0.99"
# bamtools filter -in ${DATA}.ccs.bam -out ${DATA}.lowfi_reads.bam -tag "rq":"<0.99"

# Then convert to fastq:
# ml SMRT-Link/8.0.0.80529-cli-tools-only
# bam2fastq ${DATA}.hifi_reads.bam -o ${DATA}.hifi_reads 
# bam2fastq ${DATA}.lowfi_reads.bam -o ${DATA}.lowfi_reads

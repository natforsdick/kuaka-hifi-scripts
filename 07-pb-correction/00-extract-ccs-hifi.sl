#!/bin/bash -e

#SBATCH -A ga03186
#SBATCH -J pbccs-hifi
#SBATCH --time 18:00:00 # initial test: 10hr
#SBATCH --mem 20G 
#SBATCH --cpus-per-task 42
#SBATCH --error=%x.%j.err
#SBATCH --out=%x.%j.out
#SBATCH --profile=task

# 00-extract-ccs-hifi.sl
# Nat Forsdick, 2022-04-08
# script for generating ccs and hifi reads from PacBio bam files

##########
# MODULES
module purge
module load Anaconda3
source activate pacbio # env containing pbccs, extract-hifi, zmwfilter tools
##########

##########
# PARAMS
INDIR=/nesi/nobackup/ga03186/kaki-pacbio-data/
INBAM=m54349U_210221_005741.subreads.bam
OUTDIR=/nesi/nobackup/ga03186/kaki-pacbio-data/extract/
##########

mkdir -p ${OUTDIR}
cd $OUTDIR

filename=$(basename "$INBAM")
filename=${filename%.*.*}

# Generating CCS reads from raw PacBio subreads file

echo "Starting ccs at"
date
ccs ${INDIR}${INBAM} ${filename}.ccs.1.bam --chunk 1/6 -j 20
ccs ${INDIR}${INBAM} ${filename}.ccs.2.bam --chunk 2/6 -j 20
ccs ${INDIR}${INBAM} ${filename}.ccs.3.bam --chunk 3/6 -j 20
ccs ${INDIR}${INBAM} ${filename}.ccs.4.bam --chunk 4/6 -j 20
ccs ${INDIR}${INBAM} ${filename}.ccs.5.bam --chunk 5/6 -j 20
ccs ${INDIR}${INBAM} ${filename}.ccs.6.bam --chunk 6/6 -j 20
echo "completed ccs at"
date

#module load SAMtools

# Merging chunked CCS reads 
#samtools merge -@24 ${filename}.ccs.bam ${filename}.ccs.bam

# extract HiFi reads (>=Q20)
#extracthifi ${filename}.ccs.bam ${filename}.hifi.bam
# convert bam to fq and compress
#samtools bam2fq ${filename}.hifi.bam | gzip -c > ${filename}.hifi.fastq.gz


#!/bin/bash -e
#SBATCH -A ga03186
#SBATCH -J TGSGapClose
#SBATCH --time=01:30:00
#SBATCH --cpus-per-task=14
#SBATCH --mem=80G
#SBATCH --output %x.%j.out
#SBATCH --error %x.%j.err
#SBATCH --profile=task

# TGS-GapCloser.sl
# Nat Forsdick, 2022-08-18
# Running TGS-GapCloser for kakī 
# Takes 2 params: 1: full path to assembly, and 2: output directory.

##########
# PARAMS #
TGSGapCloser=/nesi/nobackup/ga03048/modules/TGS-GapCloser/TGS-GapCloser.sh
HIFIDIR=/nesi/project/ga03186/data/kuaka-pacbio/AGRF_CAGRF22029575_PacBio/
HIFI=kuaka-hifi
HIFIIN=$HIFIDIR${HIFI}.fastq.gz
ASM=/nesi/nobackup/ga03186/kuaka-genome/03-purge-dups/01-kuaka-hifiasm-p_ctg-purged.fa
OUTDIR=/nesi/nobackup/ga03186/kuaka-genome/04-gapclose/
##########

##########
# ENV 
ml purge
ml Racon/1.5.0-GCC-11.3.0

RACON=/opt/nesi/mahuika/Racon/1.5.0-GCC-11.3.0/bin/racon
##########

cd $OUTDIR

# First need to convert HiFi fastq to fasta
if [ ! -e ${HIFI}.fasta ]; then
zcat $HIFIIN | sed -n '1~4s/^@/>/p;2~4p' > ${HIFI}.fasta
fi

$TGSGapCloser \
        --scaff $ASM \
        --reads ${HIFI}.fasta \
        --output ${OUTDIR}kuaka-asm1-p-pur-TGC \
        --minmap_arg '-x asm20' \
        --racon $RACON \
        --tgstype pb \
        --thread 32 \
        >${OUTDIR}pipe.log 2>${OUTDIR}pipe.err

#!/bin/bash -e
#SBATCH -A ga03186
#SBATCH -J TGSGapClose
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=14
#SBATCH --mem=44G
#SBATCH --output %x.%j.out
#SBATCH --error %x.%j.err
#SBATCH --profile=task

# TGS-GapCloser.sl
# Nat Forsdick, 2022-08-18
# Running TGS-GapCloser for kuaka
# Takes 2 params: 1: full path to assembly, and 2: output directory.

##########
# PARAMS #
TGSGapCloser=/nesi/nobackup/ga03048/modules/TGS-GapCloser/TGS-GapCloser.sh
HIFIDIR=/nesi/project/ga03186/data/kuaka-pacbio/AGRF_CAGRF22029575_PacBio/
HIFI=kuaka-hifi
HIFIIN=$HIFIDIR${HIFI}.fastq.gz
ASM=/nesi/nobackup/ga03186/kuaka-genome/05-scaffolding/05_yahs/01-kuaka-hifiasm-p_ctg-purged-yahs_scaffolds_final.fa
OUTDIR=/nesi/nobackup/ga03186/kuaka-genome/gapclose/
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
        --output ${OUTDIR}01-kuaka-hifiasm-p_ctg-purged-yahs-TGS \
        --minmap_arg '-x asm20' \
        --racon $RACON \
        --tgstype pb \
        --thread 32 \
        >${OUTDIR}pipe.log 2>${OUTDIR}pipe.err

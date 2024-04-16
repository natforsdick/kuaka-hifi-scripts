#!/bin/bash -e

###################
# BUSCO - protein annotation
# Nat Forsdick
###################

# Load modules
module purge
module load BUSCO/5.6.1-gimkl-2022a

export NUMEXPR_MAX_THREADS=8
OUTDIR=/nesi/nobackup/ga03186/kuaka-genome/annotation/vPelecanoides/
INDIR=/nesi/nobackup/ga03186/kuaka-genome/annotation/vPelecanoides/
samplist=Pelecanoides-urinatrix-genomic-prot.fasta
INDB=/nesi/nobackup/ga03186/kuaka-genome/annotation/vPelecanoides/busco_downloads/lineages/aves_odb10

filename=$(basename "$samplist")
filename=${filename%.*}

cd $OUTDIR

echo "Starting BUSCO for ${samplist}"
# -f = force, -r = restart
busco -i ${INDIR}${samplist} -o BUSCO_${filename} -f --offline -l ${INDB} -m prot -c 18
echo "Finished BUSCO for ${samp}"

# To make BUSCO plots:
# ml BUSCO/5.2.2-gimkl-2020a R/4.1.0-gimkl-2020a
# generate_plot.py -wd ./



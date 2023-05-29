#!/bin/bash -e

###################
# BUSCO - INSECTS - 2020-12-13
# Nat Forsdick
###################

# Load modules
module purge
module load BUSCO/5.2.2-gimkl-2020a
#cp -r $AUGUSTUS_CONFIG_PATH ./MyAugustusConfig

# Set up environment
#export AUGUSTUS_CONFIG_PATH=/nesi/project/ga03048/scripts/QC/MyAugustusConfig
#export BUSCO_CONFIG_FILE="/nesi/project/ga03048/scripts/config.ini"

OUTDIR=/nesi/nobackup/ga03186/kuaka-genome/05-scaffolding/05b-Dovetail-OmniC/yahs/
INDIR=/nesi/nobackup/ga03186/kuaka-genome/05-scaffolding/05b-Dovetail-OmniC/yahs/
samplist='01-kuaka-hifiasm-p_ctg-purged-DT-yahs_scaffolds_final.fa' #'weta-hic-hifiasm.cns.fa'
INDB=/nesi/nobackup/ga03186/kaki-hifi-asm/asm-stats/busco_downloads/lineages/aves_odb10

for samp in $samplist
do

filename=$(basename "$samp")
filename=${filename%.*}

cd $OUTDIR

echo "Starting BUSCO for ${samp}"
# -f = force, -r = restart
	busco -i ${INDIR}${samp} -o BUSCO5_${filename} -f --offline -l ${INDB} -m geno -c 24
	echo "Finished BUSCO for ${samp}"
done

# To make BUSCO plots:
# ml BUSCO/5.2.2-gimkl-2020a R/4.1.0-gimkl-2020a
# generate_plot.py -wd ./



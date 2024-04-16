#!/bin/bash -e

#SBATCH -A ga03186
#SBATCH -J BUSCO # job name (shows up in the queue)
#SBATCH -c 8
#SBATCH --mem=5G
#SBATCH --partition=large
#SBATCH --time=00:40:00 #Walltime (HH:MM:SS)
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=forsdickn@landcareresearch.co.nz
#SBATCH --output %x.%j.out #
#SBATCH --error %x.%j.err 

###################
# BUSCO - protein annotation
# Nat Forsdick
###################

# Load modules
module purge
module load BUSCO/5.6.1-gimkl-2022a

export NUMEXPR_MAX_THREADS=10
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


#!/bin/bash -e 
#SBATCH -A ga03048
#SBATCH -J meryl-merge
#SBATCH -c 28
#SBATCH --mem=2G
#SBATCH --time=00:03:00
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err

##########
# PARAMS #
##########
genome=kaki
outdir=/nesi/nobackup/ga03186/kaki-hifi-asm/asm-stats/merqury/
##########

cd $outdir

echo "merging"
date
# 2. Merge
meryl union-sum threads=22 memory=$SLURM_MEM_PER_NODE output ${genome}.meryl ${outdir}*.meryl


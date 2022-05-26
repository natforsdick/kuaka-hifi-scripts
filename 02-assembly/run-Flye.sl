#!/bin/bash -e
#SBATCH -A ga03186
#SBATCH -J flye
#SBATCH --time=06:00:00
#SBATCH -c 46
#SBATCH --mem=50G
#SBATCH --mail-user=forsdickn@landcareresearch.co.nz
#SBATCH --mail-type=FAIL,END
#SBATCH --output %x.%j.out
#SBATCH --error %x.%j.err
#SBATCH --profile=task

# run-Flye.sl
# Nat Forsdick, 2021-08-31
# Script to trial Flye genome assembly of HiFi reads

ml purge
ml load Flye/2.9-gimkl-2020a-Python-3.8.2

INDIR=/nesi/project/ga03186/data/JF_PacBio-kaki-Steeves-Order260/processed/
OUTDIR=/nesi/nobackup/ga03186/kaki-hifi-asm/asm4-flye/

#mkdir -p $OUTDIR
cd $OUTDIR

# -g = expected genome size
# Can use --resume to resume previous run
flye --pacbio-hifi ${INDIR}m54349U_210221_005741.fastq -o $OUTDIR -g 1.2g -t $SLURM_CPUS_PER_TASK --resume

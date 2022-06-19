#!/bin/bash -e

#SBATCH --job-name=make-paf
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --time=00:15:00
#SBATCH --mem=10G
#SBATCH --ntasks=1
#SBATCH --profile=task 
#SBATCH --account=ga03186
#SBATCH --cpus-per-task=36

# Purge_dups pipeline
# Created by Sarah Bailey, UoA
# Modified by Nat Forsdick, 2021-08-24

# step 01: align HiFi sequencing data to the assembly and generate a paf file
# Takes one parameter - PRI or ALT

#########
# MODULES
module purge
module load minimap2/2.20-GCC-9.2.0 
#########

#########
# PARAMS
INDIR=/nesi/nobackup/ga03186/kuaka-genome/02-assembly/asm1-hifiasm/
OUTDIR=/nesi/nobackup/ga03186/kuaka-genome/02-assembly/asm1-hifiasm/purge_dups/
IN_DATA=/nesi/project/ga03186/data/kuaka-pacbio/AGRF_CAGRF22029575_PacBio/
DATA=kuaka-hifi.fastq.gz
PRE=kuaka-hifiasm # PREFIX
PRI=p_ctg
ALT=a_ctg
R1=01- # Designate cutoffs round - either default (01) or modified (02) and whether Primary or Alternate assembly
R2=02-
#########

mkdir -p $OUTDIR
cd $OUTDIR

if [ "$1" == "PRI" ]; then
  minimap2 -x map-hifi -t 24 ${INDIR}${PRE}.${PRI}.fa ${IN_DATA}${DATA} | gzip -c - > ${R1}${PRE}-${PRI}-mapped.paf.gz
 
elif [ "$1" == "ALT" ]; then
  minimap2 -x map-hifi -t 24 ${INDIR}${R1}${PRE}.${ALT}.hap-merged.fa ${IN_DATA}${DATA} | gzip -c - > ${R1}${PRE}-${ALT}-merged-mapped.paf.gz

else
  minimap2 -x map-hifi -t 24 ${INDIR}${PRE}.fasta ${IN_DATA}${DATA} | gzip -c - > ${R1}${PRE}-mapped.paf.gz

fi

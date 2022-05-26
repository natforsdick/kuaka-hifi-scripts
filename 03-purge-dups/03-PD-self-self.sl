#!/bin/bash -e

#SBATCH --job-name=self-self
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --time=00:03:00
#SBATCH --mem=10G
#SBATCH --ntasks=1
#SBATCH --profile=task 
#SBATCH --account=ga03186
#SBATCH --cpus-per-task=16

# Purge_dups pipeline
# Created by Sarah Bailey, UoA
# Modified by Nat Forsdick, 2021-08-24

# step 03: do a self-self alignment

#########
# MODULES
module purge
module load minimap2/2.20-GCC-9.2.0 
#########

#########
# PARAMS
OUTDIR=/nesi/nobackup/ga03186/kaki-hifi-asm/asm5-masurca/purge_dups/
PRE=primary.genome.scf # PREFIX
PRI=p_ctg
ALT=a_ctg
R1=01- # Designate cutoffs round - either default (01) or modified (02) and whether Primary or Alternate assembly
R2=02-
#########

cd $OUTDIR
# -x asm5: intra-specific asm-to-asm alignment
if [ "$1" == "PRI" ]; then 
  minimap2 -x asm5 -t $SLURM_CPUS_PER_TASK -DP ${R1}${PRE}${PRI}.split ${R1}${PRE}${PRI}.split | gzip -c - > ${R1}${PRE}${PRI}.split.self.paf.gz
elif [ "$1" == "ALT" ]; then
  minimap2 -x asm5 -t $SLURM_CPUS_PER_TASK -DP ${R1}${PRE}${ALT}.split ${R1}${PRE}${ALT}.split | gzip -c - > ${R1}${PRE}${ALT}.split.self.paf.gz
else
  minimap2 -x asm5 -t $SLURM_CPUS_PER_TASK -DP ${R1}${PRE}.split ${R1}${PRE}.split | gzip -c - > ${R1}${PRE}.split.self.paf.gz
fi 

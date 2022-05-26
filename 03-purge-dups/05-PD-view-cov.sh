#!/bin/bash -e

# Purge_dups pipeline
# Created by Sarah Bailey, UoA
# Modified by Nat Forsdick, 2021-08-24

# step 05: to view the coverage distribution

##########
# PARAMS
PURGE_DUPS=/nesi/nobackup/ga03186/purge_dups/scripts/
OUTDIR=/nesi/nobackup/ga03186/kaki-hifi-asm/asm5-masurca/purge_dups/
PRE=primary.genome.scf # PREFIX
PRI=p_ctg
ALT=a_ctg
R1=01- # Designate cutoffs round - either default (01) or modified (02) and whether Primary or Alternate assembly
R2=02-
ASMSTATS=/nesi/project/ga03186/scripts/Assemblathon_scripts/assemblathon_stats.pl
##########

##########
# MODULES
ml purge
ml load Python
##########

cd ${OUTDIR}

if [ "$1" == "PRI" ]; then
 
  if [ "$2" == "R1" ]; then
  
    mv purged.fa ${R1}${PRE}${PRI}-purged.fa
    mv hap.fa ${R1}${PRE}${PRI}-hap.fa
    python3 ${PURGE_DUPS}hist_plot.py -c ${R1}${PRE}${PRI}-cutoffs ${R1}${PRE}${PRI}-PB.stat ${R1}${PRE}${PRI}-PB.cov.png

    # Run assemblathon stats
    $ASMSTATS ${R1}${PRE}${PRI}-purged.fa > ${R1}${PRE}${PRI}-purged.stats

  elif [ "$2" == "R2" ]; then
    mv purged.fa ${R2}${PRE}${PRI}-purged.fa
    mv hap.fa ${R2}${PRE}${PRI}-hap.fa
    python3 ${PURGE_DUPS}hist_plot.py -c ${R2}${PRE}${PRI}-cutoffs ${R1}${PRE}${PRI}-PB.stat ${R2}${PRE}${PRI}-PB.cov.png

    $ASMSTATS ${R2}${PRE}${PRI}-purged.fa > ${R2}${PRE}${PRI}-purged.stats
  fi

elif [ "$1" == "ALT" ]; then
  mv purged.fa ${R1}${PRE}${ALT}-purged.fa
  mv hap.fa ${R1}${PRE}${ALT}-hap.fa
  if [ "$2" == "R1" ]; then
    python3 ${PURGE_DUPS}hist_plot.py -c ${R1}${PRE}${ALT}-cutoffs ${R1}${PRE}${ALT}-PB.stat ${R1}${PRE}${ALT}-PB.cov.png

    # Run assemblathon stats
    $ASMSTATS ${R1}${PRE}${ALT}-purged.fa > ${R1}${PRE}${ALT}-purged.stats

  elif [ "$2" == "R2" ]; then
    mv purged.fa ${R2}${PRE}${ALT}-purged.fa
    mv hap.fa ${R2}${PRE}${ALT}-hap.fa
    python3 ${PURGE_DUPS}hist_plot.py -c ${R2}${PRE}${ALT}-cutoffs ${R1}${PRE}${ALT}-PB.stat ${R2}${PRE}${ALT}-PB.cov.png

    $ASMSTATS ${R2}${PRE}${ALT}-purged.fa > ${R2}${PRE}${ALT}-purged.stats
  fi

else
  if [ "$1" == "R1" ]; then

    mv purged.fa ${R1}${PRE}-purged.fa
    mv hap.fa ${R1}${PRE}-hap.fa
    python3 ${PURGE_DUPS}hist_plot.py -c ${R1}${PRE}-cutoffs ${R1}${PRE}-PB.stat ${R1}${PRE}-PB.cov.png

    # Run assemblathon stats
    $ASMSTATS ${R1}${PRE}-purged.fa > ${R1}${PRE}-purged.stats

  elif [ "$1" == "R2" ]; then
    mv purged.fa ${R2}${PRE}-purged.fa
    mv hap.fa ${R2}${PRE}-hap.fa
    python3 ${PURGE_DUPS}hist_plot.py -c ${R2}${PRE}-cutoffs ${R1}${PRE}-PB.stat ${R2}${PRE}-PB.cov.png

    $ASMSTATS ${R2}${PRE}-purged.fa > ${R2}${PRE}-purged.stats
  fi
fi

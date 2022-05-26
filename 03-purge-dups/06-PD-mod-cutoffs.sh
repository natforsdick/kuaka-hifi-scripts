#!/bin/bash -e

# Purge_dups pipeline
# Created by Sarah Bailey, UoA
# Modified by Nat Forsdick, 2021-08-24

# step 06: modify cutoffs

##########
# PARAMS
PURGE_DUPS=/nesi/nobackup/ga03186/purge_dups/bin/
OUTDIR=/nesi/nobackup/ga03186/kaki-hifi-asm/asm5-masurca/purge_dups/
PRE=primary.genome.scf # PREFIX
PRI=p_ctg
ALT=a_ctg
R1=01-
R2=02- # Designate cutoffs round - either default (01) or modified (02) and whether Primary or Alternate assembly
CUTOFFS="-l2 -m15 -u60"
##########

cd ${OUTDIR}
echo $CUTOFFS
${PURGE_DUPS}calcuts ${CUTOFFS} ${R1}${PRE}-PB.stat > ${R2}${PRE}-cutoffs

# Following this, you need to run steps 04-07 with $ROUND modified for new cutoffs.

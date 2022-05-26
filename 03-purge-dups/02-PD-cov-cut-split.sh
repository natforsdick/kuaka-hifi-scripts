#!/bin/bash -e

# Purge_dups pipeline
# Created by Sarah Bailey, UoA
# Modified by Nat Forsdick, 2021-08-24

# step 02: calculate read depth histogram and base-level read depth, generate default cutoffs, split the assembly
# Takes one argument: PRI or ALT

#########
# PARAMS
INDIR=/nesi/nobackup/ga03186/kaki-hifi-asm/asm5-masurca/CA.mr.99.17.15.0.02/
OUTDIR=/nesi/nobackup/ga03186/kaki-hifi-asm/asm5-masurca/purge_dups/
PURGE_DUPS=/nesi/nobackup/ga03186/purge_dups/bin/
PRE=primary.genome.scf # PREFIX
PRI=p_ctg
ALT=a_ctg
R1=01- # Designate cutoffs round - either default (01) or modified (02) and whether Primary or Alternate assembly
R2=02-
#########

cd $OUTDIR

if [ "$1" == "PRI" ]; then 
# step 02a: Produce PB.base.cov and PB.stat files
${PURGE_DUPS}pbcstat ${R1}${PRE}-${PRI}-mapped.paf.gz

mv PB.stat ${R1}${PRE}-${PRI}-PB.stat
mv PB.base.cov ${R1}${PRE}-${PRI}-PB.base.cov
mv PB.cov.wig ${R1}${PRE}-${PRI}-PB.cov.wig

## step 02b: generate default cutoffs
${PURGE_DUPS}calcuts ${R1}${PRE}-${PRI}-PB.stat > ${R1}${PRE}-${PRI}-cutoffs 2> ${R1}${PRE}-${PRI}-calcults.log

## step 02c: split the assembly 
${PURGE_DUPS}split_fa ${INDIR}${PRE}.${PRI}.fa > ${R1}${PRE}-${PRI}.split

elif [ "$1" == "ALT" ]; then
  ${PURGE_DUPS}pbcstat ${R1}${PRE}-${ALT}-merged-mapped.paf.gz

  mv PB.stat ${R1}${PRE}-${ALT}-PB.stat
  mv PB.base.cov ${R1}${PRE}-${ALT}-PB.base.cov
  mv PB.cov.wig ${R1}${PRE}-${ALT}-PB.cov.wig
  
  ${PURGE_DUPS}calcuts ${R1}${PRE}-${ALT}-PB.stat > ${R1}${PRE}-${ALT}-cutoffs 2> ${R1}${PRE}-${ALT}-calcults.log

  ${PURGE_DUPS}split_fa ${INDIR}${R1}${PRE}.${ALT}.hap-merged.fa > ${R1}${PRE}-${ALT}.split

else
# step 02a: Produce PB.base.cov and PB.stat files
${PURGE_DUPS}pbcstat ${R1}${PRE}-mapped.paf.gz

mv PB.stat ${R1}${PRE}-PB.stat
mv PB.base.cov ${R1}${PRE}-PB.base.cov
mv PB.cov.wig ${R1}${PRE}-PB.cov.wig

## step 02b: generate default cutoffs
${PURGE_DUPS}calcuts ${R1}${PRE}-PB.stat > ${R1}${PRE}-cutoffs 2> ${R1}${PRE}-calcults.log

## step 02c: split the assembly
${PURGE_DUPS}split_fa ${INDIR}${PRE}.fasta > ${R1}${PRE}.split
fi

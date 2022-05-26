#!/bin/bash -e

# Purge_dups pipeline
# Created by Sarah Bailey, UoA
# Modified by Nat Forsdick, 2021-08-24

# Purge_dups pipeline
# Purge haplotigs and overlaps
# Get purged primary and haplotig sequences from draft assembly
# Takes two arguments: PRI/ALT and R1/R2

##########
# PARAMS
##########
INDIR=/nesi/nobackup/ga03186/kaki-hifi-asm/asm5-masurca/CA.mr.99.17.15.0.02/
OUTDIR=/nesi/nobackup/ga03186/kaki-hifi-asm/asm5-masurca/purge_dups/
PURGE_DUPS=/nesi/nobackup/ga03186/purge_dups/bin/
PRE=primary.genome.scf # PREFIX
PRI=p_ctg
ALT=a_ctg
R1=01-
R2=02- # Designate cutoffs round - either default (01) or modified (02) and whether Primary or Alternate assembly
##########

cd $OUTDIR

# R1 version
if [ "$1" == "PRI" ]; then
  if [ "$2" == "R1" ]; then
    # Step 04a: Purge haplotigs and overlaps
    ${PURGE_DUPS}purge_dups -2 -T ${R1}${PRE}-${PRI}-cutoffs -c ${R1}${PRE}-${PRI}-PB.base.cov ${R1}${PRE}-${PRI}.split.self.paf.gz > ${R1}${PRE}-${PRI}-dups.bed 2> ${R1}${PRE}-${PRI}-purge_dups.log

    # Step 04b: Get purged primary and haplotig sequences from draft assembly
    ${PURGE_DUPS}get_seqs -e ${R1}${PRE}-${PRI}-dups.bed ${INDIR}${PRE}.${PRI}.fa
  
  elif [ "$2" == "R2" ]; then
    ${PURGE_DUPS}purge_dups -2 -T ${R2}${PRE}-${PRI}-cutoffs -c ${R1}${PRE}-${PRI}-PB.base.cov ${R1}${PRE}-${PRI}.split.self.paf.gz > ${R2}${PRE}-${PRI}-dups.bed 2> ${R2}${PRE}-${PRI}-purge_dups.log

    ${PURGE_DUPS}get_seqs -e ${R2}${PRE}-${PRI}-dups.bed ${INDIR}${PRE}.${PRI}.fa
  fi

elif [ "$1" == "ALT" ]; then
  if [ "$2" == "R1" ]; then
    ${PURGE_DUPS}purge_dups -2 -T ${R1}${PRE}-${ALT}-cutoffs -c ${R1}${PRE}-${ALT}-PB.base.cov ${R1}${PRE}-${ALT}.split.self.paf.gz > ${R1}${PRE}-${ALT}-dups.bed 2> ${R1}${PRE}-${ALT}-purge_dups.log

    ${PURGE_DUPS}get_seqs -e ${R1}${ALT}-${PRI}-dups.bed ${INDIR}${R1}${PRE}.${ALT}.hap-merged.fa

  elif [ "$2" == "R2" ]; then
    ${PURGE_DUPS}purge_dups -2 -T ${R2}${PRE}${ALT}-cutoffs -c ${R1}${PRE}${ALT}-PB.base.cov ${R1}${PRE}${ALT}.split.self.paf.gz > ${R2}${PRE}${ALT}-dups.bed 2> ${R2}${PRE}${ALT}-purge_dups.log

    ${PURGE_DUPS}get_seqs -e ${R2}${PRE}${ALT}-dups.bed ${INDIR}${R2}${PRE}.${ALT}.hap-merged.fa
  fi

else
if [ "$1" == "R1" ]; then
    # Step 04a: Purge haplotigs and overlaps
    ${PURGE_DUPS}purge_dups -2 -T ${R1}${PRE}-cutoffs -c ${R1}${PRE}-PB.base.cov ${R1}${PRE}.split.self.paf.gz > ${R1}${PRE}-dups.bed 2> ${R1}${PRE}-purge_dups.log

    # Step 04b: Get purged primary and haplotig sequences from draft assembly
    ${PURGE_DUPS}get_seqs -e ${R1}${PRE}-dups.bed ${INDIR}${PRE}.fasta

  elif [ "$1" == "R2" ]; then
    ${PURGE_DUPS}purge_dups -2 -T ${R2}${PRE}-cutoffs -c ${R1}${PRE}-PB.base.cov ${R1}${PRE}.split.self.paf.gz > ${R2}${PRE}-dups.bed 2> ${R2}${PRE}-purge_dups.log

    ${PURGE_DUPS}get_seqs -e ${R2}${PRE}-dups.bed ${INDIR}${PRE}.fasta
  fi

fi

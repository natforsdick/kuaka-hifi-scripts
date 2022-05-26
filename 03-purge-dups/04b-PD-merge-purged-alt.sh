#!/bin/bash -e

# Purge_dups pipeline
# Created by Sarah Bailey, UoA
# Modified by Nat Forsdick, 2021-08-24

# step 04b: merge purged haplotigs with alternate draft assembly

##########
# PARAMS
##########
INDIR=/nesi/nobackup/ga03186/kaki-hifi-asm/asm3-hic-hifiasm-p/
OUTDIR=/nesi/nobackup/ga03186/kaki-hifi-asm/asm3-hic-hifiasm-p/purge_dups/
PURGE_DUPS=/nesi/nobackup/ga03186/purge_dups/bin/
PRE=asm3-hic-hifiasm-p # PREFIX
PRI=p_ctg
ALT=a_ctg
R1=01P- # Designate cutoffs round - either default (01) or modified (02) and whether Primary or Alternate assembly
R2=02P-
##########

cd $OUTDIR

# Rename outputs
if [ "$1" == "PRI" ]; then

mv purged.fa ${R1}${PRE}${PRI}-purged.fa
mv hap.fa ${R1}${PRE}${PRI}-hap.fa

# Step 05: Merge purged haplotigs (hap.fa) with the alternate draft assembly
cat ${R1}${PRE}-${PRI}-hap.fa ${INDIR}${PRE}.${ALT}.fa > ${INDIR}${R1}${PRE}.${PRI}.hap-merged.fa

elif [ "$1" == "ALT" ]; then

mv purged.fa ${R1}${PRE}${ALT}-purged.fa
mv hap.fa ${R1}${PRE}${ALT}-hap.fa

# Step 05: Merge purged haplotigs (hap.fa) with the alternate draft assembly
cat ${R1}${PRE}-${ALT}-hap.fa ${INDIR}${PRE}.${ALT}.fa > ${INDIR}${R1}${PRE}.${ALT}.hap-merged.fa
fi

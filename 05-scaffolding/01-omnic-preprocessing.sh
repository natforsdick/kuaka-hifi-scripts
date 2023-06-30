#!/bin/bash -e

ml purge 
ml SAMtools/1.15.1-GCC-11.3.0 BWA/0.7.17-GCC-11.3.0

REF=01-kuaka-hifiasm-p_ctg-purged # prefix
OUTDIR=/nesi/nobackup/ga03186/kuaka-genome/05-scaffolding/05b-Dovetail-OmniC/
TMPDIR=/nesi/nobackup/ga03186/tmp-omnic

mkdir $TMPDIR

cd $OUTDIR

samtools faidx $REF.fa

cut -f1,2 $REF.fa.fai > $REF.genome

bwa index $REF.fa

#!/bin/bash -e

INDIR=/nesi/nobackup/ga03186/kuaka-genome/ref-genomes/Pelecanoides-urinatrix/
OUTDIR=/nesi/nobackup/ga03186/kuaka-genome/annotation/vPelecanoides/
FASTA=Pelecanoides-urinatrix-genomic
GFF=Pelecanoides-urinatrix-genomic.gff

cd $OUTDIR
ml purge
ml AGAT/1.0.0-gimkl-2022a-Perl-5.34.1-R-4.2.1

agat_sp_extract_sequences.pl -g ${INDIR}${GFF} -f  ${INDIR}${FASTA}.fna -p -o ${FASTA}-prot.fasta
ml purge

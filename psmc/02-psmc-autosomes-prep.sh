#!/bin/bash -e

# we need to work out which scaffolds are sex-linked, and remove these
# using D-GENIES to identify scaffolds that are syntenic to the Z chromosome in the chicken genome assembly
# Chicken genome assembly GRCg6a has Z chromosome represented by CM000122.
# We need to extract a list of Z-chrom scaffolds, and a list of the autosomal scaffolds, which will be what we base our inferences on
# Extracted the D-GENIES output as an association table (tsv file). Now going to remove scaffolds if they have substantial match (>50%) to Z chrom

# get match length
OUTDIR=/nesi/nobackup/ga03186/kuaka-genome/psmc/
cd $OUTDIR

awk 'BEGIN {OFS="\t"}
{$1=$1}
NR == 1 {
  print $0, "match_length"
  next
} 
{
  print $0, $6-$5
}' 01-kuaka-hifiasm-p_ctg-purged-clean-omnic-mapped.PT-yahsNMC_scaffolds_final_GCA_000002315.5_GRCg6a_genomic_assoc.tsv > 01-kuaka-hifiasm-p_ctg-purged-clean-omnic-mapped.PT-yahsNMC_scaffolds_final_GCA_000002315.5_GRCg6a_genomic_assoc_length.tsv

# divide by length of query scaffold to get match proportion
awk 'BEGIN {OFS="\t"}
{$1=$1}
NR == 1 {
  print $0, "proportion"
  next
} 
{
  print $0, $10/$4*100
}' 01-kuaka-hifiasm-p_ctg-purged-clean-omnic-mapped.PT-yahsNMC_scaffolds_final_GCA_000002315.5_GRCg6a_genomic_assoc_length.tsv > 01-kuaka-hifiasm-p_ctg-purged-clean-omnic-mapped.PT-yahsNMC_scaffolds_final_GCA_000002315.5_GRCg6a_genomic_assoc_proportion.tsv

# extract those scaffolds with >50% match to the Z-chrom
grep "CM000122.5" 01-kuaka-hifiasm-p_ctg-purged-clean-omnic-mapped.PT-yahsNMC_scaffolds_final_GCA_000002315.5_GRCg6a_genomic_assoc_proportion.tsv | mawk '$11 > 50' | cut -f1 > chromosome_scaffolds_Z.txt
wc -l chromosome_scaffolds_Z.txt

# add > to start of line for scaff names
sed 's/^/>/g' chromosome_scaffolds_Z.txt > chromosome_scaffolds_Z2.txt
mv chromosome_scaffolds_Z2.txt chromosome_scaffolds_Z.txt
# extract remaining (non-Z) scaff names
REF=01-kuaka-hifiasm-p_ctg-purged-clean-omnic-mapped.PT-yahsNMC_scaffolds_final.fa
grep "^>" ${REF} | grep -v -f chromosome_scaffolds_Z.txt > chromosome_scaffolds_auto.txt

# Lets get the length of each scaffold of the reference file
module purge
module load bioawk/1.0

bioawk -c fastx '{print ">" $name ORS length($seq)}' ${REF} | paste - - > kuaka-scaffold-lengths_ensembl.txt

# For downstream analysis we need to make bed files
cut -f1 chromosome_scaffolds_Z.txt | grep -f - kuaka-scaffold-lengths_ensembl.txt | sed 's,>,,' | sed 's,\.1,\.1\t0,' > chromosome_scaffolds_Z.bed
cut -f1 chromosome_scaffolds_auto.txt | grep -f - kuaka-scaffold-lengths_ensembl.txt| sed 's,>,,' | sed 's,\.1,\.1\t0,' > chromosome_scaffolds_auto.bed

# and add in second column zeroes to make proper bed format
awk 'OFS="\t"{print $1, "0", $2}' chromosome_scaffolds_auto.bed > chromosome_scaffolds_auto2.bed
awk 'OFS="\t"{print $1, "0", $2}' chromosome_scaffolds_Z.bed > chromosome_scaffolds_Z2.bed

mv chromosome_scaffolds_auto2.bed chromosome_scaffolds_auto.bed
mv chromosome_scaffolds_Z2.bed chromosome_scaffolds_Z.bed

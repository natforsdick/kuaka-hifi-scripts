#!/bin/bash -e

OUTDIR=/nesi/nobackup/ga03186/kuaka-genome/05-scaffolding/05b-Dovetail-OmniC/mapped-masked/
REFDIR=/nesi/nobackup/ga03186/kuaka-genome/05-scaffolding/05b-Dovetail-OmniC/
REF=01-kuaka-hifiasm-p_ctg-purged # excl .fa/.fasta
BAMDIR=/nesi/nobackup/ga03186/kuaka-genome/05-scaffolding/05b-Dovetail-OmniC/
BAM=01-kuaka-hifiasm-p_ctg-purged-clean-omnic-mapped.PT # excl .bam
SDUST=/nesi/project/ga03186/bin/minimap-sdust/minimap/sdust

mkdir $OUTDIR
cd $OUTDIR

$date

# masking genome
$SDUST ${REFDIR}${REF}.fa > ${OUTDIR}${REF}-masked.bed
$date
echo regions masked from $REF

# then remove mappings that intersect with masked regions
ml purge && ml BEDTools/2.30.0-GCC-11.3.0 SAMtools/1.16.1-GCC-11.3.0

# -v = only report those with no overlap -f = minimum overlap required of a, -F = minimum overlap required of b, -e = minimum fraction of a OR b
bedtools intersect -a ${BAMDIR}${BAM}.bam -b ${OUTDIR}${REF}-masked.bed -ubam -v -f 0.8 -F 0.8 -e > ${OUTDIR}${BAM}.masked.bam
$date
echo intersections removed

# collect stats
regions=$(wc -l ${OUTDIR}${REF}-masked.bed)
raw=$(samtools view -c ${BAMDIR}${BAM}.bam)
masked=$(samtools view -c ${OUTDIR}${BAM}.masked.bam)
perc_retained=`echo "scale=4;($masked/$raw)*100" | bc`

echo "total masked regions: $regions"
echo "mapped: $raw"
echo "remaining: $masked"        
echo "% retained: $perc_retained"
$date


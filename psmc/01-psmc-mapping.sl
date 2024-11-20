#!/bin/bash -e
#SBATCH -A ga03186
#SBATCH -J psmc-mapping
#SBATCH --cpus-per-task=24
#SBATCH --mem=24G
#SBATCH -t 24:00:00
#SBATCH --out %x.%j.out
#SBATCH --err %x.%j.err

REFDIR=/nesi/nobackup/ga03186/kuaka-genome/05-scaffolding/05b-Dovetail-OmniC/all-data-yahs/
REF=01-kuaka-hifiasm-p_ctg-purged-clean-omnic-mapped.PT-yahsNMC_scaffolds_final.fa
DATADIR=/nesi/nobackup/ga03186/kuaka-pop-gen/output/03-merged/
OUTDIR=/nesi/nobackup/ga03186/kuaka-genome/psmc/

module load BWA/0.7.18-GCC-12.3.0 SAMtools/1.19-GCC-12.3.0 picard/2.26.10-Java-11.0.4

# index reference
cd $OUTDIR
if [ ! -e $REF ]; then
ln -s $REFDIR$REF $REF
fi

if [ ! -e $REF.sa ]; then
bwa index $REF
fi

# want to run eveything after this as a loop, once for KWH, once for CDP
if [ ! -e CDP_D206823_map_kuaka_ref.bam ]; then
echo mapping KWH
# map to ref - KWH
bwa mem -M -t 12 -R "@RG\tID:A01193\tSM:D206824\tLB:IlluminaWGS\tPL:ILLUMINA" $REF \
${DATADIR}D206824_S7_val_1.fq ${DATADIR}D206824_S7_val_2.fq |\
samtools view -@ 12 -bh - |\
samtools sort -@ 12 -T tmp -o KWH_D206824_map_kuaka_ref.bam

# map to ref - CDP
echo mapping CDP
bwa mem -M -t 12 -R "@RG\tID:A01193\tSM:D206823\tLB:IlluminaWGS\tPL:ILLUMINA" $REF \
${DATADIR}D206823_S67_val_1.fq ${DATADIR}D206823_S67_val_2.fq |\
samtools view -@ 12 -bh - |\
samtools sort -@ 12 -T tmp -o CDP_D206823_map_kuaka_ref.bam
fi

for bam in *_map_kuaka_ref.bam
do
    filename=$(basename "$bam") 
    filename=${filename%.*}
    # filter mapped data
    echo sorting $filename
    samtools view -@ 16 -bh -F 4 -q 30 ${bam} |
    samtools sort -@ 16 -T ${filename}_filtered_temp -o ${filename}_filtered_sorted.bam

    # remove PCR duplicates
    echo finding dups $filename
    picard MarkDuplicates \
        --MAX_FILE_HANDLES_FOR_READ_ENDS_MAP 900 \
        -I ${filename}_filtered_sorted.bam \
        -O ${filename}_filtered_sorted_rmdup.bam \
        --ASSUME_SORTED TRUE \
        --REMOVE_DUPLICATES true \
        --METRICS_FILE ${filename}.rmdup.metrics.txt \
        --TMP_DIR ./ \
        --VALIDATION_STRINGENCY SILENT

    # index final bam
    echo indexing $filename
    samtools index -@ 24 ${filename}_filtered_sorted_rmdup.bam

    # collect alignment stats
    echo collecting stats for $filename
    samtools flagstat -@ 24 ${filename}_filtered_sorted_rmdup.bam > ${filename}_filtered_sorted_rmdup.flagstat

    picard CollectAlignmentSummaryMetrics -R $REF -I ${filename}_filtered_sorted_rmdup.bam -O ${filename}_filtered_sorted_rmdup.stats.txt
    echo completed processing $filename
done

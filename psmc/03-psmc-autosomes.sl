#!/bin/bash -e
#SBATCH -A ga03186
#SBATCH -J psmc-auto
#SBATCH --cpus-per-task=6
#SBATCH --mem=6G
#SBATCH -t 15:00:00
#SBATCH --out %x.%j.out
#SBATCH --err %x.%j.err

# prepping consensus fq for autosomal regions as input for PSMC

REF=01-kuaka-hifiasm-p_ctg-purged-clean-omnic-mapped.PT-yahsNMC_scaffolds_final.fa
OUTDIR=/nesi/nobackup/ga03186/kuaka-genome/psmc/
cd $OUTDIR

for bam in *_filtered_sorted_rmdup.bam
do
    filename=$(basename "$bam")
    filename=${filename%.*}

    echo "extracting autosomal regions for $filename"
    ml purge
    module load SAMtools/1.19-GCC-12.3.0
# extract mapped autosomal regions
    samtools view -@ 24 -b -L chromosome_scaffolds_auto.bed ${filename}.bam > ${filename}_auto.bam
done &&

wait

# create consensus fq file
ml purge
module load BCFtools/1.19-GCC-11.3.0
echo "creating consensus for KWH"
bcftools mpileup --threads 12 -C50 -f ${REF} KWH_D206824_map_kuaka_ref_filtered_sorted_rmdup_auto.bam | bcftools call --threads 12 -c - | vcfutils.pl vcf2fq -d 13 -D 90 | gzip > KWH_D206824_map_kuaka_ref_diploid.fq.gz
echo "creating consensus for CDP"
bcftools mpileup --threads 12 -C50 -f ${REF} CDP_D206823_map_kuaka_ref_filtered_sorted_rmdup_auto.bam | bcftools call --threads 12 -c - | vcfutils.pl vcf2fq -d 16 -D 105 | gzip > CDP_D206823_map_kuaka_ref_diploid.fq.gz
ml purge

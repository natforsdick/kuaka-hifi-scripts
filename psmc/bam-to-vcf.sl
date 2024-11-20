#!/bin/bash -e
#SBATCH -A ga03186
#SBATCH -J psmc-mapping
#SBATCH --cpus-per-task=8
#SBATCH --mem=14G
#SBATCH -t 8:00:00
#SBATCH --out %x.%j.out
#SBATCH --err %x.%j.err

REFDIR=/nesi/nobackup/ga03186/kuaka-genome/05-scaffolding/05b-Dovetail-OmniC/all-data-yahs/
REF=01-kuaka-hifiasm-p_ctg-purged-clean-omnic-mapped.PT-yahsNMC_scaffolds_final.fa
DATADIR=/nesi/nobackup/ga03186/kuaka-pop-gen/output/03-merged/
OUTDIR=/nesi/nobackup/ga03186/kuaka-genome/psmc/
samplist="CDP_D206823_map_kuaka_ref_filtered_sorted_rmdup_auto KWH_D206824_map_kuaka_ref_filtered_sorted_rmdup_auto"

ml BCFtools/1.19-GCC-11.3.0
cd $OUTDIR
for samp in $samplist
do
    echo processing $samp
    bcftools mpileup -Ob --threads 8 -o ${samp}.bcf -f $REFDIR$REF ${samp}.bam
    bcftools call -vmO v --threads 8 -o ${samp}.vcf ${samp}.bcf
done

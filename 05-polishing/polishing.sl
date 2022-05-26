#!/bin/bash -e
#SBATCH -A ga03186
#SBATCH -J kaki-polish
#SBATCH --cpus-per-task=36
#SBATCH --mem=32G
#SBATCH --partition=large
#SBATCH --time=02:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=forsdickn@landcareresearch.co.nz
#SBATCH --output %x.%j.out
#SBATCH --error %x.%j.err
#SBATCH --profile=task

#######
# PARAMS
asmdir=/nesi/nobackup/ga03186/kaki-hifi-asm/asm3-hic-hifiasm-p/03-polishing/
fo=01P-asm3-hic-hifiasm-p-p_ctg-purged.cns
datadir=/nesi/project/ga03186/data/JF_PacBio-kaki-Steeves-Order260/processed/
outdir=/nesi/nobackup/ga03186/kaki-hifi-asm/asm-hic-hifiasm-p/03-polishing/
######

ml purge
ml wtdbg/2.5-GCC-9.2.0 BWA/0.7.17-GCC-9.2.0 \
SAMtools/1.13-GCC-9.2.0 minimap2/2.20-GCC-9.2.0

if [ ! -e ${outdir} ]; then
	mkdir -p $outdir
fi

cd $outdir

date
echo "minimap2 -I 24G -t $SLURM_CPUS_PER_TASK -ax map-hifi -r2k ${asmdir}${fo}.fa ${datadir}*.fastq.gz |\
samtools sort -@4 -o ${fo}.bam"
minimap2 -I 24G -t $SLURM_CPUS_PER_TASK -ax map-hifi -r2k ${asmdir}${fo}.fa ${datadir}*.fastq.gz |\
samtools sort -@4 -o ${fo}.bam

date
echo "samtools view -F0x900 ${fo}.bam | wtpoa-cns -t $SLURM_CPUS_PER_TASK \
-d ${asmdir}${fo}.fa -i - -fo ${fo}.cns.fa"
samtools view -F0x900 ${fo}.bam | wtpoa-cns -t $SLURM_CPUS_PER_TASK \
-d ${asmdir}${fo}.fa -i - -fo ${fo}.pol2.fa

date
echo "bwa index ${fo}.pol2.fa"
bwa index ${fo}.pol2.fa

# If you wanted to do this with short-read data (probably not recommended with HiFi, as you
# may incur mismapping across repetitive regions from Illumina)≈
#bwa mem -t 24 ${fo}.cns.fa ${srdir}/*fastq.gz | \
#samtools sort -O SAM | wtpoa-cns -t 24-x sam-sr \
#-d ${fo}.cns.fa -i - -fo ${fo}.srp.fa 
# wtpoa-cns calls in short-read data for polishing - would want to polish with CCS first



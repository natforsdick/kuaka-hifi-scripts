#!/bin/bash -e
#SBATCH --account=ga03186
#SBATCH --job-name=fastp 
#SBATCH --cpus-per-task=12 
#SBATCH --mem=16G
#SBATCH --time=2:00:00 
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=forsdickn@landcareresearch.co.nz
#SBATCH --output %x.%j.out # CHANGE number for new run
#SBATCH --error %x.%j.err #  CHANGE number for new run

##########
# PARAMS
OUTDIR=/nesi/nobackup/ga03186/kuaka-genome/05-scaffolding/
ASSEMBLY=/nesi/nobackup/ga03186/kuaka-genome/03-purge-dups/01-kuaka-hifiasm-p_ctg-purged.fa
APREFIX=kuaka_omnic
HIC_DIR=/nesi/project/ga03186/data/kuaka-HiC/
HIC_RAW1=${HIC_DIR}EXT041-07_L001_ds.75a4cbe42df64c0598d205c5a4e312df/Kuaka-OmniC-D206851_S1_L001_R
READ1=${HIC_RAW1}1_001.fastq.gz
READ2=${HIC_RAW1}2_001.fastq.gz
HIC_RAW2=${HIC_DIR}EXT041-07_L002_ds.dc52a69ac6504f0ca25e48595345a88b/Kuaka-OmniC-D206851_S1_L002_R
READ3=${HIC_RAW2}1_001.fastq.gz
READ4=${HIC_RAW2}2_001.fastq.gz
############

ml purge && module load fastp/0.23.2-GCC-11.3.0

### Clean HiC Reads with fastp.###
echo processing $READ1
fastp \
-i ${READ1} \
-o ${HIC_DIR}${APREFIX}_clean1_R1.fastq.gz \
-I ${READ2} \
-O ${HIC_DIR}${APREFIX}_clean1_R2.fastq.gz \
--trim_front1 15 \
--trim_front2 15 \
--qualified_quality_phred 20 \
--length_required 50 \
--thread 12

echo processing $READ3
fastp \
-i ${READ3} \
-o ${HIC_DIR}/${APREFIX}_clean2_R1.fastq.gz \
-I ${READ4} \
-O ${HIC_DIR}/${APREFIX}_clean2_R2.fastq.gz \
--trim_front1 15 \
--trim_front2 15 \
--qualified_quality_phred 20 \
--length_required 50 \
--thread 12

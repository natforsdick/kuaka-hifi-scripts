#!/bin/bash -e
#SBATCH -J bwa-HiFi
#SBATCH -A ga03186
#SBATCH --time=10:00:00 
#SBATCH --mem=6G  
#SBATCH --cpus-per-task=32
#SBATCH --out=%x.%j.out
#SBATCH --err=%x.%j.err
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=forsdickn@landcareresearch.co.nz
#SBATCH --profile=task

############################
# bwa-HiFi.sl
# Mapping HiFi to ref
# Nat Forsdick, 2021-11-02
############################

###########
# MODULES
module purge
module load  BWA/0.7.17-GCC-9.2.0 SAMtools/1.10-GCC-9.2.0
module list
###########

###########
# PARAMS

REFDIR=/nesi/nobackup/ga03186/kaki-hifi-asm/asm3-hic-hifiasm-p/01-purge-dups/
REFFILE=01P-asm3-hic-hifiasm-p-p_ctg-purged
REF=${REFDIR}${REFFILE}

INDIR1=/nesi/project/ga03186/data/JF_PacBio-kaki-Steeves-Order260/processed/
QUERY1=m54349U_210221_005741
OUTDIR1=/nesi/nobackup/ga03186/kaki-hifi-asm/asm3-hic-hifiasm-p/02-bwa/
###########

###########
# MKDIR
if [ ! -e ${OUTDIR1}  ]; then
    mkdir -p $OUTDIR1
    fi
###########

###########
# Index reference
#if [ ! -f ${REF}.fa.amb  ]; then
#    cd $REFDIR
#    echo "The reference has not been indexed. Indexing now"
# Added -b to speed up the indexing - tells it how large the block size should be when processing, and is apparently \
#    optimal when set to GenomeSize / 8. 
#    bwa index -a bwtsw ${REF}.fa -b 100000000
#    else
#        echo "BWA index file found" 
#        fi

###########
# Mapping
cd $OUTDIR1

echo Processing $QUERY1

echo "bwa mem -x pacbio -t 18 ${REF}.fa ${INDIR1}${QUERY1}.fastq"
bwa mem -x pacbio -t $SLURM_CPUS_PER_TASK ${REF}.fa ${INDIR1}${QUERY1}.fastq | samtools sort | samtools view -q 30 -Sb > ${QUERY1}_bwa.bam

samtools index ${QUERY1}_bwa.bam


# Now let's grab some mapping stats:
echo "Getting stats"
map=$(samtools view -F4 -c ${QUERY1}_bwa.bam)
unmap=$(samtools view -f4 -c ${QUERY1}_bwa.bam)
total=$(($map + $unmap))
perc_mapped=`echo "scale=4;($map/$total)*100" | bc`

echo "${QUERY1}_bwa.bam" >> ${QUERY1}_bwa_mapping_stats.txt
echo "mapped $map" >> ${QUERY1}_bwa_mapping_stats.txt
echo "% mapped $perc_mapped" >> ${QUERY1}_bwa_mapping_stats.txt
echo "unmapped $unmap" >> ${QUERY1}_bwa_mapping_stats.txt

echo completed $QUERY1


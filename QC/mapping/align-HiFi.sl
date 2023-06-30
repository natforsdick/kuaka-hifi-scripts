#!/bin/bash -e

#SBATCH --job-name=minimap-HiFi
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --time=01:30:00 #45m
#SBATCH --mem=56G # 6G
#SBATCH --ntasks=1
#SBATCH --profile=task 
#SBATCH --account=ga03186
#SBATCH --cpus-per-task=36
#SBATCH --profile=task

# align-HiFi.sl
# Align HiFi reads to a reference using minimap
# Nat Forsdick, 2021-10-28
 
# Based on subsamp of 60,000 HiFi reads taking 2 min to map (incl. initial indexing), 
# should take 45 min to map all HiFi data.

#########
# MODULES
module purge
module load minimap2/2.20-GCC-9.2.0 SAMtools/1.13-GCC-9.2.0
#########

#########
# PARAMS
#########
DIR=/nesi/nobackup/ga03186/kuaka-genome/05-mapping/
REFDIR=/nesi/nobackup/ga03186/kuaka-genome/03-purge-dups/
REF=01-kuaka-hifiasm-p_ctg-purged
HIFIDIR=/nesi/project/ga03186/data/kuaka-pacbio/AGRF_CAGRF22029575_PacBio/
HIFI=kuaka-hifi.fastq.gz
#########

cd $DIR

echo Indexing ${REF}

# To index reference genome the first time you run this - can then just call the index ref.mmi following this
minimap2 -t $SLURM_CPUS_PER_TASK -d ${REF}.mmi ${REFDIR}${REF}.fa

echo Aligning ${HIFI} against ${REF}
# To map HiFi reads to assembly
minimap2 -ax map-hifi -t $SLURM_CPUS_PER_TASK ${REF}.mmi ${HIFIDIR}${HIFI} | samtools sort -@ $SLURM_CPUS_PER_TASK -O BAM -o ${REF}-hifi.bam -

echo getting stats
samtools coverage ${REF}-hifi.bam -o ${REF}-hifi-cov.txt
samtools stats -@ 12 ${REF}-hifi.bam > ${REF}-hifi-stats.txt

echo plotting stats
plot-bamstats -p ${REF}-hifi-stats ${REF}-hifi-stats.txt

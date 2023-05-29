#!/bin/bash -e 
#SBATCH -A ga03186
#SBATCH -J hic_qc 
#SBATCH -c 32 # hicqc only needs 2 cpu
#SBATCH --mem=28G # hicqc requires <500M 
#SBATCH --partition=large 
#SBATCH --time=00:30:00 #hicqc only takes <5 min 
#SBATCH --mail-type=FAIL 
#SBATCH --mail-user=forsdickn@landcareresearch.co.nz 
#SBATCH --output %x.%j.out #  
#SBATCH --error %x.%j.err # 

################### 
# 02-hic-qc.sl
# Nat Forsdick, 2020-11-27
################### 

# Run Hi-C specific QC using the bam file produced from mapping Hi-C reads to assembly.
# This script takes two arguments, $1 : the path to the input directory, / 
# and $2 : the prefix of the input file. 
# e.g. to execute: 
# sbatch 02_hic_qc.sl /nesi/nobackup/ga03048/results/hic-qc-kaki/ Kaki_HiC_mapped 

################### 
# Need this version of miniconda
module purge
#module load Miniconda3/4.7.10
module load Miniconda3/4.10.3
ml SAMtools/1.13-GCC-9.2.0
################### 

#### ENVIRONMENT #### 
source /opt/nesi/CS400_centos7_bdw/Miniconda3/4.10.3/etc/profile.d/conda.sh
conda activate hic_qc

hic_qc=/nesi/nobackup/ga03048/modules/hic_qc/hic_qc.py
IN_DIR=$1
IN_BAM=$2
CPU=28
#####################

cd ${IN_DIR}

if [ ! -e ${IN_DIR}${IN_BAM}-sorted.bam ]; then
	echo "Sorting ${HIC}"
	samtools sort -n -@ $CPU ${IN_DIR}${IN_BAM}.bam -o ${IN_DIR}${IN_BAM}-sorted.bam
else
	echo "BAM index found, running QC"
fi


cd $IN_DIR
echo "running hic_qc for ${IN_BAM}"
python ${hic_qc} -b ${IN_BAM}-sorted.bam -o ${IN_BAM}-sorted.hicqc 
echo "finished running hic_qc for ${IN_BAM}"

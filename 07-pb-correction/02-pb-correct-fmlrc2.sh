#!/bin/bash -e

#SBATCH -A ga03048
#SBATCH -J fmlrc2
#SBATCH --time 1-00:00:00 # 
#SBATCH --mem 80G # 
#SBATCH --cpus-per-task 24 # 
#SBATCH	--error=%x.%A.%a.err
#SBATCH	--output=%x.%A.%a.out
#SBATCH --array=1-11%4
#SBATCH	--profile=task

# 02-pb-correct-fmlrc2.sh
# Nat Forsdick, 2022-04-07
# Step 2 for correcting raw pacbio CLR data
# Run 01-build-bwt-fmlrc2.sh first

###########
# PARAMS  #
datadir=/nesi/project/ga03048/weta/fastq/ # raw PacBio reads
illuminadir=/nesi/nobackup/ga03048/correction/trimmomatic/
outdir=/nesi/nobackup/ga03048/correction/fmlrc2/
###########

###########
# MODULES #
module purge
module load Anaconda3
source activate ropebwt2
module load rust-fmlrc/0.1.5-GCCcore-9.2.0
###########

echo "Beginning building BWT at"
date

cd ${datadir}

echo "Beginning fmlrc2 correction at"
date

#step 2: run fmlrc2 correction - $datadir holds raw PacBio reads

file=$(ls *.fastq | sed -n ${SLURM_ARRAY_TASK_ID}p)
echo $file

fmlrc2 -t $SLURM_CPUS_PER_TASK -C 14 ${outdir}weta_msbwt.npy ${file} ${outdir}${file}-fmlrc2-corr-reads.fa

echo "Completed fmlrc2 correction for ${file} at"
date

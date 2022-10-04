#!/bin/bash -e
#SBATCH --account ga03186
#SBATCH --job-name make-meryl
#SBATCH --time 00:10:00 # could need a couple of hours per fastq
#SBATCH --mem=40G # will prob need at least 24
#SBATCH --cpus-per-task=34
#SBATCH --array=1-2
#SBATCH --output %x.%A.%a.out
#SBATCH --error %x.%A.%a.err
#SBATCH --profile=task

##########
# PARAMS #
##########
k=20
genome=kaki
indir=/nesi/nobackup/ga03186/Kaki-HiSeq/trimmomatic/
outdir=/nesi/nobackup/ga03186/kaki-hifi-asm/asm-stats/merqury/
##########

cd $indir

# call samplist from file, and pass to array
SAMPLE_LIST=($(<input-fastq-list.txt))
SAMPLE=${SAMPLE_LIST[${SLURM_ARRAY_TASK_ID}+1]}

filename=$(basename "$SAMPLE") 
filename=${filename%.*} 
filename=${filename%.*}
echo $filename

meryl threads=26 memory=$SLURM_MEM_PER_NODE k=$k count output ${outdir}${filename}.meryl ${indir}${filename}.fastq.gz
date 

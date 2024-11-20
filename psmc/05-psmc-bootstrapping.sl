#!/bin/bash -e
#SBATCH -A ga03186
#SBATCH -J psmc-mapping
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH -t 15:00:00
#SBATCH --array=1-100%8
#SBATCH --out %x.%j.out
#SBATCH --err %x.%j.err

# run bootstrapping for psmc
# takes one parameter on start - the input *diploid.fq.gz file
# e.g., CDP_D206823_map_kuaka_ref_diploid.fq.gz
OUTDIR=/nesi/nobackup/ga03186/kuaka-genome/psmc/boots/
cd $OUTDIR

# set parameter space for psmc to the same as step 04

module purge
module load psmc/0.6.5-gimkl-2018b

#fastq=`cat ../samplist.txt | awk -v line=$SLURM_ARRAY_TASK_ID '{if(NR == line) print $0}'`
echo ieration $SLURM_ARRAY_TASK_ID

fastq=$1
filename=$(basename "$fastq")
filename=${filename%.fq.gz}

# perform bootstrapping
echo bootstrapping $fastq

xargs -i -n 1 -P 8 echo psmc -N30 -t5 -r5 -p "4+30*2+4+6+10" -o ${filename}_round-${SLURM_ARRAY_TASK_ID}.psmc ${filename}_split.psmcfa | sh

# combine results
#echo combining results for $fastq
#cat ../${filename}-2.psmc ${filename}_round-*.psmc > ${filename}_combined.psmc
#echo "completed bootstrapping for $filename"

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

OUTDIR=/nesi/nobackup/ga03186/kuaka-genome/psmc/boots/
cd $OUTDIR

# set parameter space for psmc
# params='-N30 -t5 -r5 -p "4+30*2+4+6+10"'

module purge
module load psmc/0.6.5-gimkl-2018b

#fastq=`cat ../samplist.txt | awk -v line=$SLURM_ARRAY_TASK_ID '{if(NR == line) print $0}'`
echo ieration $SLURM_ARRAY_TASK_ID

fastq=KWH_D206824_map_kuaka_ref_diploid.fq.gz
filename=$(basename "$fastq")
filename=${filename%.fq.gz}

# perform bootstrapping
echo bootstrapping $fastq
xargs -i -n 1 -P 8 echo psmc -N30 -t5 -r5 -p "4+30*2+4+6+10" -o ${filename}_round-${SLURM_ARRAY_TASK_ID}.psmc ${filename}_split.psmcfa | sh


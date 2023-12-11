#!/bin/bash -e
#SBATCH --account=ga03186
#SBATCH --job-name=sorting # job name (shows up in the queue)
#SBATCH --cpus-per-task=8 # mapping can use 18, subsequent processing requires 6
#SBATCH --mem=22G
#SBATCH --time=00:10:00 #Walltime (HH:MM:SS) # Total processing minimum 12 hrs
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=forsdickn@landcareresearch.co.nz
#SBATCH --output %x.%j.out # CHANGE number for new run
#SBATCH --error %x.%j.err #  CHANGE number for new run

#########
# PARAMS
PREFIX=01-kuaka-hifiasm-p_ctg-purged-clean-omnic
INDIR='/nesi/nobackup/ga03186/kuaka-genome/05-scaffolding/05b-Dovetail-OmniC/'
TMPDIR="/nesi/nobackup/ga03186/tmp-omnic-${SLURM_JOB_ID}"
CPU=16
########

cd $INDIR

echo sorting 
ml purge && ml SAMtools/1.15.1-GCC-11.3.0
samtools sort -@ $CPU ${PREFIX}_presort_marked.bam -o ${PREFIX}.bam

# index bam
echo indexing final bam
samtools index ${PREFIX}.bam
echo complete

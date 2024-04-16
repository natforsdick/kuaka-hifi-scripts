#!/bin/bash -e
#SBATCH --account=ga03186
#SBATCH --job-name=iggy-remdup # job name (shows up in the queue)
#SBATCH --cpus-per-task=8 # mapping can use 18, subsequent processing requires 6
#SBATCH --mem=4G
#SBATCH --time=01:00:00 #Walltime (HH:MM:SS) # Total processing minimum 12 hrs
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=forsdickn@landcareresearch.co.nz
#SBATCH --output %x.%j.out # CHANGE number for new run
#SBATCH --error %x.%j.err #  CHANGE number for new run

ml purge
ml SAMtools/1.15.1-GCC-11.3.0 BWA/0.7.17-GCC-11.3.0 pairtools/1.0.2-gimkl-2022a-Python-3.10.5

#########
# PARAMS
PREFIX=01-kuaka-hifiasm-p_ctg-purged-clean-omnic
INDIR='/nesi/nobackup/ga03186/kuaka-genome/05-scaffolding/05b-Dovetail-OmniC/'
TMPDIR="/nesi/nobackup/ga03186/tmp-omnic-${SLURM_JOB_ID}"
CPU=16
########

cd $INDIR
mkdir $TMPDIR

### Flag PCR Duplicates with SAMBLASTER #######################################
echo convert bam to sam
ml purge && ml SAMtools/1.15.1-GCC-11.3.0
samtools sort -n ${PREFIX}-mapped.PT.bam -o ${PREFIX}-mapped.PT.sam

echo flagging dups
ml purge && ml load samblaster/0.1.26-GCC-9.2.0
samblaster -i ${PREFIX}-mapped.PT.sam -o ${TMPDIR}/${PREFIX}_marked_byread.sam

ml purge && ml SAMtools/1.15.1-GCC-11.3.0
### Remove unmmaped and non-primary aligned reads. Sort and Index bam files.###
echo removing unmapped/non-primary reads
samtools view -S -b -h -@ $CPU -F 2316 ${TMPDIR}/${PREFIX}_marked_byread.sam > ${PREFIX}_presort_marked.bam
echo complete

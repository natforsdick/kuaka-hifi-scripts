#!/bin/bash -e
#SBATCH --account=ga03186
#SBATCH --job-name=omnic-map # job name (shows up in the queue)
#SBATCH --cpus-per-task=36 # mapping can use 18, subsequent processing requires 6
#SBATCH --mem=48G
#SBATCH --time=2-00:00:00 #Walltime (HH:MM:SS) # Total processing minimum 12 hrs
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=forsdickn@landcareresearch.co.nz
#SBATCH --output %x.%j.out # CHANGE number for new run
#SBATCH --error %x.%j.err #  CHANGE number for new run

ml purge
ml SAMtools/1.15.1-GCC-11.3.0 BWA/0.7.17-GCC-11.3.0 pairtools/1.0.2-gimkl-2022a-Python-3.10.5

#########
# PARAMS
PREFIX=01-kuaka-hifiasm-p_ctg-purged-clean-omnic
INDIR='/nesi/nobackup/ga03186/kuaka-genome/05-scaffolding/05b-Dovetail-OmniC/' 
OMNICR1=/nesi/project/ga03186/data/kuaka-HiC/kuaka_omnic_clean_R1.fastq.gz
OMNICR2=/nesi/project/ga03186/data/kuaka-HiC/kuaka_omnic_clean_R2.fastq.gz
REF=01-kuaka-hifiasm-p_ctg-purged
TMPDIR="/nesi/nobackup/ga03186/tmp-omnic-${SLURM_JOB_ID}"
CPU=32
########

cd $INDIR
mkdir $TMPDIR

# alignment -T0 = all alignments to generate stats, quality filtering later
#echo aligning
#bwa mem -5SP -T0 -t $CPU $REF.fa $OMNICR1 $OMNICR2 -o ${PREFIX}-aligned.sam

# find ligation junctions
#echo finding ligation junctions
#pairtools parse --min-mapq 40 --walks-policy 5unique --max-inter-align-gap 30 --nproc-in $CPU \
#--nproc-out $CPU --chroms-path ${REF}.genome ${PREFIX}-aligned.sam > ${PREFIX}-parsed.pairsam

# sort pairsam
echo sorting pairsam
pairtools sort --nproc $CPU --tmpdir=$TMPDIR ${PREFIX}-parsed.pairsam > ${PREFIX}-sorted.pairsam

# remove duplicates
echo removing duplicates
pairtools dedup --nproc-in $CPU --nproc-out $CPU --mark-dups --output-stats ${PREFIX}-stats.txt \
--output ${PREFIX}-dedup.pairsam ${PREFIX}-sorted.pairsam

# split .bam, .pairs
echo splitting bam
pairtools split --nproc-in $CPU --nproc-out $CPU --output-pairs ${PREFIX}-mapped.pairs \
--output-sam ${PREFIX}-unsorted.bam ${PREFIX}-dedup.pairsam

# sort bam
echo sorting bam
samtools sort -@${CPU} -T ${TMPDIR}tempfile.bam -o ${PREFIX}-mapped.PT.bam ${PREFIX}-unsorted.bam

# index bam
echo indexing final bam
samtools index ${PREFIX}-mapped.PT.bam

if [ -f ${PREFIX}-mapped.PT.bam ]
then
echo "pipeline completed"
else
echo "pipeline not complete"
fi


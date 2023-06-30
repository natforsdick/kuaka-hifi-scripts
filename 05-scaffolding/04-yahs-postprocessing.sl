#!/bin/bash
#SBATCH --account=ga03186
#SBATCH --job-name=make-juice-in # job name (shows up in the queue)
#SBATCH --cpus-per-task=4
#SBATCH --mem=26G
#SBATCH --time=01:30:00 #Walltime (HH:MM:SS)
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=forsdickn@landcareresearch.co.nz
#SBATCH --output %x.%j.out # CHANGE number for new run
#SBATCH --error %x.%j.err 

# Generating .hic file for visualisation and manual curation

##########
# PARAMS #
YAHSJUICE=/nesi/project/ga03186/scripts/Hi-C_scripts/yahs/juicer
OUTDIR=/nesi/nobackup/ga03186/kuaka-genome/05-scaffolding/05b-Dovetail-OmniC/yahs/
JUICER=/nesi/nobackup/ga03048/juicer/scripts/juicer_tools.1.9.9_jcuda.0.8.jar
REF_DIR=/nesi/nobackup/ga03186/kuaka-genome/03-purge-dups/
REF=01-kuaka-hifiasm-p_ctg-purged.fa
REFPRE=01-kuaka-hifiasm-p_ctg-purged-DT-yahs
##########

cd $OUTDIR

echo 'generating hic contact map'
$YAHSJUICE pre -a -o out_JBAT ${REFPRE}.bin ${REFPRE}_scaffolds_final.agp \
${REF_DIR}${REF}.fai > out_JBAT.log 2>&1
echo done step 1

echo 'running juicer_tools pre'
java -jar -Xmx26G $JUICER pre out_JBAT.txt out_JBAT.hic.part \
	<(cat out_JBAT.log | grep PRE_C_SIZE | awk '{print $2" "$3}')

mv out_JBAT.hic.part out_JBAT.hic

echo done step 2



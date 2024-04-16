#!/bin/bash -e

#SBATCH --account=ga03186
#SBATCH --job-name=fcsgx
#SBATCH --time=3-00:00:00
#SBATCH --gres=ssd
#SBATCH --partition=milan 
#SBATCH --mem=6G
#SBATCH --cpus-per-task=1
#SBATCH --hint=nomultithread
#SBATCH --output %x.%j.out
#SBATCH --mail-user=forsdickn@landcareresearch.co.nz
#SBATCH --mail-type=FAIL,END
#SBATCH --array=0-8

REFDIR=/nesi/nobackup/ga03186/kuaka-genome/05-scaffolding/05b-Dovetail-OmniC/all-data-yahs/
TAXID=2884279
OUTDIR=/nesi/nobackup/ga03186/kuaka-genome/contam-screening/
FCS=/nesi/project/landcare03691/fcs-gx/

# prior to running, make your filelist from your chunked fasta input:
# ls ${ref}.fa.[0-9] > ${REFDIR}filelist.txt
FILES=($(<${REFDIR}filelist.txt))
FILENAME=${FILES[${SLURM_ARRAY_TASK_ID}]}
echo $FILENAME

module purge
module load Python/3.8.2-gimkl-2020a

#export dist and scripts to dist
export PATH=/nesi/project/landcare03691/fcs-gx/dist:$PATH
export PATH=/nesi/project/landcare03691/fcs-gx/scripts:$PATH
export TMPDIR=/nesi/nobackup/ga03048/fcsgx

#Assign the number of CPUS by using the Slurm variable
export GX_NUM_CORES=12

#Path to local Database
export LOCAL_DB=/nesi/nobackup/nesi02659/LRA/resources/fcs/gxdb

#This is an attempt to cache the database to milan SSD which is pointing to $TMPDIR
#sync_files.py get --mft=$LOCAL_DB/all.manifest --dir=$TMPDIR/gxdb

export GXDB_LOC=$TMPDIR

#Check the integrity of the local database vs copied one to $TMPDIR
#python3 ${FCS}fcs.py db check --mft "$LOCAL_DB/all.manifest" --dir $TMPDIR/gxdb

mkdir -p ${OUTDIR}${FILENAME}
cd ${OUTDIR}${FILENAME}
#Run GX
echo running GX module for $FILENAME
run_gx --fasta ${REFDIR}${FILENAME} --tax-id $TAXID --gx-db "$GXDB_LOC/gxdb" --out-dir ${OUTDIR}${FILENAME}

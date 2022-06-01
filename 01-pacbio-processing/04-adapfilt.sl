#!/bin/bash -e 
#SBATCH -A ga03186 # CHANGE!! 
#SBATCH -J adapfilt # job name (shows up in the queue) 
#SBATCH -c 8
#SBATCH --mem=2G 
#SBATCH --time=02:30:00 #Walltime (HH:MM:SS) 
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err

module purge
module load BLAST/2.12.0-GCC-9.2.0 BamTools/2.5.2-GCC-11.3.0

export PATH=$PATH:/nesi/project/ga03048/modules/HiFiAdapterFilt
export PATH=$PATH:/nesi/project/ga03048/modules/HiFiAdapterFilt/DB

INDIR=/nesi/nobackup/ga03186/kuaka-genome/01-preprocessing/
HIFI=kuaka-hifi
CPU=12

cd $INDIR
echo "Beginning adapter filtering"
date
bash /nesi/project/ga03048/modules/HiFiAdapterFilt/pbadapterfilt.sh -p ${HIFI} -t $CPU

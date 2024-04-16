#!/bin/bash -e

#SBATCH -A ga03186
#SBATCH -J BUSCO # job name (shows up in the queue)
#SBATCH -c 8
#SBATCH --mem=6G
#SBATCH --partition=large
#SBATCH --time=01:00:00 #Walltime (HH:MM:SS)
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=forsdickn@landcareresearch.co.nz
#SBATCH --output %x.%j.out #
#SBATCH --error %x.%j.err #

bash ./BUSCO1.sh

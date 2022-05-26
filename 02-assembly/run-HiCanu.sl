#!/bin/bash -e
#SBATCH --account=ga03186
#SBATCH --job-name=hicanu 
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem=100M
#SBATCH --partition=large
#SBATCH --time=00:05:00 
#SBATCH --mail-type=ALL
#SBATCH --mail-user=natalie.forsdick@postgrad.otago.ac.nz
#SBATCH --output %x.%j.out 
#SBATCH --error %x.%j.err 

module purge
module load Canu/2.1.1-GCC-9.2.0
ml list

date

# -d = assembly directory
canu -p kaki-asm -d /nesi/nobackup/ga03186/kaki-hifi-asm genomeSize=1.2g useGrid=true gridOptions="--account=ga03186 --time=01:00:00" -pacbio-hifi /nesi/project/ga03186/data/JF_PacBio-kaki-Steeves-Order260/processed/m54349U_210221_005741.fastq

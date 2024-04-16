#!/bin/bash -e
#SBATCH -A ga03186
#SBATCH -J compleasm # job name (shows up in the queue)
#SBATCH -c 22
#SBATCH --mem=14GB #
#SBATCH --time=01:30:00 #Walltime (HH:MM:SS)
#SBATCH --output %x.%j.out #
#SBATCH --error %x.%j.err #

# compleasm - 2023-06-15
# Genome assembly QC - BUSCO alternative

# PARAMS
INDIR=/nesi/nobackup/ga03186/kuaka-genome/ref-genomes/
OUTDIR=/nesi/nobackup/ga03186/kuaka-genome/ref-genomes/
DB=/nesi/project/ga03186/aves_odb10
ASM=GCA_013400755.1_ASM1340075v1_genomic.fna

asm=$(basename $ASM .fna)

cd $INDIR
module purge
module load compleasm/0.2.2-gimkl-2022a

compleasm.py run -a ${INDIR}${ASM} -o ${OUTDIR}compleasm-${asm} -t 16  -l aves_odb10 -L $DB

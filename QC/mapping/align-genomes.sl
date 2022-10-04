#!/bin/bash -e

#SBATCH --job-name=minimap-aln
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --time=00:12:00
#SBATCH --mem=12G #Used 26 for whole genomes
#SBATCH --ntasks=1
#SBATCH --profile=task 
#SBATCH --account=ga03186
#SBATCH --cpus-per-task=8 # Used 32 for whole genomes

# align-genomes.sl
# Align genomes to one another using minimap
# Nat Forsdick, 2021-09-01

#########
# MODULES
module purge
module load minimap2/2.20-GCC-9.2.0 
#########

INDIR=/nesi/nobackup/ga03186/kaki-hifi-asm/asm3-hic-hifiasm-p/03-polishing/
QUERY=01P-asm3-hic-hifiasm-p-p_ctg-purged.cns
REFDIR=/nesi/nobackup/ga03186/kaki-hifi-asm/asm3-hic-hifiasm-p/01-purge-dups/
REF=01P-asm3-hic-hifiasm-p-p_ctg-purged

cd $INDIR

echo "Aligning $QUERY against reference genome"

# To index reference genome the first time you run this - can then just call the index ref.mmi following this
if [ ! -e ${REF}.mmi ]; then
	echo "Index file of reference does not exist: creating index"
	minimap2 -t $SLURM_CPUS_PER_TASK -d ${REF}.mmi ${REFDIR}${REF}.fa
fi
echo "Aligning $QUERY to $REF"
minimap2 -ax asm5 -t $SLURM_CPUS_PER_TASK ${REF}.mmi ${INDIR}${QUERY}.fa > ${REF}-${QUERY}.paf

echo "Collecting stats"
$HOME/bin/k8 /nesi/project/ga03186/HiFi-scripts/paftools.js stat ${REF}-${QUERY}.paf > ${REF}-${QUERY}-stat.out

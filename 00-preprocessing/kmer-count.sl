#!/bin/bash -e
#SBATCH -A ga03186
#SBATCH -J jellyfish
#SBATCH --time=00:40:00
#SBATCH -c 36
#SBATCH --mem=40G
#SBATCH --mail-user=forsdickn@landcareresearch.co.nz
#SBATCH --mail-type=FAIL
#SBATCH --output %x.%j.out
#SBATCH --error %x.%j.err
#SBATCH --profile=task

# kmer-count.sl
# 2022-05-27, Nat Forsdick
# Counting kmers for genome size estimation.
# Can use GenomeScope web assessment to estimate genome size.
# Requires jellyfish kmer histogram report as input.
# https://genome.umd.edu/docs/JellyfishUserGuide.pdf

##########
# PARAMS
DATADIR=/nesi/nobackup/ga03186/frog-genome/01-preprocessing/
READS=frog-hifi-cell1
mer=21
CPU=20
##########

ml purge
ml Jellyfish/2.3.0-gimkl-2020a
cd ${DATADIR}

# To generate: -C = canonical (necessary), -m = kmers, -s = hash size, -t = threads

echo 'Beginning jellyfish bloom'
date
jellyfish bc -C -m ${mer} -s 5G -t ${CPU} -o ${READS}-${mer}mer.bc ${READS}.fastq

echo 'Beginning jellyfish count'
jellyfish count -C -m ${mer} -s 5G -t ${CPU} --bc ${READS}-${mer}mer.bc -o ${READS}-${mer}mer-counts.jf ${READS}.fastq

echo 'Beginning jellyfish histo'
date
jellyfish histo -t ${CPU} ${READS}-${mer}mer-counts.jf > ${READS}-${mer}mer-histo.out

echo 'Completed kmer-count'
date

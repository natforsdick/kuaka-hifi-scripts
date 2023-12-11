#!/bin/bash -e

#SBATCH --account ga03186
#SBATCH --job-name make-fasta
#SBATCH --cpus-per-task=4
#SBATCH --mem 6G
#SBATCH --time 00:15:00 # When starting this job, set to 04:00:00, then change to 00:30:00 for final mem-intensive step
#SBATCH --output %x.%j.out
#SBATCH --error %x.%j.err

OUTDIR=/nesi/nobackup/ga03186/kuaka-genome/asm2-hifiasm-hic/
OUTPRE=kuaka-hifiasm2.asm.hic

cd $OUTDIR
awk '/^S/{print ">"$2;print $3}' ${OUTPRE}.p_ctg.gfa > ${OUTPRE}.p_ctg.fa
awk '/^S/{print ">"$2;print $3}' ${OUTPRE}.a_ctg.gfa > ${OUTPRE}.a_ctg.fa

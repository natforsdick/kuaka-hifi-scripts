#!/bin/bash -e

OUTDIR=/nesi/nobackup/ga03186/kuaka-genome/psmc/boots/
cd $OUTDIR

# combine results of bootstrapping
samplist=*diploid.psmc

ml purge; ml psmc

for fastq in $samplist
do
    echo combining results for $fastq
    cat ../${filename}.psmc ${filename}_round-*.psmc > ${filename}_combined.psmc
    echo plotting
    psmc_plot.pl -u 7.38e-9 -g 6 -R ${filename}_boot ${filename}_combined.psmc
    psmc_plot.pl -u 7.38e-9 -g 6 -R ${filename}_single ${filename}*psmc.psmc # need to check this input filename
done

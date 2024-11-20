#!/bin/bash -e
#SBATCH -A ga03186
#SBATCH -J psmc
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G # 10 when including psmcfa step
#SBATCH -t 2:20:00 # 3 when including psmcfa step
#SBATCH --out %x.%j.out
#SBATCH --err %x.%j.err

# converting consensus to psmc format, and testing parameter space

OUTDIR=/nesi/nobackup/ga03186/kuaka-genome/psmc/
cd $OUTDIR

module purge
module load psmc/0.6.5-gimkl-2018b

#for fastq in *_diploid.fq.gz
for fastq in KWH_D206824_map_kuaka_ref_diploid.fq.gz
do
    filename=$(basename "$fastq") 
    filename=${filename%.fq.gz}
    # Let's convert the diploid genome to PSMC suitable format
    if [ ! -e  ${filename}.psmcfa ]; then  
        echo "converting $filename"  
        fq2psmcfa -q20 ${fastq} > ${filename}.psmcfa
    else
        echo "$filename psmcfa file exists"
    fi

    # Running PSMC with standard parameters to get initial idea about how the parameters are working.
    echo "running psmc for $filename"
# -N max iterations, -t max 2N0 coalescent time, r initial theta/rho ratio - based on those used in Nadachowska-Brzyska et al 2015
    psmc -N30 -t5 -r5 -p "1+3+30*2+6+6" -o ${filename}-v3.psmc ${filename}.psmcfa

    # split to run smaller chunks for bootstrapping
#    echo splitting $fastq
#    splitfa ${filename}.psmcfa > ./boots/${filename}_split.psmcfa
done

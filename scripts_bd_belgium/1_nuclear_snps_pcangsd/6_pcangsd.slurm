#!/bin/bash -l
#SBATCH -A lp_kmma
#SBATCH -M genius
#SBATCH -p bigmem
#SBATCH --mem-per-cpu=19000mb
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=36
#SBATCH --job-name=pcangsd
#SBATCH --output=pcangsd.log
#SBATCH --error=pcangsd.error
ml GSL/2.7-GCC-12.3.0

cd $VSC_SCRATCH/angsd/bd_intercept_jan_miss90_dp2dp15_properpairs_q30
INPUT=bd_intercept.sorted.beagle.gz #must be .gz

pcangsd -b $INPUT -t $SLURM_CPUS_PER_TASK --maf 0.05 --iter 5000 --tree --tree-samples sample_ids.txt -o bd_intercept


#!/bin/bash -l
#SBATCH -A lp_kmma
#SBATCH -M wice
#SBATCH -p hugemem
#SBATCH --time=6:00:00
#SBATCH --mem-per-cpu=100000mb
#SBATCH --cpus-per-task=8
#SBATCH --array=35
#SBATCH --job-name=angsd
#SBATCH --output=angsd_bd_%A_%a.log
#SBATCH --error=angsd_bd_%A_%a.error
ml GSL/2.7-GCC-12.3.0

# Path to the regions file
REGIONS_FILE=$VSC_DATA/bdrefgenome/intervals.angsd

# Read the region for this task
I=$SLURM_ARRAY_TASK_ID
REGION=$(sed -n "$((I+1))p" $REGIONS_FILE) #print line I + 1

# Path to the genome
GENOME=$VSC_DATA/bdrefgenome/GCF_023373825.1_ASM2337382v1_genomic.fna

# Path to directory with BAM files
ALIDIR=$VSC_SCRATCH/bam/bdglobal/bam

# Path to output directory
OUTDIR=$VSC_SCRATCH/angsd/bd_intercept_jan_miss90_dp2dp15_properpairs_q30

# Create output directory
mkdir -p $OUTDIR >& /dev/null

# Log for this region
LOG=$OUTDIR/angsd-${REGION}.log

# Get list of bam files on a temporary file
BAMLIST=$(mktemp)
ls $ALIDIR/*.bam > $BAMLIST

# Call SNPs and genotypes with ANGSD
if [ -f $OUTDIR/angsd-${REGION}.beagle.gz ]; then
    echo "Output file for" $REGION "already exists. Skipping."
else

angsd \
    -nThreads $SLURM_CPUS_PER_TASK \
    -bam $BAMLIST \
    -ref $GENOME \
    -r $REGION \
    -out $OUTDIR/angsd-${REGION} \
    -minMapQ 30 -minQ 30 -minInd 486 -setMinDepth 2160 -setMinDepthInd 2 -setMaxDepth 8100 -setMaxDepthInd 15 -minmaf 0.05 \
    -baq 1 -uniqueOnly 1 -remove_bads 1 -SNP_pval 1e-6 \
    -doMaf 3 -doMajorMinor 1 -gl 1 -doPlink 2 -doGeno -4 --ignore-RG 0 -doPost 1 -doCounts 1 -dumpcounts 1 -dobcf 1 -doDepth 0 \
    -doGlf 2 -only_proper_pairs 1  # Output in Beagle format
fi


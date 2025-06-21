#!/bin/bash -l
#SBATCH -A lp_kmma
#SBATCH -M wice
#SBATCH -p batch
#SBATCH --job-name=mosdepth
#SBATCH --output=logs/mosdepth_%A_%a.out
#SBATCH --error=logs/mosdepth_%A_%a.err
#SBATCH --array=100-137
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=3400mb
#SBATCH --time=0:10:00

# Read the BAM file corresponding to this SLURM array task
BAM_DIR=$VSC_SCRATCH/bam
cd $BAM_DIR
BAM_FILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" bam_list.txt)
PREFIX=$(basename "$BAM_FILE" .bam)  # Extract sample name as prefix

# Check if the BAM file exists before running
if [[ ! -f "$BAM_FILE" ]]; then
    echo "Error: BAM file $BAM_FILE not found" >&2
    exit 1
fi

# Run mosdepth in fast mode (keeps per-site depth)

if [[ -f $PREFIX.flagstat.txt ]]; then
    echo "Output for $BAM_FILE already present" >&2
else
   samtools index "$BAM_FILE"
   mosdepth --fast-mode "$PREFIX" "$BAM_FILE"
   samtools flagstat "$BAM_FILE" > "$PREFIX".flagstat.txt
fi

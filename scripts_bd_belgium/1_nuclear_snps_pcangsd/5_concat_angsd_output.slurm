#!/bin/bash -l
#SBATCH -A lp_kmma
#SBATCH -M genius
#SBATCH -p batch
#SBATCH --time=03:00:00
#SBATCH --mem-per-cpu=3400mb
#SBATCH --cpus-per-task=36
#SBATCH --job-name=concat_beagle
#SBATCH --output=concat_beagle.log
#SBATCH --error=concat_beagle.error

ml GSL/2.7-GCC-12.3.0
# Path to input files !!!!!
INPUT_DIR=$VSC_SCRATCH/angsd/bd_intercept_jan_miss90_dp2dp15_properpairs_q30


GENOME=$VSC_DATA/bdrefgenome/GCF_023373825.1_ASM2337382v1_genomic.fna
REGIONS_FILE=$VSC_DATA/bdrefgenome/intervals.angsd


# Define output directories
OUTDIR=$VSC_SCRATCH/angsd/bd_intercept_jan_miss90_dp2dp15_properpairs_q30
PREFIX=bd_intercept


# Create necessary directories
mkdir -p $OUTDIR
cd $OUTDIR

BAMLIST=$(mktemp)
ls $ALIDIR/*.bam > $BAMLIST

TPED=$PREFIX.tped
TFAM=$PREFIX.sorted.tfam
BEAGLE=$PREFIX.beagle
MAF=$PREFIX.mafs
POS=$PREFIX.pos

BASE=$(basename $BEAGLE .beagle)
SORTEDBEAGLE=$BASE.sorted.beagle

# Function to extract columns with the same header
xtractcol() {
  local file=$1
  shift
  cut -f$(head -n1 "$file" | tr \\t \\n | grep -Fxn "${@/#/-e}" | cut -d: -f1 | paste -sd,) "$file"
}
for file in $INPUT_DIR/angsd*bcf; do
  bcftools index -f --threads $SLURM_CPUS_PER_TASK $file
done

# Concatenate bcf files
bcftools concat -a -d all --threads $SLURM_CPUS_PER_TASK -Ou -o $PREFIX.bcf $INPUT_DIR/angsd*bcf
bcftools index $PREFIX.bcf

rm $BEAGLE


# Concatenate regional Beagle files
first=true
while IFS= read -r line; do
  if $first; then
    # For the first line, decompress without using tail
    zcat "$INPUT_DIR/angsd-$line.beagle.gz" >> "$BEAGLE"
    bcftools query -l "$INPUT_DIR/angsd-$line.bcf" > sample_ids.txt
    first=false
  else
    # For subsequent lines, skip the header line
    zcat "$INPUT_DIR/angsd-$line.beagle.gz" | tail -n +2 >> "$BEAGLE"
  fi
done < "$REGIONS_FILE"

# Sort concatenated Beagle and remove duplicate lines
# Extract the header (first line)
cat $BEAGLE | head -n 1 > $OUTDIR/$BASE.sorted.beagle
cat $BEAGLE | sort --parallel $SLURM_CPUS_PER_TASK -k 1 -V -u | grep -v marker >> $OUTDIR/$BASE.sorted.beagle
# Remove uncompressed unsorted Beagle
rm $BEAGLE

pigz $OUTDIR/$BASE.sorted.beagle

exit

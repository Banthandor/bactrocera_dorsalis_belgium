#!/bin/bash -l
#SBATCH -A lp_kmma
#SBATCH -M wice
#SBATCH -p batch
#SBATCH --job-name=minimap
#SBATCH -o minimap_array_%a.out
#SBATCH --error minimap_array_%a.error
#SBATCH --time=1:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=72
#SBATCH --mem-per-cpu=3400mb
#SBATCH --array=0-99

# Define input and output directories

output_dir=$VSC_SCRATCH/bam/tephritidae
FILES=($VSC_SCRATCH/fastq/batch_oct_2/*R1.fastq.gz)
cd $VSC_SCRATCH/fastq_tephritidae
FILE=${FILES[$SLURM_ARRAY_TASK_ID]}
OUTFILE=$output_dir/$(basename ${FILE} .R1.fastq.gz).markdup.bam
base=$(basename $FILE .R1.fastq.gz)
echo ${FILE}

conda activate samtools

REF=$VSC_DATA/bdrefgenome/GCF_023373825.1_ASM2337382v1_genomic.fna

# Create output directory if it doesn't exist

mkdir -p "$output_dir"

# Check if the output file already exists
if [ -f $OUTFILE ]; then
    echo "Output file already exists. Skipping."
else
    echo "mapping $FILE"
    echo "${SLURM_ARRAY_TASK_ID}"
    # Run minimap2 to map the FASTQ file to the reference genome, convert to BAM, and sort

    srun fastp --in1 $FILE --in2 "$base".R2.fastq.gz --qualified_quality_phred 20 --correction --length_required 75 --stdout --stdin | \
               /vsc-hard-mounts/leuven-data/356/vsc35633/minimap2-2.28_x64-linux/minimap2 -t 36 -R '@RG\tID:'$base'\tSM:'$base -ax sr "$REF" - | \
               samtools sort -@18 -m3300m -n -u - | \
               samtools fixmate -m -u - - | \
               samtools sort -@18 -m3300m - | \
               samtools markdup --output-fmt BAM -s - $OUTFILE

    # Index the marked duplicates BAM file
    samtools index $OUTFILE
fi

echo "Mapping and duplicate marking completed successfully."
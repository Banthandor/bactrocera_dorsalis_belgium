#!/bin/bash -l
#SBATCH -A lp_kmma
#SBATCH -M wice
#SBATCH -p batch
#SBATCH --job-name=subsample
#SBATCH --output=subsample.out
#SBATCH --error=subsample.err
#SBATCH --cpus-per-task=36
#SBATCH --mem-per-cpu=3400mb
#SBATCH --time=1:00:00

cd $VSC_SCRATCH/bam/

# Step 1: Extract sample IDs and depth values from *mosdepth.summary.txt files
output_file="depth_values.txt"

# Write header to the output file
echo -e "Sample\tDepth" > "$output_file"

# Process each *mosdepth.summary.txt file
for file in *mosdepth.summary.txt; do
    # Extract the sample ID (everything before the first '.')
    sample_id=$(basename "$file" | cut -d '.' -f1)

    # Extract the depth value from the last line, fourth column
    depth=$(tail -1 "$file" | cut -f4)

    # Save the sample ID and depth value to the output file
    echo -e "$sample_id\t$depth" >> "$output_file"
done

echo "Sample IDs and depth values extracted and saved to $output_file."

# Step 2: Calculate downsampling fractions
input_file="depth_values.txt"
downsampling_file="downsampling_fractions.txt"

# Target depth
target_depth=10

# Write header to the downsampling file
echo -e "Sample\tDepth\tDownsampling_Fraction" > "$downsampling_file"

# Read the depth values file line by line (skip the header)
tail -n +2 "$input_file" | while IFS=$'\t' read -r sample depth; do
  # Calculate the downsampling fraction
  downsampling_fraction=$(echo "scale=3; $target_depth / $depth" | bc)

  # Append the result to the downsampling file
  echo -e "$sample\t$depth\t$downsampling_fraction" >> "$downsampling_file"
done

echo "Downsampling fractions calculated and saved to $downsampling_file."

# Print the results to the terminal
echo "Contents of $downsampling_file:"
cat "$downsampling_file"

# Step 3: Downsample samples with depth above 14
echo "Downsampling samples with depth above 14..."

# Create a directory to store original BAM files
original_bam_dir="original_bam_files"
mkdir -p "$original_bam_dir"

# Read the downsampling fractions file line by line (skip the header)
tail -n +2 "$downsampling_file" | while IFS=$'\t' read -r sample depth downsampling_fraction; do
  # Check if the depth is above 14
  if (( $(echo "$depth > 14" | bc -l) )); then
    # Define input BAM file (check for both ${sample}.bam and ${sample}.markdup.bam)
    if [[ -f "${sample}.bam" ]]; then
      input_bam="${sample}.bam"
    elif [[ -f "${sample}.markdup.bam" ]]; then
      input_bam="${sample}.markdup.bam"
    else
      echo "Error: BAM file for sample $sample not found (neither ${sample}.bam nor ${sample}.markdup.bam exists). Skipping."
      continue
    fi

    # Define output BAM file
    output_bam="${sample}_downsampled.bam"

    # Check if the output BAM file already exists
    if [[ -f "$output_bam" ]]; then
      echo "Output file $output_bam already exists. Skipping downsampling for sample $sample."
      continue
    fi

    # Downsample the BAM file using samtools subsample
    echo "Downsampling $input_bam to $output_bam with fraction $downsampling_fraction..."
    samtools view --threads $SLURM_CPUS_PER_TASK --subsample "$downsampling_fraction" -b -o "$output_bam" "$input_bam"

    # Index the downsampled BAM file
    samtools index "$output_bam"

    echo "Downsampled $input_bam to $output_bam and indexed."

    # Move the original BAM file to the new directory
    echo "Moving original BAM file $input_bam to $original_bam_dir..."
    #mv "$input_bam" "$original_bam_dir/"
  else
    echo "Skipping $sample (depth = $depth, not above 14)."
  fi
done

echo "Downsampling completed for samples with depth above 14."
echo "Original BAM files moved to $original_bam_dir."
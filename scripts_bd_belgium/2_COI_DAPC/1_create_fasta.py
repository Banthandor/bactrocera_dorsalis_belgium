import pandas as pd
import subprocess
from Bio import SeqIO

# Define the input and output file names, the input file is a tsv file from bold containing both sequence and metadata
input_file = "bactrocera_dorsalis_bold_acc_15012025.tsv"
output_fasta = "bactrocera_dorsalis_coi_sequences_all.fasta"

# Load the tab-separated data into a pandas DataFrame
df = pd.read_csv(input_file, sep='\t')

# Step 1: Generate fasta
with open(output_fasta, 'w') as fasta_out:
    for index, row in df.iterrows():
        if pd.isna(row['nuc']) or row['nuc'] == 'nan' or row['nuc'].strip() == "":
            print(f"Warning: Missing or invalid sequence data for {row['processid']}, skipping entry.")
            continue
        header = f">{row['processid']} [Species: {row['species']}] [Country: {row['country/ocean']}]"
        fasta_out.write(f"{header}\n")
        fasta_out.write(f"{row['nuc']}\n")

print(f"FASTA file '{output_fasta}' has been generated.")
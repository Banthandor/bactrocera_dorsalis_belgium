# Genomic Tracing of *Bactrocera dorsalis* in Belgium

This repository contains all scripts and core output files used in the analyses for the manuscript:

> **"Genomic tracing reveals multiple independent occurrences of *Bactrocera dorsalis* in Belgium"**

---

## 📁 Repository Structure

### `scripts/1_nuclear_snps_pcangsd/`
SLURM scripts and output files for population structure analysis based on nuclear SNPs using ANGSD and PCAngsd.

- `1_fastp_minimap.slurm` – raw read processing and mapping
- `2_mosdepth.slurm` – depth of coverage estimation of bam files
- `3_mosdepth_downsample.slurm` – Downsampling
- `4_angsd.slurm` – Genotype likelihood estimation
- `5_concat_angsd_output.slurm` – Merging ANGSD outputs
- `6_pcangsd.slurm` – PCAngsd inference
- `7_pcangsd_output/` – Tree and covariance matrix from PCAngsd

### `scripts/2_COI_DAPC/`
Scripts and data for COI haplotype-based Discriminant Analysis of Principal Components (DAPC).

- `1_create_fasta.py` – Generates FASTA from BOLD data + metadata file
- `2_dapc.R` – DAPC model construction
- `3_predict.R` – Posterior cluster prediction for new query samples
- `bactrocera_dorsalis_coi_bold_wgs_dapc.fasta` – filtered COI alignment used for DAPC

### `scripts/3_mitogenome_tree/`
Pipeline to reconstruct and analyze a phylogeny based on full mitochondrial genomes.

- `1_reconstruct_mitogenome.slurm` – mapping reads to mitogenome reference and consensus sequence generation
- `2_mafft_fasttree.slurm` – Multiple alignment and tree inference
- `bd.mitogenomes.fasta` – Input sequences
- `bd.mitogenomes.mafft.filtered.trimmed.fasta` – Filtered alignment
- `*.tree.newick` – Phylogenetic tree of mitogenomes

---

## 📄 How to Use

This repository is meant for transparency and reproducibility. The scripts are designed to run in an HPC environment (SLURM-based), and assume access to input data described in the manuscript.

---

## 📜 Citation

If you use these scripts or workflows, please cite the article:

> *Vanbergen et al., "Genomic tracing reveals multiple independent occurrences of Bactrocera dorsalis in Belgium", [Journal Name], [Year]*

---



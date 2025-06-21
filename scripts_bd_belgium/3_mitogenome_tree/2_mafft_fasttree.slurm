#!/bin/bash -l
#SBATCH -A lp_kmma
#SBATCH -M genius
#SBATCH -p batch
#SBATCH --job-name=fastree_boot
#SBATCH -o fasttree.out
#SBATCH --error fasttree.error
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=36
#SBATCH --mem-per-cpu=3400mb


#FASTA=bd.mitogenomes.fasta
ALN=bd.mitogenomes.mafft.phy
FILTEREDALN=bd.mitogenomes.mafft.filtered.trimmed.phy
FILTEREDFASTA=bd.mitogenomes.mafft.filtered.trimmed.correctlabelformat.fasta 
# we use .fasta format because because FastTree cannot work with relaxed phylip format. correctlabelformat: ID suffixes are removed.
BOOT=bd.mitogenomes.mafft.filtered.trimmed.phy.bootstraps
BOOTTREE=bd.mitogenomes.mafft.filtered.trimmed.phy.bootstraps.tree
TREE=bd.mitogenomes.mafft.filtered.trimmed.correctlabelformat.fasta.tree
FINAL=bd.mitogenomes.mafft.filtered.trimmed.correctlabelformat.bootstrapped.tree
cd $VSC_SCRATCH/fastq/bd

# align sequences
mafft --thread -1 --phylipout --auto $FASTA > $ALN

# manually filter ALN to produce FILTEREDALN (.phy format) and FILTEREDFASTA (saved as .fasta)

# generate 500 bootstrap replicates with random seed number 1
touch outfile
printf "${FILTEREDALN}\nR\n200\nY\n1\nR\n" | seqboot
mv outfile $BOOT

# infer tree based on real sequences
FastTreeMP -nt -log $TREE.log < $FILTEREDFASTA > $TREE

# infer bootstrapped trees
FastTreeMP -nt -n 200 -log $BOOTTREE.log < $BOOT > $BOOTTREE

# add traditional boostrap values to the tree
perl -I $VSC_HOME $VSC_HOME/CompareToBootstrap.pl -tree $TREE -boot $BOOTTREE > $FINAL

# print all sample IDs of the final tree to use for metadata file
cat $FINAL | tr '(),:;' '\n' | grep -E "^[A-Za-z0-9_/-]+$" | sort -u
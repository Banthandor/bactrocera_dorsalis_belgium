# ==============================================================================
# 📌 group assignment of Bactrocera dorsalis COI haplotypes based on DAPC
# ------------------------------------------------------------------------------
# this script performs cluster assignment of COI-5P haplotypes using as input the trained DAPC model 
# AND either
# unaligned query sequences together with the reference alignment that was used 
# to train the model (if calling the alignment from R)  
# OR
# the (un)aligned query sequences and in a separate file the alignment containing
# both the training alignment and the query sequences with identical query sequence IDs for both files.  
#
#
# ⚠️ NOTE FOR WINDOWS USERS:
# If you are using Windows and do not have MAFFT installed, we recommend using
# the MAFFT online alignment service to add your sequences:
#
# 🔗 https://mafft.cbrc.jp/alignment/server/add_sequences.html
#
# Upload your existing alignment and the new sequences there. After downloading
# the updated alignment, place it in the appropriate file path below. 
# Make sure to tick the "Sequence title: Same as input" and "Keep alignment length" options.
#
# ==============================================================================


# --------- Load required libraries ---------
library(ape)
library(adegenet)
library(ips) # required only if doing the alignment with the mafft() function in R


setwd("path/to/working/directory")

# Add new sequences using MAFFT via ips (skip if the query sequences are already aligned to the reference alignment):
mafft_executeable <- read.FASTA("/user/local/bin/mafft") # path to mafft executable

alignment <- read.FASTA("bactrocera_dorsalis_coi_aln.fasta")  # existing aligned COI sequences

new <- read.FASTA("new_queries.fasta") # FASTA file with one or more query sequences

aligned_updated <- mafft(x = alignment, y = new, add = "add", method = "auto", thread = -1, exec = mafft_executable)

alignment <- aligned_updated

dapc_model <- readRDS("dapc_final_model.rds") 

# If the alignment was done beforehand (not in R), provide the alignment including the query sequences here:

alignment <- read.FASTA("updated_aln.fasta") 

new <- read.FASTA("new_queries.fasta") # FASTA file with one or more query sequences that are also in 'alignment'

dapc_model <- readRDS("dapc_final_model.rds")



# --------- Create genind object from alignment updated with query sequence  ---------

gen <- DNAbin2genind(alignment)

X <- scaleGen(gen, NA.method = "mean", scale = FALSE)


# --------- Predict group membership for each new sample ---------

new_sample_ids <- names(new)

# Extract the new individual as a matrix
new_sample_indices <- which(rownames(X) %in% new_sample_ids) # if the output of this command is not as expected, 
# make sure the query sample ids in "new" and "alignment" are matching exactly

pred <- predict(dapc_model, X[new_sample_indices, , drop = FALSE])
head(pred$assign)
posteriors <- as.data.frame(round(pred$posterior, 2))

colnames(posteriors) <- paste0("(P | cluster ",colnames(posteriors), ")")

# Add a new column with "query1", "query2", ..., "queryN"
posteriors$Query_id <- as.factor(paste0("query ", seq_len(nrow(posteriors))))

# Assuming df is a data frame where columns 1:10 are the posterior probabilities for clusters
posteriors$Cluster <- as.factor(apply(posteriors[, 1:10], 1, function(x) which.max(x))) 

# Load regional sample proportions per cluster normalized to account for 
# differences in sample size per geographic region 
cluster_proportions_normalized <- readRDS("cluster_proportions_per_region_normalized.rds")

# Restructure dataset to have regional normalized proportions as columns and clusters as rows
cluster_proportions_normalized_wide <-cluster_proportions_normalized %>%
  select(Cluster, Region, PropNormalized) %>%
  pivot_wider(names_from = Region, values_from = PropNormalized)

# join query assignments with cluster_proportions_normalized_wide
posteriors <- posteriors %>%
  left_join(cluster_proportions_normalized_wide, by = "Cluster") 

# add query sample ids as rownames
rownames(posteriors) <- new_sample_ids %>%
  na.omit()

# reorder columns for clarity
posteriors <- posteriors[, c(c("Query_id", "Cluster"),
                             setdiff(names(posteriors),
                                     c("Query_id", "Cluster")))]

# Combine query labels by cluster
query_labels <- posteriors[order(posteriors$Query_id),] %>%
  group_by(Cluster) %>%
  summarise(Label = paste(Query_id, collapse = "\n"), .groups = "drop")

# Create y position for placing labels above bars
label_positions <- cluster_proportions_normalized_wide %>%
  group_by(Cluster) %>%
  summarise(y = 1.2) %>%
  left_join(query_labels, by = "Cluster") %>%
  na.omit()

# Draw plot with Normalized regional composition of COI haplotype clusters. 
# Each bar represents one COI haplotype cluster.
# For each cluster, the proportions across regions sum to 1. 
# Regional frequencies were based on normalized counts by total sampling effort per region 
# to account for unequal sample sizes.
# The regions correspond with sampled countries in the table at the end of the script
# See the original paper for more details
ggplot(cluster_proportions_normalized, aes(x = Cluster, y = PropNormalized, fill = Region)) +
  geom_bar(stat = "identity") +
  geom_text(data = label_positions, aes(x = Cluster, y = y, label = Label), inherit.aes = FALSE, size = 3, angle = 0, hjust = 0.5) +
  theme_minimal() +
  labs(x = "Cluster", y = "Normalized proportion") +
  scale_x_discrete(labels = custom_cluster_labels) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0, 1.35)) +
  scale_fill_manual(values = continent_colors) +
  theme_classic() +
  theme(axis.line.y.right = element_blank(), axis.line.y.left = element_blank(), axis.ticks.x.bottom = element_blank() , axis.line.x.bottom = element_blank())




# --------- Result interpretation ---------

# pred contains the cluster the invidiuals were assigned to with posterior probability.
# The posterior probabillities must be taken into account
# We recommend that users familiarize themselves with DAPC by reading
# this tutorial ( especially 3.5 'Interpreting Group Memberships) : https://adegenet.r-forge.r-project.org/files/tutorial-dapc.pdf
# or any other information source on DAPC theory.



# Sampled countries per region (see paper for sample sizes per country in supplementary figures)
# ----------------------------
# Country                          | Region
# ---------------------------------|-----------
# Kenya                            | Africa
# Ethiopia                         | Africa
# Tanzania                         | Africa
# Burundi                          | Africa
# Uganda                           | Africa
# Burkina Faso                     | Africa
# Malawi                           | Africa
# South Africa                     | Africa
# Democratic Republic of the Congo | Africa
# Sudan                            | Africa
# Benin                            | Africa
# Mali                             | Africa
# Guinea                           | Africa
# Cote d'Ivoire                    | Africa
# Nigeria                          | Africa
# Senegal                          | Africa
# Mayotte                          | Africa
# Comoros                          | Africa
# Madagascar                       | Africa
# Cameroon                         | Africa
# Mozambique                       | Africa
# Ghana                            | Africa
# Togo                             | Africa
# Liberia                          | Africa
# India                            | South Asia
# Sri Lanka                        | South Asia
# Nepal                            | South Asia
# Bangladesh                       | South Asia
# Pakistan                         | South Asia
# Bhutan                           | South Asia
# China                            | East Asia
# Taiwan                           | East Asia
# Japan                            | East Asia
# Thailand                         | Southeast Asia
# Vietnam                          | Southeast Asia
# Malaysia                         | Southeast Asia
# Cambodia                         | Southeast Asia
# Myanmar                          | Southeast Asia
# Indonesia                        | Southeast Asia
# Philippines                      | Southeast Asia
# Laos                             | Southeast Asia
# Papua New Guinea                 | Pacific
# French Polynesia                 | Pacific
# United States (Hawaii/California)| Pacific
# Reunion                          | Mascarenes
# Mauritius                        | Mascarenes
# Required libraries
library(ape)
library(adegenet)
library(dplyr)
library(tidyr)
library(ggplot2)
library(treemap)
library(tibble)
library(spider)
library(stringr)

# 1. Load alignment
setwd("path/to/working/directory")

# Step 1: Load the alignment
alignment <- read.dna("bactrocera_dorsalis_coi_bold_wgs_nooutgroup_folmer.aln.trimmed.rmvfolmer.fasta", format = "fasta")

# First round of missingness filtering
miss <- checkDNA(alignment, gapsAsMissing = T)

# Total number of sites
total_sites <- ncol(alignment)

# Percentage missing
miss_pct <- 100 * miss / total_sites
plot(miss_pct,
     pch = 16,
     xlab = "Individual index",
     ylab = "% Missing")

# Set threshold, e.g., 10%
keep <- miss_pct >= 85 & miss_pct <= 99

# Filter alignment
alignment_filtered <- alignment[keep, ]

miss <- checkDNA(alignment_filtered, gapsAsMissing = T)

# Total number of sites
total_sites <- ncol(alignment_filtered)

# Percentage missing
miss_pct <- 100 * miss / total_sites

plot(miss_pct,
     pch = 16,
     col = "steelblue",
     xlab = "Individual index",
     ylab = "% Missing",
     main = "Missing Data per Individual")

# Set threshold, e.g., 10%
keep <- miss_pct <= 90

# Filter alignment
alignment_filtered <- alignment_filtered[keep, ]

# Save filtered alignment
write.dna(alignment_filtered,
          file = "bactrocera_dorsalis_coi_bold_wgs_nooutgroup_folmer.aln.trimmed.rmvfolmer.filtered.fasta",
          format = "fasta",
          nbcol = -1,      # writes full sequences on one line
          colsep = "")     # no space between bases

# Save filtered alignment

# manually removed 
#>GBMTG993-16/8-626 [Species: Bactrocera dorsalis] [Country: nan] #rRNA
#>GBART68647-23/1-522 [Species: Bactrocera dorsalis] [Country: China] #rRNA
#>GBART68648-23/1-522 [Species: Bactrocera dorsalis] [Country: China] # not COI

# Step 1: Load the alignment

alignment <- read.dna("bactrocera_dorsalis_coi_bold_wgs_nooutgroup_folmer.aln.trimmed.rmvfolmer.filtered.rmvoutliers.fasta", format = "fasta")

# First round of missingness filtering
miss <- checkDNA(alignment, gapsAsMissing = F)

# Total number of sites
total_sites <- ncol(alignment)

# Percentage missing
miss_pct <- 100 * miss / total_sites
plot(miss_pct,
     pch = 16,
     xlab = "Individual index",
     ylab = "% Missing")

# Set threshold, e.g., 10%
keep <- miss_pct <= 2

# Filter alignment
alignment_filtered <- alignment[keep, ]

alignment <- alignment_filtered



# Define country to continent mapping
country_to_continent <- function(cntry) {
  if (cntry %in% c("Kenya", "Ethiopia","Tanzania", "Burundi", "Uganda", "Burkina Faso", "Malawi", "South Africa", 
                   "Democratic Republic of the Congo", "Sudan", "Benin", "Mali", 
                   "Guinea", "Cote d'Ivoire", "Nigeria", "Senegal", "Ethiopia",
                   "Mayotte", "Comoros", "Madagascar" ,"Cameroon", "Mozambique", "Ghana", "Togo", "Liberia")) {
    return("Africa")
  } else if (cntry %in% c("India", "Sri Lanka", "Nepal", "Bangladesh", "Pakistan", "Bhutan")) {
    return("South Asia")
  } else if (cntry %in% c("China", "Taiwan", "Japan")) {
    return("East Asia")
  } else if (cntry %in% c("Thailand", "Vietnam", "Malaysia", "Cambodia", 
                          "Myanmar", "Indonesia", "Philippines", "Laos")) {
    return("Southeast Asia")
  } else if (cntry %in% c("Papua New Guinea", "French Polynesia", "United States")) {
    return("Pacific")
  } else if (cntry %in% c("Reunion", "Mauritius")) {
    return("Mascarenes")
  } else if (cntry %in% c("nan", "Unrecoverable", "NA")) {  
  } else {
    return("Other")
  }
}

# Define prefix-to-country mapping from Table S1
prefix_to_country <- c(
  HNLY = "China", JSZJ = "China", HBWH = "China", ZJZJ = "China", HNSY = "China", 
  JXNC = "China", SCLS = "China", GZDY = "China", YNBS = "China", YNJH = "China",
  GXNN = "China", FJZZ = "China", GDGZ = "China", HNDZ = "China", TWPD = "Taiwan", TWTY = "Taiwan",
  LAVT = "Laos", VNTG = "Vietnam", MMYG = "Myanmar", THBK = "Thailand",
  THPK = "Thailand", PHDM = "Philippines", PHDV = "Philippines", MYKL = "Malaysia", 
  IDJI = "Indonesia", PGPM = "Papua New Guinea",
  PKMR = "Pakistan", INPJ = "India", INAS = "India", INBH = "India", INTG = "India", INTN = "India",
  LKAD = "Sri Lanka", BDKG = "Bangladesh", NPSH = "Nepal",
  SDSG = "Sudan", SNBG = "Senegal", MLBM = "Mali", GHND = "Ghana", CIBD = "Cote d'Ivoire",
  NGNS = "Nigeria", ETMS = "Ethiopia", UGEB = "Uganda", KENB = "Kenya", CDKB = "Democratic Republic of the Congo",
  BIBJ = "Burundi", MWZB = "Malawi", ZAMP = "South Africa", MGAT = "Madagascar",
  USHW = "United States", GSP = "Reunion",
  MYRB = "B. carambolae", IDBD = "B. carambolae", SRPM = "B. carambolae",
  FAVV = "Belgium", Bd = "Belgium", GBAAY = "Italy", GBAAW = "Oman", FFIPM = "Austria"
)

# Extract headers
headers <- rownames(alignment)

# Extract country and define continent
country <- str_extract(headers, "(?<=\\[Country: ).+?(?=\\])")

# Extract prefix: characters before the first digit
prefix <- sub("^([^0-9]+).*", "\\1", headers)

# Fill in missing countries from prefix
for (i in seq_along(country)) {
  if (is.na(country[i]) && prefix[i] %in% names(prefix_to_country)) {
    country[i] <- prefix_to_country[[prefix[i]]]
  }
}

continent_raw <- sapply(country, country_to_continent)

# Get unique countries that were categorized as "Other"
countries_other <- unique(country[continent_raw == "Other"])
countries_other

# Find and print matching sequence IDs
ids_other <- headers[country %in% countries_other]
ids_other

# Rename specific sequences using mapping table
mapping <- tibble::tribble(
  ~old_id, ~new_id,
  "Bd_2024_1_MKDN240014528-1A_227WT5LT4_L8", "Bd_2024_1",
  "Bd_2024_2_MKDN240014885-1A_22C3W7LT4_L7", "Bd_2024_2",
  "Bd_2024_3_MKDN240017087-1A_22GVLGLT4_L4", "Bd_2024_3",
  "Bd_2024_4_MKDN240018839-1A_22GV3TLT4_L4", "Bd_2024_4",
  "Bd_B1_MKDN240000608-1A_22GYCNLT3_L2", "Bd_2023_7",
  "Bd_B2_MKDN240000609-1A_22GYCNLT3_L2", "Bd_2023_4",
  "Bd_B3_MKDN240000610-1A_22GYCNLT3_L2", "Bd_2023_5",
  "Bd_B4_MKDN240000611-1A_22GYCNLT3_L2", "Bd_2023_6",
  "BE249", "Bd_2023_1",
  "BE313", "Bd_2023_2",
  "BE316", "Bd_2023_3",
  "FAVV_11_MKDN240017098-1A_22GVLGLT4_L4", "FAVV_11",
  "FAVV_12_MKDN240017099-1A_22GVLGLT4_L6", "FAVV_12",
  "FAVV_13_MKDN240017100-1A_22GVLGLT4_L6", "FAVV_13",
  "FAVV_14_MKDN240017101-1A_22GVLGLT4_L7", "FAVV_14",
  "FAVV_1_MKDN240017088-1A_22GVLGLT4_L4", "FAVV_1",
  "FAVV_2_MKDN240017089-1A_22GVLGLT4_L4", "FAVV_2",
  "FAVV_3_MKDN240017090-1A_22GVLGLT4_L4", "FAVV_3",
  "FAVV_4_MKDN240017091-1A_22GVLGLT4_L4", "FAVV_4",
  "FAVV_5_MKDN240017092-1A_22GVLGLT4_L4", "FAVV_5",
  "FAVV_6_MKDN240017093-1A_22C3W5LT4_L7", "FAVV_6",
  "FAVV_7_MKDN240017094-1A_22GVLGLT4_L4", "FAVV_7",
  "GBAAY34563-24", "GBAAY34563-24/Italy",
  "GBAAY34565-24", "GBAAY34565-24/Italy",
  "GBAAY34588-24", "GBAAY34588-24/Italy",
  "GBAAY34813-24", "GBAAY34813-24/Italy",
  "FFIPM495-22",   "FFIPM495-22/Austria",
  "GBMNC2473-20", "GBMNC2473-20/Australia",
  "GBAAW3233-24", "GBAAW3233-24/Oman"
)

# Create named replacement vector from mapping
replacement_vec <- setNames(mapping$new_id, mapping$old_id)

# Get current rownames (these are your alignment labels)
alignment_ids <- rownames(alignment)

# Extract prefix from alignment IDs to match mapping keys
alignment_prefix <- stringr::str_extract(alignment_ids, "^[^/\\s\\[]+")

# Perform the replacement based on prefix
alignment_ids_clean <- ifelse(
  alignment_prefix %in% names(replacement_vec),
  replacement_vec[alignment_prefix],
  alignment_ids
)

# Assign updated names back to the alignment
rownames(alignment) <- alignment_ids_clean

write.dna(alignment,
          file = "bactrocera_dorsalis_coi_bold_wgs_nooutgroup_folmer.aln.trimmed.rmvfolmer.filtered.max2percentmiss.fasta",
          format = "fasta",
          nbcol = -1,      # writes full sequences on one line
          colsep = "")     # no space between bases

# Get current and cleaned names again
original_ids <- alignment_ids
new_ids <- rownames(alignment)

# Compare and print changes
changed <- original_ids != new_ids
data.frame(
  original = original_ids[changed],
  updated  = new_ids[changed]
)

# Convert to genind
gen <- DNAbin2genind(alignment)

set.seed(123)
diffNgroup_grp <- find.clusters(gen, n.pca = 80, n.clust = NULL,
                                , max.n.clust = 40,
                                n.iter = 1e6, n.start = 500, scale = FALSE,
                                choose.n.clust = F, criterion ="diffNgroup")

plot(diffNgroup_grp$Kstat, type="o", xlab="number of clusters (K)", ylab="BIC",
     col="blue")
points(10, diffNgroup_grp$Kstat[10], pch="O", cex=2)

k.prelim <- nlevels(diffNgroup_grp$grp)
k.prelim

# Impute and scale
X <- scaleGen(gen, NA.method = "mean", scale = FALSE)

# Cross-validate DAPC
set.seed(123)
xval <- xvalDapc(
  X,               # genetic data matrix
  grp = diffNgroup_grp$grp,                  # cluster assignments
  n.pca = c(10, 15, 20, 25, 30, 35, 40,45 , 50, 55, 60, 70, 80, 90, 100), 
  n.pca.max = 100,                # max PCs to test
  training.set = 0.9,             # 90% training data
  result = "groupMean",           # method for prediction
  center = TRUE, scale = FALSE,
  n.rep = 30                      # number of replicates
)

# Check optimal number of PCs
xval$`Number of PCs Achieving Highest Mean Success`

# Construct final dapc model
set.seed(123)
dapc_final <- dapc(X, grp = diffNgroup_grp$grp, n.pca = as.numeric(xval$`Number of PCs Achieving Highest Mean Success`), n.da = k.prelim - 1 , return.posterior = TRUE)

saveRDS(dapc_final, file = "dapc_final_model.rds")

dapc_final <- readRDS(file = "dapc_final_model.rds")
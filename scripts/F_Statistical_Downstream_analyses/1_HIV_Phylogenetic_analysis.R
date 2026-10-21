################################################################################
##### 2000HIV: Phylogenetic analysis
### Author(s):Asier Fernández-Pato
### Last updated: 15th September, 2026
################################################################################

#****************
# Load modules
#****************
library(dplyr)
library(readr)
library(stringr)
library(ape)
library(ggplot2)
library(ggtree)
library(ggtreeExtra)
library(tidyr)
library(ggnewscale)
library(picante)
library(phytools)

#****************
# Define functions
#****************
# Function to read the domtblout files
read_domtbl <- function(file) {
  pfam <- str_split(basename(file), "\\.")[[1]][1]
  
  read_table(file, comment = "#", col_names = FALSE) %>%
    transmute(
      protein = X1,
      pfam = pfam,
      score = X8,
      hmm_from = X16,
      hmm_to = X17,
      ali_from = X18,
      ali_to = X19,
      hmm_len = X6
    ) %>%
    mutate(coverage = (hmm_to - hmm_from + 1) / hmm_len)
}

# Set working directory 
setwd("~/Desktop/PhD/Projects/HIV 2000/")

# Read abundance table and metadata tables
Abundance_table <- read.delim("ABUNDANCE TABLE/HIV_Abundance_Table_01052025.txt")
Sample_metadata <- read.delim("METADATA/Sample_metadata_11052026.txt")
Viral_metadata <- read.delim("METADATA/HIV_Viral_Metadata_06052025.txt")


##************************************************************************
# 1. Processing HMM matches to Caudoviricetes markers
#*************************************************************************

# We generated HMMs from different PFAM entries covering different Caudoviricetes markers
# They included: TerL ATPase domain, TerL nuclease domain,  Major capsid protein and Portal protein

# Load all the hmmsearch results
files <- list.files("7_PHYLOGENY/HMMsearch_results/", pattern="*.domtblout", full.names=TRUE)
results <- bind_rows(lapply(files, read_domtbl))
results$genome <- sub("_[^_]+$", "", results$protein)

# Filter results to select only significant hits according to out thresholds
# Terl: 10 and others 25. Coverage >= 0.6
terl <- c("PF03237","PF03354","PF17289","PF20441","PF20454")

results <- results %>%
  mutate(threshold = ifelse(pfam %in% terl, 10, 25)) %>%
  filter(score >= threshold, coverage >= 0.6)

# Filter results to select only best hit per protein
best_hits <- results %>%
  group_by(protein) %>%
  slice_max(score, n = 1, with_ties = FALSE) %>%
  ungroup()

# Also, in case of the same genome having multiple proteins matching the same PFAM, select one (highest score) 
best_hits <- best_hits %>%
  group_by(genome, pfam) %>%
  slice_max(score, n = 1, with_ties = FALSE) %>%
  ungroup()

write_tsv(best_hits, "7_PHYLOGENY/HMM_best_hits.tsv")

##************************************************************************
# 2. Visualizing Caudoviricetes phylogenetic tree
#*************************************************************************

# Load the tree
tree <- read.tree("7_PHYLOGENY/tree.nwk") 

################################
# Prepare metadata for the plot (20,130 phages included)
################################

metadata_plot <- Viral_metadata %>%
  dplyr::rename(label = Virus_ID) %>%
  dplyr::select(label, Length, Quality, vOTU_DB_composition, Lifestyle,
                Bacterial_family_host)

# Quality binary
metadata_plot$Quality_binary <- ifelse(
  metadata_plot$Quality %in% c("Complete", "High-quality"),
  "High-quality", "No"
)

# Top 10 families
top10_families <- metadata_plot %>%
  dplyr::count(Bacterial_family_host, sort = TRUE) %>%
  slice_head(n = 10) %>%
  pull(Bacterial_family_host)

# Update Bacterial_family_host
metadata_plot <- metadata_plot %>%
  mutate(
    Bacterial_family_host = ifelse(
      Bacterial_family_host %in% top10_families,
      Bacterial_family_host,
      "Other"
    ),
    Bacterial_family_host = ifelse(
      is.na(Bacterial_family_host),
      "Unassigned",
      Bacterial_family_host
    ),
    Bacterial_family_host = factor(
      Bacterial_family_host,
      levels = c(top10_families, "Other", "Unassigned")
    )
  )

# Ensure only tips present in tree
metadata_plot <- metadata_plot %>% filter(label %in% tree$tip.label)

# Lifestyle: replace NA with Unassigned and factor
metadata_plot$Lifestyle[is.na(metadata_plot$Lifestyle)] <- "Unassigned"
metadata_plot$Lifestyle <- factor(
  metadata_plot$Lifestyle,
  levels = c("Virulent", "Temperate", "Unassigned")
)

# Quality binary factor
metadata_plot$Quality_binary <- factor(
  metadata_plot$Quality_binary,
  levels = c("High-quality", "No")
)

################################
# Define color palettes
################################

# Top families
family_levels <- levels(metadata_plot$Bacterial_family_host)
family_colors <- setNames(
  c("#6FCFC0", "#6B9FE6", "#E178B3", "#E5986A",
    "#7FB8EB", "#D4706A", "#E6C276", "#A18FE6",
    "#7FD18C", "grey90", "white"), family_levels)

# Lifestyle
lifestyle_colors <- c("Virulent" = "#FFD54F",  
  "Temperate" = "#5C6BC0", "Unassigned" = "white")

quality_colors <- c("High-quality" = "#546E7A",
  "No" = "white")

novelty_colors <- c("HIV" = "#D16C94", "OTHER_DBs" = "grey90")

################################
# Generate tree layers
################################

# Base circular tree
p <- ggtree(tree, layout = "circular", size = 0.2, color = "grey30") +
  theme(
    legend.position = "right",
    plot.background = element_blank()
  )

# Quality
p1 <- p +
  new_scale_fill() +
  geom_fruit(
    data = metadata_plot,
    geom = geom_tile,
    mapping = aes(y = label, fill = Quality_binary),
    width = 0.7,
    offset = 0.065,
    alpha = 0.8
  ) +
  scale_fill_manual(
    values = quality_colors,
    guide = guide_legend(override.aes = list(alpha = 0.8))
  )

# Novelty
p2 <- p1 +
  new_scale_fill() +
  geom_fruit(
    data = metadata_plot,
    geom = geom_tile,
    mapping = aes(y = label, fill = vOTU_DB_composition),
    width = 0.7,
    offset = 0.065,
    alpha = 0.8
  ) +
  scale_fill_manual(
    values = novelty_colors,
    guide = guide_legend(override.aes = list(alpha = 0.8))
  )

# Lifestyle
p3 <- p2 +
  new_scale_fill() +
  geom_fruit(
    data = metadata_plot,
    geom = geom_tile,
    mapping = aes(y = label, fill = Lifestyle),
    width = 0.7,
    offset = 0.065,
    alpha = 0.8
  ) +
  scale_fill_manual(
    values = lifestyle_colors,
    guide = guide_legend(override.aes = list(alpha = 0.8))
  )

# Host family
p4 <- p3 +
  new_scale_fill() + 
  geom_fruit(
    data = metadata_plot,
    geom = geom_tile,
    mapping = aes(y = label, fill = Bacterial_family_host),
    width = 0.7,
    offset = 0.065,
    alpha = 0.8
  ) +
  scale_fill_manual(
    values = family_colors,
    guide = guide_legend(override.aes = list(alpha = 0.8))
  )

# Genome length
p_final <- p4 +
  geom_fruit(
    data = metadata_plot,
    geom = geom_col,
    mapping = aes(y = label, x = Length),
    orientation = "y",
    pwidth = 0.4,
    offset = 0.045,
    color = "grey40",
    alpha = 0.8
  ) +
  geom_treescale(x = 10, y = 50, width = 2.5)

# Save the tree
ggsave("7_PHYLOGENY/Circular_tree_plot.pdf", plot = p_final, width = 10, height = 10, units = "in")

##************************************************************************
# 3. Estimating Faith's phylogenetic diversity for each sample
#*************************************************************************

# Transpose abundance table
Abundance_table_t <- t(Abundance_table)

# Ensure only vOTUs present in tree are present
Abundance_table_tree <- data.frame(Abundance_table_t[, colnames(Abundance_table_t) %in% tree$tip.label])

# Estimate Compute Faith's PD
faith_pd <- picante::pd(
  Abundance_table_tree,
  tree,
  include.root = FALSE
)

faith_pd$SR <- NULL

# Add to sample metadata
Sample_metadata <- Sample_metadata %>%
  mutate(SampleID = rownames(faith_pd)) %>%
  left_join(
    data.frame(
      SampleID = rownames(faith_pd),
      Faith_PD = faith_pd$PD
    ),
    by = "SampleID"
  )

##************************************************************************
# 4. Estimating association of distances with metadata
#*************************************************************************

# As these are not continuous variables, Pagel's λ cannot be used
# We decide to use a permutation-based approach

##################
# A. Lifestyle
##################

lifestyle <- setNames(Viral_metadata$Lifestyle, Viral_metadata$Virus_ID)
lifestyle <- lifestyle[tree$tip.label]
lifestyle <- lifestyle[!is.na(lifestyle)]

tree_lifestyle <- drop.tip(tree, setdiff(tree$tip.label, names(lifestyle)))
lifestyle <- lifestyle[tree_lifestyle$tip.label]

# Phylogenetic distances
phylo_dist <- cophenetic(tree_lifestyle)

# Function to calculate mean within-group phylogenetic distance
calc_stat <- function(x) {
  mean(unlist(lapply(split(names(x), x), function(tips) {
    if (length(tips) < 2) return(NULL)
    d <- phylo_dist[tips, tips]
    d[lower.tri(d)]
  })))
}

# Observed test statistic
observed_stat <- calc_stat(lifestyle)

# Permutation test
num_permutations <- 1000
permuted_stats <- numeric(num_permutations)

set.seed(123)

for (i in 1:num_permutations) {
  permuted_lifestyle <- lifestyle
  permuted_lifestyle[] <- sample(lifestyle)
  permuted_stats[i] <- calc_stat(permuted_lifestyle)
}

# One-sided p-value
# smaller distances = stronger phylogenetic clustering
p_value <- (sum(permuted_stats <= observed_stat) + 1) /
  (num_permutations + 1)

observed_stat
p_value

##################
# B. Bacterial host
##################

host_family <- setNames(Viral_metadata$Bacterial_family_host, Viral_metadata$Virus_ID)
host_family <- host_family[tree$tip.label]
host_family <- host_family[!is.na(host_family)]

tree_host <- drop.tip(tree, setdiff(tree$tip.label, names(host_family)))
host_family <- host_family[tree_host$tip.label]

# Phylogenetic distances
phylo_dist_host <- cophenetic(tree_host)

# Function to calculate mean within-group phylogenetic distance
calc_stat_host <- function(x) {
  mean(unlist(lapply(split(names(x), x), function(tips) {
    if (length(tips) < 2) return(NULL)
    d <- phylo_dist_host[tips, tips]
    d[lower.tri(d)]
  })))
}

# Observed test statistic
observed_stat_host <- calc_stat_host(host_family)

# Permutation test
num_permutations <- 1000
permuted_stats_host <- numeric(num_permutations)

set.seed(123)

for (i in 1:num_permutations) {
  permuted_host <- host_family
  permuted_host[] <- sample(host_family)
  permuted_stats_host[i] <- calc_stat_host(permuted_host)
}

# One-sided p-value
p_value_host <- (sum(permuted_stats_host <= observed_stat_host) + 1) /
  (num_permutations + 1)

observed_stat_host
p_value_host


##************************************************************************
# Write results
##************************************************************************

# Select only phenotypes of interest
phenotypes_keep <- c(
  "SampleID", "CENTER", "COHORT", "ETHNICITY", "SEX_BIRTH", "GEN",
  "AGE", "DATE_VISIT", "COVID_GROUP", "BMI_BASELINE", "HIV_DURATION",
  "CD4_NADIR", "CD4CD8_PRECART", "CART_DURATION", "CD4_LATEST", "CD8_LATEST",
  "Steatosis_yesno", "Fibrosis_yesno", "Total_million_avg", "MSM", "MSM_derived",
  "log10_counts", "EC_yesno", "EC_lost", "EC_vir", "EC_combined",
  "EC_ART", "INR_YESNO", "Faith_PD"
)

Sample_metadata <- Sample_metadata[, phenotypes_keep]

write.table(Sample_metadata, file = "METADATA/Sample_metadata_11052026.txt",
            sep = "\t", row.names = T, quote = FALSE)




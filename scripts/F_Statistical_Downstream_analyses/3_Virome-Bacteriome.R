################################################################################
##### 2000HIV: Virome-Bacteriome analysis
### Author(s):Asier Fernández-Pato
### Last updated: 16th September, 2026
################################################################################

#****************
# Load modules
#****************
library(dplyr)
library(stringr)
library(vegan)
library(purrr)
library(ggplot2)


#****************
# Define functions
#****************

# Function to process Metaphlan4 abundance table
Clean_MP4 = function(Microbial){
  unknown_rel = Microbial$UNCLASSIFIED
  Microbial %>% select(-UNCLASSIFIED) -> M2
  lapply(colnames(M2)[2:dim(M2)[2]], FUN=function(x){ str_split(x,"\\|")[[1]] -> y ; trimws(y[length(y)]) }) %>% as_vector() -> New_names
  colnames(M2) = c("ID", New_names)
  M2 = M2[, !grepl("(^$|g__$|s__$)", colnames(M2)), drop=FALSE] # remove unassigned taxa
  colnames(M2) = make.unique(colnames(M2), sep="_dup")           # avoid duplicate names
  lapply(str_split(colnames(select(M2, -ID)), "__"), FUN=function(x){ x[1] }) %>% as_vector() -> Taxonomies
  Taxonomies %>% unique() -> Taxonomy_Levels
  Bugs = colnames(select(M2, -ID)) ; tibble(Bug = Bugs, Level = Taxonomies) -> Bug_array
  M3 = tibble(ID = M2$ID)
  for (T_L in Taxonomy_Levels){
    filter(Bug_array, Level == T_L) -> Bugs
    Back_to_100(select(M2, Bugs$Bug)) -> Normalized
    M3 = left_join(M3, mutate(Normalized, ID = M2$ID), by="ID")
  }
  return(list(unknown_rel, M2, M3, Bug_array))
}


# Function to recalculate the relative abundance after removing the unclassified
Back_to_100 = function(Microbial){
  #We can go back to 100%, but this function should be applied per taxonomic level
  apply(Microbial,1, FUN= function(x){ sum(x) -> y ;  return(y) }) %>% as_vector() -> Total_per_sample
  Microbial/Total_per_sample -> Normalized
  return(as_tibble(Normalized))
}


# Set working directory 
setwd("~/Desktop/PhD/Projects/HIV 2000/")

# Read abundance table and metadata tables
Abundance_table <- read.delim("ABUNDANCE TABLE/HIV_Abundance_Table_01052025.txt")
Abundance_table_host <- read.delim("5_HOST_PREDICTION/phage_abundance_by_host_genus_relative_abundance.tsv")
Abundance_table_bact <- read.delim("ABUNDANCE TABLE/results_metaphlan4.txt")
Sample_metadata <- read.delim("METADATA/Sample_metadata_excluded_11052026.txt")
Viral_metadata <- read.delim("METADATA/HIV_Viral_Metadata_06052025.txt")
Name_matching <-  readRDS("METADATA/name_trans.rds")


##************************************************************************
# 0. Preprocessing of the abundance tables
#*************************************************************************

# Make sure all tables are the same format
rownames(Sample_metadata) <- Sample_metadata$SampleID
rownames(Abundance_table_host) <- Abundance_table_host$Host_genus
Abundance_table_host$Host_genus <- NULL

# Remove excluded participants
Abundance_table <- Abundance_table [, colnames(Abundance_table) %in% Sample_metadata$SampleID]
Abundance_table_host <- Abundance_table_host [, colnames(Abundance_table_host) %in% Sample_metadata$SampleID]

# Get genus-level abundances of bacteria

# Adapt colnames and transpose
colnames(Abundance_table_bact) <- sub("_[^_]*$", "", colnames(Abundance_table_bact))
rownames(Abundance_table_bact) <- Abundance_table_bact[, 1]
Abundance_table_bact_t <- t(Abundance_table_bact[, -1])
Abundance_table_bact <- as.data.frame(Abundance_table_bact_t)

# Change rownames ID
rownames(Abundance_table_bact) <- Name_matching$sample_name[
  match(rownames(Abundance_table_bact), Name_matching$sample_name2)
]

# Select only those samples present in Sample metadata
bacterial_abundances <- Abundance_table_bact[rownames(Abundance_table_bact) %in% Sample_metadata$SampleID,]

# Process the table to use it as input of Clean_MP4 function
tax <- bacterial_abundances
tax$ID <- rownames(bacterial_abundances)
tax$ID=make.unique(tax$ID, sep = '_')
tax <- tax[,c(ncol(tax), 1:(ncol(tax) - 1))] # bringing ID to the 1st column

# Process Metaphlan4 table
List_microbial = Clean_MP4(tax) # This removes the unclassified column and shortens the taxa name. 
Unknown_fraction = List_microbial[[1]] # unknown fractions
Microbial_abundance = List_microbial[[2]] # relative abundances without re-assessment
Microbial_abundance_reassigned = List_microbial[[3]] # relative abundances after re-assessment
Taxonomy_level = List_microbial[[4]] # Microbes and their specific taxonomic levels

# Get species level abundance table
microbes <-as.data.frame(Microbial_abundance_reassigned)
row.names(microbes)<-microbes$ID
microbes$ID=NULL
species_abundance <- microbes[,grep("s__",colnames(microbes))]
species_abundance <- species_abundance[,-which(colSums(species_abundance)==0)] #remove those species not present
colnames(species_abundance) <- sub(".*s__", "", colnames(species_abundance))
species_abundance <- species_abundance[, colnames(species_abundance) != ""]

# Get genus level abundance table
genus_abundance <- microbes[, grepl("g__[^;]+$", colnames(microbes)), drop = FALSE]
colnames(genus_abundance) <- sub(".*g__", "", colnames(genus_abundance))
genus_abundance <- genus_abundance[, colSums(genus_abundance) != 0, drop = FALSE]

# Reorder according to Sample metadata
genus_abundance <- genus_abundance[match(Sample_metadata$SampleID, rownames(genus_abundance)), , drop = FALSE]
species_abundance <- species_abundance[match(Sample_metadata$SampleID, rownames(species_abundance)), , drop = FALSE]

##************************************************************************
# 1. Procrustes analysis
#*************************************************************************

# Reorder abundance tables
viral_ab <- as.data.frame(t(Abundance_table))                  
viral_ab_host <- as.data.frame(t(Abundance_table_host))       
bacteria_species_ab <- species_abundance                       
bacteria_genus_ab <- genus_abundance  

############################
# 1. Overall bacteriome–virome concordance
############################

# Prefilter taxa present in >=5% of samples
bacteria_species_f <- bacteria_species_ab[, colMeans(bacteria_species_ab > 0) >= 0.05, drop = FALSE]
viral_f <- viral_ab[ , colMeans(viral_ab > 0) >= 0.05, drop = FALSE]

# Remove samples with no viral signal
viral_f <- viral_f[rowSums(viral_f, na.rm = TRUE) > 0, , drop = FALSE]

# Convert viral RPKM to relative abundance within each sample
viral_f <- viral_f / rowSums(viral_f, na.rm = TRUE)

# Keep only samples present in both datasets
common_samples <- intersect(rownames(viral_f),
  rownames(bacteria_species_f))

viral_ab_general <- viral_f[common_samples, , drop = FALSE]
bacteria_species_general <- bacteria_species_f[rownames(viral_ab_general), , drop = FALSE]

# Bray-Curtis distances
viral_dist <- vegdist(viral_ab_general, method = "bray")
bact_species_dist <- vegdist(bacteria_species_general, method = "bray")

# NMDS
viral_nmds <- metaMDS(viral_ab_general, distance = "bray", k = 2)
bact_species_nmds <- metaMDS(bacteria_species_general, distance = "bray", k = 2)

# Procrustes
proc_general <- protest(
  scores(bact_species_nmds, "sites"),
  scores(viral_nmds, "sites"),
  permutations = 10000
)

# Extract aligned coordinates
bact_coords <- proc_general$X
viral_coords <- proc_general$Yrot

# Plot dataframe
plot_nmds <- rbind(
  data.frame(NMDS1 = bact_coords[, 1], NMDS2 = bact_coords[, 2],
    Community = "Bacteriome", Sample_ID = rownames(bact_coords)),
  data.frame(NMDS1 = viral_coords[, 1], NMDS2 = viral_coords[, 2],
    Community = "Virome", Sample_ID = rownames(viral_coords))
)

# Remove outliers within each community
plot_nmds <- plot_nmds %>%
  group_by(Community) %>%
  mutate(
    centroid_NMDS1 = mean(NMDS1),
    centroid_NMDS2 = mean(NMDS2),
    distance_from_centroid = sqrt(
      (NMDS1 - centroid_NMDS1)^2 +
        (NMDS2 - centroid_NMDS2)^2
    ),
    Q1 = quantile(distance_from_centroid, 0.25),
    Q3 = quantile(distance_from_centroid, 0.75),
    IQR = Q3 - Q1,
    upper_bound = Q3 + 1.5 * IQR
  ) %>%
  filter(distance_from_centroid <= upper_bound) %>%
  ungroup()

# Plot
pdf("Plots/Virome_bacteriome_overall_concordance.pdf", width = 3.2, height = 3.2)
ggplot(plot_nmds, aes(NMDS1, NMDS2, colour = Community)) +
  geom_point(shape = 1, size = 2.2, stroke = 0.75, alpha = 0.45) +
  scale_colour_manual(values = c("Virome" = "#3F4A8A", "Bacteriome" = "#9B4F72")) +
  labs(x = "NMDS1", y = "NMDS2", colour = NULL) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.75),
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 13),
    axis.title.y = element_text(margin = margin(r = 8)),
    axis.title.x = element_text(margin = margin(t = 8)),
    legend.position = c(0.22, 0.91),
    legend.text = element_text(size = 10),
    legend.key = element_blank()
  )
dev.off()


############################
# 2. Host-level concordance
############################

# Keep only genera represented in both 
common_genera <- intersect(colnames(bacteria_genus_ab), colnames(viral_ab_host))
bacteria_genus_host <- bacteria_genus_ab[, common_genera, drop = FALSE]
viral_host_common <- viral_ab_host[, common_genera, drop = FALSE]

# Prefilter genera present in >=5% of samples in both 
keep_genus <- colMeans(bacteria_genus_host > 0) >= 0.05 &
  colMeans(viral_host_common > 0) >= 0.05

bacteria_genus_host <- bacteria_genus_host[, keep_genus, drop = FALSE]
viral_host_common <- viral_host_common[, keep_genus, drop = FALSE]

# Keep samples with non-zero abundance in both communities
keep <- rowSums(viral_host_common) > 0 &
  rowSums(bacteria_genus_host) > 0

viral_host_common <- viral_host_common[keep, , drop = FALSE]
bacteria_genus_host <- bacteria_genus_host[rownames(viral_host_common), , drop = FALSE]

# Bray-Curtis distances
viral_host_dist <- vegdist(viral_host_common, method = "bray")
bact_genus_dist <- vegdist(bacteria_genus_host, method = "bray")

# NMDS
viral_host_nmds <- metaMDS(viral_host_common, distance = "bray", k = 2)
bact_genus_nmds <- metaMDS(bacteria_genus_host, distance = "bray", k = 2)

# Procrustes
proc_host <- protest(
  scores(bact_genus_nmds, "sites"),
  scores(viral_host_nmds, "sites"),
  permutations = 10000
)

# Plot NMDS results

# Extract aligned coordinates
bact_coords <- proc_host$X
viral_coords <- proc_host$Yrot

# Plot dataframe
plot_nmds <- rbind(
  data.frame(NMDS1 = bact_coords[, 1], NMDS2 = bact_coords[, 2],
    Community = "Bacteriome",Sample_ID = rownames(bact_coords)),
  data.frame(NMDS1 = viral_coords[, 1], NMDS2 = viral_coords[, 2],
    Community = "Virome", Sample_ID = rownames(viral_coords))
)

# Remove outliers within each community
plot_nmds <- plot_nmds %>%
  group_by(Community) %>%
  mutate(
    centroid_NMDS1 = mean(NMDS1),
    centroid_NMDS2 = mean(NMDS2),
    distance_from_centroid = sqrt(
      (NMDS1 - centroid_NMDS1)^2 +
        (NMDS2 - centroid_NMDS2)^2
    ),
    Q1 = quantile(distance_from_centroid, 0.25),
    Q3 = quantile(distance_from_centroid, 0.75),
    IQR = Q3 - Q1,
    upper_bound = Q3 + 1.5 * IQR
  ) %>%
  filter(distance_from_centroid <= upper_bound) %>%
  ungroup()


pdf("Plots/Virome_bacteriome_host_concordance.pdf", width = 3.2, height = 3.2)
ggplot(plot_nmds, aes(NMDS1, NMDS2, colour = Community)) +
  geom_point(shape = 1, size = 2.2, stroke = 0.75, alpha = 0.45) +
  scale_colour_manual(
    values = c("Virome" = "#3F4A8A", "Bacteriome" = "#9B4F72")
  ) +
  labs(x = "NMDS1", y = "NMDS2", colour = NULL) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.75),
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 13),
    axis.title.y = element_text(margin = margin(r = 8)),
    axis.title.x = element_text(margin = margin(t = 8)),
    legend.position = c(0.22, 0.91),
    legend.text = element_text(size = 10),
    legend.key = element_blank()
  )
dev.off()

#*************************************************************************
# 2. Host-level correlations
#*************************************************************************

# Spearman correlations between viral and bacterial host groups
cor_results <- expand.grid(
  Viral_host = colnames(viral_host_common),
  Bacterial_host = colnames(bacteria_genus_host),
  stringsAsFactors = FALSE
)

cor_results$rho <- NA
cor_results$p <- NA

for(i in seq_len(nrow(cor_results))) {
  
  test <- cor.test(
    viral_host_common[, cor_results$Viral_host[i]],
    bacteria_genus_host[, cor_results$Bacterial_host[i]],
    method = "spearman",
    exact = FALSE
  )
  
  cor_results$rho[i] <- unname(test$estimate)
  cor_results$p[i] <- test$p.value
}

# FDR correction
cor_results$FDR <- p.adjust(cor_results$p, method = "BH")

# Same-host correlations
same_host <- subset(
  cor_results,
  Viral_host == Bacterial_host
)

# Sort by correlation
same_host <- same_host[order(-same_host$rho), ]

# Significant same-host associations
same_host_sig <- subset(same_host, FDR < 0.05)

# Summary
nrow(same_host_sig)
sum(same_host_sig$rho > 0)
sum(same_host_sig$rho < 0)

###################
# Generate a plot 
###################

plot_df <- same_host %>%
  arrange(rho) %>%
  mutate(
    Host = factor(Viral_host, levels = Viral_host),
    Column = rep(c("A", "B"), each = ceiling(n()/2))
  )

pdf("Plots/Host_level_correlation_coefficients.pdf", width = 5.5, height = 5.5)
ggplot(plot_df, aes(Host, rho, fill = rho > 0)) +
  geom_col(width = 0.8) +
  geom_hline(yintercept = 0, colour = "black", linewidth = 0.5) +
  coord_flip() +
  facet_wrap(~Column, ncol = 2, scales = "free_y") +
  scale_fill_manual(values = c("TRUE" = "#3F8F98", "FALSE" = "#A84F45"), guide = "none") +
  labs(x = NULL, y = "Spearman ρ") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank(), strip.text = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.75),
    axis.text = element_text(size = 10), axis.title = element_text(size = 13)
  )
dev.off()


##************************************************************************
# Write results
##************************************************************************

# Save correlation results
write.table(same_host,
  file = "BACTERIA_VIRUS_CORRELATION/Same_host_all_correlations.txt",
  sep = "\t", row.names = FALSE, quote = FALSE)

# Save species and genus-level bacterial abundances
species_abundance <- as.data.frame(t(species_abundance))
genus_abundance <- as.data.frame(t(genus_abundance))
species_abundance <- species_abundance[, Sample_metadata$SampleID]
genus_abundance <- genus_abundance[, Sample_metadata$SampleID]

write.table(species_abundance,
  file = "ABUNDANCE TABLE/Bacterial_species_abundance_11052026.txt",
  sep = "\t", row.names = TRUE, quote = FALSE)

write.table(genus_abundance,
  file = "ABUNDANCE TABLE/Bacterial_genus_abundance_11052026.txt",
  sep = "\t", row.names = TRUE, quote = FALSE)




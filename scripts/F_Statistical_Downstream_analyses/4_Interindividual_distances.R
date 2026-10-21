################################################################################
##### 2000HIV: Interindividual viral distances
### Author(s):Asier Fernández-Pato
### Last updated: 16th September, 2026
################################################################################

#****************
# Load modules
#****************
library(dplyr)
library(stringr)
library(vegan)

# Set working directory 
setwd("~/Desktop/PhD/Projects/HIV 2000/")

# Read abundance table and metadata tables
Abundance_table <- read.delim("ABUNDANCE TABLE/HIV_Abundance_Table_01052025.txt")
Sample_metadata <- read.delim("METADATA/Sample_metadata_excluded_11052026.txt")
Viral_metadata <- read.delim("METADATA/HIV_Viral_Metadata_06052025.txt")

# Remove excluded participants
Abundance_table <- Abundance_table [, colnames(Abundance_table) %in% Sample_metadata$SampleID]

# Filter vOTUs with <1% prevalence
prevalence <- rowMeans(Abundance_table > 0, na.rm = TRUE)
Abundance_table <- Abundance_table[prevalence >= 0.01, , drop = FALSE]

##************************************************************************
# 1. Between-individual distance comparison according to sexual behaviour
#*************************************************************************

# Comparison of between-individual distances according to phenotypes
# First, calculate dissimilarities in virome composition between infant samples
Ab_t <- t(Abundance_table)
Ab_t_nonzero <- Ab_t[rowSums(Ab_t) > 0, ] # remove 2 samples with no viruses detected
dist_vOTUs <- vegdist(Ab_t_nonzero, method = "bray")

# Transform distance matrix into a long format
dist_mat <- as.matrix(dist_vOTUs)
dist_mat[lower.tri(dist_mat, diag = TRUE)] <- NA  
distances <- reshape2::melt(dist_mat, na.rm = TRUE)
colnames(distances) <- c("Sample1", "Sample2", "Distance")

# Add metadata to the long-formatted table
distances$Center1 <- Sample_metadata$CENTER[match(distances$Sample1, rownames(Sample_metadata))]
distances$Center2 <-  Sample_metadata$CENTER[match(distances$Sample2, rownames(Sample_metadata))]
distances$Sex1 <- recode(Sample_metadata$SEX_BIRTH[match(distances$Sample1, rownames(Sample_metadata))],
  `0` = "Man", `1` = "Woman") 
distances$Sex2 <- recode(Sample_metadata$SEX_BIRTH[match(distances$Sample2, rownames(Sample_metadata))],
                         `0` = "Man", `1` = "Woman") 
distances$Sex_behaviour1 <-  Sample_metadata$MSM_derived[match(distances$Sample1, rownames(Sample_metadata))]
distances$Sex_behaviour2 <-  Sample_metadata$MSM_derived[match(distances$Sample2, rownames(Sample_metadata))]
distances$EC1 <-  Sample_metadata$EC_yesno[match(distances$Sample1, rownames(Sample_metadata))]
distances$EC2 <-  Sample_metadata$EC_yesno[match(distances$Sample2, rownames(Sample_metadata))]
distances$INR1 <-  Sample_metadata$INR_YESNO[match(distances$Sample1, rownames(Sample_metadata))]
distances$INR2 <-  Sample_metadata$INR_YESNO[match(distances$Sample2, rownames(Sample_metadata))]

# Generate a Behaviour_pair variable
distances <- distances %>%
  filter(!is.na(Sex_behaviour1) & !is.na(Sex_behaviour2)) %>%
  mutate(
    Behaviour_pair = recode(
      paste0(pmin(Sex_behaviour1, Sex_behaviour2), "-",
             pmax(Sex_behaviour1, Sex_behaviour2)),
      "msm-msm" = "MSM-MSM",
      "msw-msw" = "MSW-MSW",
      "wsm-wsm" = "WSM-WSM",
      "msm-msw" = "MSM-MSW",
      "msm-wsm" = "MSM-WSM",
      "msw-wsm" = "MSW-WSM"
    )
  )

# Generate a metadata only for men 
distances_men <- distances %>%
  filter(Sex_behaviour1 != "wsm" & Sex_behaviour2 != "wsm")


######################
# Statistical tests
######################

set.seed(123)

######################
# A) MSW higher similarity to women than MSM
######################

# Observed difference
obs_means_A <- tapply(distances$Distance, distances$Behaviour_pair, mean)
obs_diff_A <- obs_means_A["MSM-WSM"] - obs_means_A["MSW-WSM"]

# Permutation
n_perm <- 1000
perm_diff_A <- numeric(n_perm)

for(i in 1:n_perm){
  
  # Shuffle behaviour labels across all individuals
  shuffled_labels <- sample(Sample_metadata$MSM_derived)
  names(shuffled_labels) <- Sample_metadata$SampleID
  
  # Re-assign shuffled labels
  Shuf1 <- shuffled_labels[distances$Sample1]
  Shuf2 <- shuffled_labels[distances$Sample2]
  
  # Define pair types
  pair_shuf <- paste0(pmin(Shuf1, Shuf2), "-", pmax(Shuf1, Shuf2))
  
  # MSM-WSM minus MSW-WSM
  perm_diff_A[i] <- mean(distances$Distance[pair_shuf == "msm-wsm"]) -
    mean(distances$Distance[pair_shuf == "msw-wsm"])
}

# One-sided p-value
p_value_A <- mean(perm_diff_A >= obs_diff_A)
p_value_A

######################
# B) Within-group dissimilarity MSM vs MSW
######################

# Observed difference: positive = MSM-MSM more similar than MSW-MSW
obs_means <- tapply(distances_men$Distance, distances_men$Behaviour_pair, mean)
obs_diff <- obs_means["MSW-MSW"] - obs_means["MSM-MSM"]

# Permutation
perm_diff <- numeric(n_perm)

for(i in 1:n_perm){
  
  shuffled_labels <- sample(distances_men$Sex_behaviour1,
                            length(distances_men$Sex_behaviour1),
                            replace = FALSE)
  
  Shuf1 <- shuffled_labels[match(distances_men$Sample1, Sample_metadata$SampleID)]
  Shuf2 <- shuffled_labels[match(distances_men$Sample2, Sample_metadata$SampleID)]
  
  type_shuf <- ifelse(Shuf1 == "msm" & Shuf2 == "msm", "msm-msm",
                      ifelse(Shuf1 == "msw" & Shuf2 == "msw", "msw-msw", "msm-msw"))
  
  # MSW-MSW minus MSM-MSM
  perm_diff[i] <- mean(distances_men$Distance[type_shuf == "msw-msw"]) -
    mean(distances_men$Distance[type_shuf == "msm-msm"])
}

# One-sided p-value
p_value <- mean(perm_diff >= obs_diff)
p_value

######################
# Plots
######################

# Subset relevant pairs
distances_A <- distances %>%
  filter(Behaviour_pair %in% c("MSM-WSM", "MSW-WSM"))

# For the plots, represent similarity instead of Distance (1-Distance)
distances_A$Similarity <- 1 - distances_A$Distance
distances_men$Similarity <- 1 - distances_men$Distance

# Downsample for plotting 
distances_A_groups <- distances_A %>% group_by(Behaviour_pair) %>% group_split()
distances_men_groups <- distances_men %>% group_by(Behaviour_pair) %>% group_split()

distances_A_sample <- purrr::map_df(distances_A_groups, ~ {
  n_sample <- min(25000, nrow(.x))
  .x[sample(seq_len(nrow(.x)), n_sample), ]
})

distances_men_sample <- purrr::map_df(distances_men_groups, ~ {
  n_sample <- min(25000, nrow(.x))
  .x[sample(seq_len(nrow(.x)), n_sample), ]
})

# Plot: MSM vs MSW to women
pdf('Plots/Virome_similarity_MSM_vs_MSW_to_women.pdf', width = 3, height = 3.2)
ggplot(distances_A_sample, aes(x = Behaviour_pair, y = Similarity, fill = Behaviour_pair)) +
  geom_jitter(alpha = 0.4, color = "lightgrey", size = 1.5, width = 0.2) +
  geom_boxplot(width = 0.6, color = "black", outlier.color = NA, size = 0.8) +
  labs(x = 'Behaviour pair', y = 'Virome similarity (1 - Bray distance)') +
  scale_fill_manual(values = c("MSM-WSM" = "#8A9A5B", "MSW-WSM" = "#C98B4A")) +
  scale_y_sqrt() +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.75),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    axis.title.y = element_text(margin = margin(r = 10)),
    axis.title.x = element_text(margin = margin(t = 10)),
    legend.position = "none"
  )
dev.off()


distances_men_sample_2cat <- distances_men_sample %>%
  filter(Behaviour_pair %in% c("MSM-MSM", "MSW-MSW"))

pdf('Plots/Virome_similarity_men_sexual_behaviour.pdf', width=2.9, height=3.2)
ggplot(distances_men_sample_2cat, aes(x = Behaviour_pair, y = 1 - Distance, fill = Behaviour_pair)) +
  geom_jitter(alpha = 0.4, color = "lightgrey", size = 1.8, width = 0.2) +
  geom_boxplot(width = 0.6, color = "black", outlier.color = NA, size = 0.9) +
  labs(x = 'Behaviour pair', y = 'Virome similarity (1 - Bray distance)') +
  scale_fill_manual(values = c("MSM-MSM" = "#7A1F3D", "MSW-MSW" = "#6E8CA6")) +
  scale_y_sqrt() +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.75),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    axis.title.y = element_text(margin = margin(r = 10)),
    axis.title.x = element_text(margin = margin(t = 10)),
    legend.position = "none"
  )
dev.off()




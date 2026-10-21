################################################################################
##### 2000HIV: Strain sharing analysis - Strain similarity and sharing events
### Author(s):Asier Fernández-Pato
### Last updated: 14th September, 2026
################################################################################

#****************
# Load modules
#****************
library(dplyr)
library(scales)
library(ggplot2)


# Set working directory 
setwd("~/Desktop/PhD/Projects/HIV 2000/")

# Read abundance table and metadata tables
Abundance_table <- read.delim("ABUNDANCE TABLE/HIV_Abundance_Table_01052025.txt")
Sample_metadata <- read.delim("METADATA/Sample_metadata_excluded_11052026.txt")
Viral_metadata <- read.delim("METADATA/HIV_Viral_Metadata_06052025.txt")

# Load also inStrain results
inStrain_results <- read.delim("8_STRAIN_TRANSMISSION/inStrain_results.txt")
inStrain_results_men <- read.delim("8_STRAIN_TRANSMISSION/inStrain_results_men.txt")

##************************************************************************
# 1. Exploring viral strain distances and sharing: Sex behaviour
#*************************************************************************

######################################################################
# 1.1. Check if MSM have more similar strains than non-MSM men
######################################################################

# Perform statistical test (use permutation-based test due to lack of independence (repeated measures within infants)
num_permutations <- 10000
permuted_stats <- numeric(num_permutations)
set.seed(123)

# extract the observed test statistic: difference in medians
observed_stat <- median(inStrain_results$popANI[inStrain_results$Sex_behaviour_pair == "msm-msm"]) - 
  median(inStrain_results$popANI[inStrain_results$Sex_behaviour_pair == "msw-msw"])

# perform the permutation test
for (i in 1:num_permutations) {
  permuted_labels <- sample(inStrain_results$Sex_behaviour_pair)
  permuted_stat <- median(inStrain_results$popANI[permuted_labels == "msm-msm"]) - 
    median(inStrain_results$popANI[permuted_labels == "msw-msw"])
  permuted_stats[i] <- permuted_stat
}

# one-sided p-value
p_value <- mean(permuted_stats >= observed_stat)

# Generate boxplot of strain similarity in msm-msm vs msw-msw
MSM_strain_data <- inStrain_results[inStrain_results$Sex_behaviour_pair %in% c("msm-msm", "msw-msw"), ]

pdf("8_STRAIN_TRANSMISSION/Plots/MSM_MSM_vs_MSW_MSW_relatedness.pdf", width = 3.2, height = 3.2)
ggplot(MSM_strain_data, aes(x = popANI, fill = Sex_behaviour_pair)) +
  geom_density(alpha = 0.5, color = "black", size = 0.3) +
  geom_vline(data = aggregate(popANI ~ Sex_behaviour_pair, MSM_strain_data, median),
             aes(xintercept = popANI, color = Sex_behaviour_pair),
             linetype = "dashed", size = 0.8, show.legend = FALSE) +
  labs(x = "Strain similarity (popANI)", y = "Density") +
  scale_fill_manual(values = c("msm-msm" = "#7A1F3D", "msw-msw" = "#6E8CA6")) +
  scale_color_manual(values = c("msm-msm" = "#7A1F3D", "msw-msw" = "#6E8CA6")) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.75),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    axis.title.y = element_text(margin = margin(r = 10)),
    axis.title.x = element_text(margin = margin(t = 10)),
    legend.position = c(0.5, 0.75),
    legend.background = element_rect(fill = "white", color = "black", size = 0.3),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
  )
dev.off()


###########################################################################
# 1.2. Check proportion of within sex behavior group pairs with strain sharing
###########################################################################

# Generate all possible pairs
filtered_metadata <- Sample_metadata %>%
  select(SampleID, MSM_derived) %>%
  filter(!is.na(MSM_derived))

all_pairs <- t(combn(filtered_metadata$SampleID, 2)) %>%
  as.data.frame() %>%
  setNames(c("Sample1", "Sample2")) %>%
  mutate(
    Pair_ID = paste(pmin(Sample1, Sample2), pmax(Sample1, Sample2), sep = "_"),
    Sex_behaviour1 = filtered_metadata$MSM_derived[match(Sample1, filtered_metadata$SampleID)],
    Sex_behaviour2 = filtered_metadata$MSM_derived[match(Sample2, filtered_metadata$SampleID)],
    Sex_behaviour_pair = paste(pmin(Sex_behaviour1, Sex_behaviour2),
                               pmax(Sex_behaviour1, Sex_behaviour2), sep = "-"),
    Sharing = Pair_ID %in% inStrain_results$Pair_ID[inStrain_results$Sharing == "Yes"]
  )


# Calculate observed proportions
sharing_summary <- all_pairs %>%
  filter(Sex_behaviour_pair %in% c("msm-msm", "msw-msw")) %>%
  group_by(Sex_behaviour_pair) %>%
  summarise(Proportion = mean(Sharing), .groups = "drop")

observed_stat <- sharing_summary$Proportion[
  sharing_summary$Sex_behaviour_pair == "msm-msm"
] -
  sharing_summary$Proportion[
    sharing_summary$Sex_behaviour_pair == "msw-msw"
  ]

# Individual-level permutation test
msm_ids <- filtered_metadata$SampleID[filtered_metadata$MSM_derived == "msm"]
all_ids <- filtered_metadata$SampleID

num_permutations <- 10000
permuted_stats <- numeric(num_permutations)
set.seed(123)

for (i in 1:num_permutations) {
  
  permuted_msm <- sample(all_ids, length(msm_ids))
  permuted_msw <- setdiff(all_ids, permuted_msm)
  
  permuted_stats[i] <- mean(all_pairs$Sharing[all_pairs$Sample1 %in% permuted_msm &
                        all_pairs$Sample2 %in% permuted_msm]) -
    mean(all_pairs$Sharing[all_pairs$Sample1 %in% permuted_msw &
                          all_pairs$Sample2 %in% permuted_msw])
}

# One-sided p-value
p_value <- mean(permuted_stats >= observed_stat)

observed_stat
p_value

# Generate a plot comparing the proportions
pdf("8_STRAIN_TRANSMISSION/Plots/MSM_MSM_vs_MSW_MSW_sharing_proportion.pdf",
    width = 2.8, height = 3.6)
ggplot(sharing_summary,aes(x = Sex_behaviour_pair, y = Proportion, fill = Sex_behaviour_pair)) +
  geom_bar(stat = "identity", width = 0.8) +
  geom_text(aes(label = paste0(round(Proportion * 100, 1), "%")),
            vjust = -0.5, size = 4.5) +
  scale_fill_manual(values = c("msm-msm" = "#7A1F3D","msw-msw" = "#6E8CA6")) +
  labs(x = "Group",y = "Proportion of pairs\nwith strain sharing") +
  scale_y_continuous(
    labels = scales::percent_format(),
    limits = c(0, max(sharing_summary$Proportion) * 1.2),
    expand = expansion(mult = c(0, 0.02))
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.75),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    legend.position = "none"
  )
dev.off()

##************************************************************************
# 2. Exploring vOTU sharing
#*************************************************************************

###########################################################################
# 2.1. Explore if vOTU sharing is was higher among MSM individuals
###########################################################################

# Convert abundance table to presence/absence 
Abundance_pa <- Abundance_table
Abundance_pa[Abundance_pa > 0] <- 1
votu_ids <- rownames(Abundance_pa)
Abundance_pa_matrix <- as.matrix(Abundance_pa)

# Generate all pairwise sample combinations
sample_pairs <- t(combn(colnames(Abundance_pa_matrix), 2))
sample_pairs <- as.data.frame(sample_pairs)
colnames(sample_pairs) <- c("Sample1", "Sample2")

# Calculate pairwise vOTU overlap
pairwise_votu_sharing <- sample_pairs %>%
  rowwise() %>%
  mutate(
    # Number of shared vOTUs
    shared_vOTUs = sum(Abundance_pa_matrix[, Sample1] == 1 & Abundance_pa_matrix[, Sample2] == 1),
    
    # Total number of detected vOTUs across both individuals
    total_vOTUs = sum(Abundance_pa_matrix[, Sample1] == 1 | Abundance_pa_matrix[, Sample2] == 1),
    
    # Jaccard similarity
    sharing_fraction = shared_vOTUs / total_vOTUs
    
  ) %>%
  ungroup()

# Add metadata
pairwise_votu_sharing <- pairwise_votu_sharing %>%
  mutate(
    Sex_behaviour1 = Sample_metadata$MSM_derived[match(Sample1, Sample_metadata$SampleID)],
    Sex_behaviour2 = Sample_metadata$MSM_derived[match(Sample2, Sample_metadata$SampleID)],
    Sex_behaviour_pair = paste(
      pmin(Sex_behaviour1, Sex_behaviour2),
      pmax(Sex_behaviour1, Sex_behaviour2),
      sep = "-"
    )
  )

# Keep only msm-msm and msw-msw comparisons
pairwise_votu_sharing_filtered <- pairwise_votu_sharing %>%
  filter(Sex_behaviour_pair %in% c("msm-msm", "msw-msw"))

# Summarise overlap statistics
votu_summary <- pairwise_votu_sharing_filtered %>%
  group_by(Sex_behaviour_pair) %>%
  summarise(
    median_sharing = median(sharing_fraction, na.rm = TRUE),
    mean_sharing = mean(sharing_fraction, na.rm = TRUE),
    median_shared_vOTUs = median(shared_vOTUs, na.rm = TRUE),
    mean_shared_vOTUs = mean(shared_vOTUs, na.rm = TRUE),
    n_pairs = n(),
    .groups = "drop"
  )

votu_summary

####################################################
# Permutation test: MSM vs MSW vOTU sharing fraction
####################################################

num_permutations <- 10000
set.seed(123)

# Observed difference in medians
observed_stat <- median(pairwise_votu_sharing_filtered$sharing_fraction[
    pairwise_votu_sharing_filtered$Sex_behaviour_pair == "msm-msm"], na.rm = TRUE) -
  median(pairwise_votu_sharing_filtered$sharing_fraction[
    pairwise_votu_sharing_filtered$Sex_behaviour_pair == "msw-msw"],na.rm = TRUE)

# Permutation test
permuted_stats <- numeric(num_permutations)

for (i in 1:num_permutations) {
  permuted_labels <- sample(pairwise_votu_sharing_filtered$Sex_behaviour_pair)
  
  permuted_stats[i] <- median(pairwise_votu_sharing_filtered$sharing_fraction[
      permuted_labels == "msm-msm"],na.rm = TRUE) -
    median(pairwise_votu_sharing_filtered$sharing_fraction[permuted_labels == "msw-msw"],
      na.rm = TRUE)
}

# One-sided p-value
p_value_votu <- mean(permuted_stats >= observed_stat,na.rm = TRUE)
p_value_votu

####################################################
# Plot
####################################################

votu_plot <- pairwise_votu_sharing_filtered %>% filter(is.finite(sharing_fraction))

pdf("8_STRAIN_TRANSMISSION/Plots/vOTU_sharing_MSM_vs_MSW.pdf", width = 2.8, height = 3.4)
ggplot(votu_plot, aes(Sex_behaviour_pair, sharing_fraction, fill = Sex_behaviour_pair)) +
  geom_violin(alpha = 0.35, color = NA, trim = FALSE) +
  geom_boxplot(width = 0.2, color = "black", linewidth = 0.4,
               outlier.shape = NA, alpha = 0.85) +
  scale_fill_manual(values = c("msm-msm" = "#7A1F3D", "msw-msw" = "#6E8CA6")) +
  scale_y_continuous(labels = scales::percent_format(),
                     limits = c(0, 0.10), breaks = seq(0, 0.10, 0.02)) +
  labs(x = "Sex behaviour pair",y = "vOTU sharing (Jaccard similarity)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.75),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.title.x = element_text(margin = margin(t = 10)),
        legend.position = "none")
dev.off()


##************************************************************************
# 3. Explore strain sharing events across viral taxa
#*************************************************************************

# Retrieve list of shared vOTUs and their number of sharing events
inStrain_results_sharing <- inStrain_results[inStrain_results$Sharing == "Yes",]
shared_vOTUs <- unique(unlist(inStrain_results_sharing[, "Virus_ID"])) # 1443
sharing_events_vOTU <- data.frame(table(inStrain_results_sharing$Virus_ID))
colnames(sharing_events_vOTU) <- c("Virus_ID", "N_sharing_events")

# Add to viral metadata
Viral_metadata <- left_join(Viral_metadata, sharing_events_vOTU, by = "Virus_ID")
Viral_metadata$N_sharing_events[is.na(Viral_metadata$N_sharing_events)] <- 0

# Exclude HIV_91378, flagged as a potential contaminant (~20% of all sharing events)
sharing_events_vOTU <- sharing_events_vOTU %>%
  filter(Virus_ID != "HIV_91378")

# Sort vOTUs by number of sharing events
sharing_events_ranked <- sharing_events_vOTU %>%
  arrange(desc(N_sharing_events)) %>%
  mutate(
    Rank = row_number(),
    Cumulative_fraction = cumsum(N_sharing_events) / sum(N_sharing_events)
  )

# Add host information to vOTUs
sharing_events_ranked$Bacterial_genus_host <- 
  Virus_metadata$Bacterial_genus_host[
    match(sharing_events_ranked$Virus_ID, Virus_metadata$Virus_ID)
  ]

sharing_events_ranked$Bacterial_family_host <- 
  Virus_metadata$Bacterial_family_host[
    match(sharing_events_ranked$Virus_ID, Virus_metadata$Virus_ID)
  ]

# Plot cumulative fraction of sharing events vs number of vOTUs
pdf("8_STRAIN_TRANSMISSION/Plots/core_transmissible_virome_cumulative.pdf", width = 4.5, height = 3.8)
ggplot(sharing_events_ranked, aes(x = Rank, y = Cumulative_fraction)) +
  geom_line(size = 1, color = "#800020") +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey40") +
  scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1)) +
  labs( x = "vOTU rank", y = "Cumulative fraction of sharing events") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.75),
    axis.text = element_text(size = 13),
    axis.title = element_text(size = 15)
  )
dev.off()

# Plot the top host genus viral groups with more sharing events
top_hosts_sharing <- Viral_metadata %>%
  filter(N_sharing_events > 0) %>%
  mutate(
    Bacterial_genus_host = ifelse(
      is.na(Bacterial_genus_host),
      "Unknown",
      Bacterial_genus_host
    )
  ) %>%
  group_by(Bacterial_genus_host) %>%
  summarise(
    Total_sharing_events = sum(N_sharing_events),
    .groups = "drop"
  ) %>%
  arrange(desc(Total_sharing_events)) %>%
  slice_head(n = 15)

pdf('8_STRAIN_TRANSMISSION/Plots/top15_host_genera_sharing_events.pdf', width = 10, height = 9)
ggplot(top_hosts_sharing,aes(
    x = reorder(Bacterial_genus_host, Total_sharing_events),
    y = Total_sharing_events,
    fill = Bacterial_genus_host)) +
  geom_bar(stat = "identity", width = 0.9) +
  geom_text(aes(label = paste0("n=", Total_sharing_events)),
    hjust = -0.1, color = "black", size = 4.3) +
  coord_flip() +
  labs(x = "Predicted bacterial host genus",
    y = "Number of strain-sharing events", fill = "Host genus") +
  scale_fill_manual(
    values = scales::alpha(c(
      "#1B4965", "#5B8E7D", "#C17C74", "#7B6D8D", "#D4A373",
      "#6B705C", "#9A78B8", "#A68A64", "#52796F", "#8D99AE",
      "#A67C8E", "#7A8F9A", "#B08B6F", "#857A9B", "#739681"
    ), 0.85)
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.75),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    axis.title.y = element_text(margin = margin(r = 10)),
    axis.title.x = element_text(margin = margin(t = 10)),
    legend.position = "none"
  )
dev.off()

##################################################################
# Statistical test: are strain sharing events enriched among temperate phages?
##################################################################

temperate_test <- Viral_metadata %>% 
  filter(!is.na(Lifestyle)) %>%
  mutate(Sharing = N_sharing_events > 0)

contingency_table_temperate <- table(
  temperate_test$Lifestyle,
  temperate_test$Sharing)

# Chi-square test
chi2_temperate <- chisq.test(contingency_table_temperate)
chi2_temperate

# Generate a plot comparing the proportions
sharing_summary <- temperate_test %>%
  group_by(Lifestyle) %>%
  summarise(
    Proportion = mean(Sharing),
    .groups = "drop"
  )

pdf("8_STRAIN_TRANSMISSION/Plots/temperate_vs_virulent_sharing_proportion.pdf", width = 2.8, height = 3.6)
ggplot(sharing_summary, aes(x = Lifestyle, y = Proportion, fill = Lifestyle)) +
  geom_bar(stat = "identity", width = 0.8) +
  geom_text(aes(label = paste0(round(Proportion * 100, 1), "%")),
    vjust = -0.5, size = 4.5) +
  scale_fill_manual(values = c("Temperate" = "#5C6BC0", "Virulent" = "#FFD54F")) + #"Temperate" = "#1B4965", "Virulent" = "#C17C74"
  labs(x = "Phage lifestyle", y = "Proportion of phages with strain sharing") +
  scale_y_continuous(labels = scales::percent_format(),
    limits = c(0, max(sharing_summary$Proportion) * 1.2),
    expand = expansion(mult = c(0, 0.02))) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.75),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    legend.position = "none"
  )
dev.off() 

##************************************************************************
# Write results
##************************************************************************

write.table(sharing_events_ranked,
  file = "8_STRAIN_TRANSMISSION/Sharing_events_ranked.txt", sep = "\t",
  row.names = FALSE, quote = FALSE)

write.table(Sample_metadata, file = "METADATA/Sample_metadata_excluded_11052026.txt",
            sep = "\t", row.names = F, quote = FALSE)


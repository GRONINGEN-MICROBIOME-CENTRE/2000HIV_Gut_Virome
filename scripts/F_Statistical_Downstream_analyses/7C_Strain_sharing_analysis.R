################################################################################
##### 2000HIV: Strain sharing analysis - Connectivity and HIV phylogeny
### Author(s):Asier Fernández-Pato
### Last updated: 15th September, 2026
################################################################################

#****************
# Load modules
#****************
library(dplyr)
library(scales)
library(ggplot2)
library(stringr)
library(emmeans)
library(igraph)
library(tidygraph)
library(ggraph)


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
# 1. Explore viral strain-sharing network connectivity
#*************************************************************************

###########################################################################
# 1.1. Check if Sex behavior and N partners associate with connectivity
###########################################################################

# Collapse strain-sharing events per pair
strain_sharing_network <- inStrain_results %>%
  filter(Sharing == "Yes") %>%
  distinct(Pair_ID, Sample1, Sample2)

# Convert to undirected edge list
network_edges_1 <- strain_sharing_network %>%
  select(Sample = Sample1, Partner = Sample2)

network_edges_2 <- strain_sharing_network %>%
  select(Sample = Sample2, Partner = Sample1)

network_edges <- bind_rows(network_edges_1, network_edges_2)

# Count number of unique sharing partners per individual
network_connectivity <- network_edges %>%
  group_by(Sample) %>%
  summarise(
    n_connections = n_distinct(Partner),
    .groups = "drop"
  )

# Add metadata
network_connectivity <- network_connectivity %>%
  mutate(
    Sex_behaviour = as.factor(Sample_metadata$MSM_derived[match(Sample, Sample_metadata$SampleID)]),
    Sexual_partner = as.factor(Sample_metadata$Sexual_partner[match(Sample, Sample_metadata$SampleID)]),
    Vaginal_sex = as.factor(Sample_metadata$Vaginal_sexual_contact[match(Sample, Sample_metadata$SampleID)]),
    Anal_sex = as.factor(Sample_metadata$Anal_sexual_contact[match(Sample, Sample_metadata$SampleID)]),
    Receptive_anal_sex = as.factor(Sample_metadata$Receptive_anal_sex[match(Sample, Sample_metadata$SampleID)]),
    Sexual_orientation = as.factor(Sample_metadata$Sexual_orientation[match(Sample, Sample_metadata$SampleID)]),
    Read_depth = Sample_metadata$read_depth[match(Sample, Sample_metadata$SampleID)]
  )

# Summary stats
connectivity_summary <- network_connectivity %>%
  group_by(Sex_behaviour) %>%
  summarise(
    median_connections = median(n_connections, na.rm = TRUE),
    mean_connections = mean(n_connections, na.rm = TRUE),
    n_individuals = n(),
    .groups = "drop"
  )

connectivity_summary


############################################################################
# Statistical test: compare connectivity according to sex behaviour
############################################################################

###################
# A) Wilcoxon test (not included)
###################

wilcox_sex_beh_connectivity <- pairwise.wilcox.test(
  network_connectivity$n_connections,
  network_connectivity$Sex_behaviour,
  p.adjust.method = "BH"
)
wilcox_sex_beh_connectivity

wilcox_partners <- pairwise.wilcox.test(
  network_connectivity$n_connections,
  network_connectivity$Sexual_partner,
  p.adjust.method = "BH"
)
wilcox_partners

###################
# B) GLM 
###################

dat <- network_connectivity %>%
  filter(!is.na(Sex_behaviour), !is.na(Sexual_partner)) %>%
  mutate(Sex_behaviour = factor(Sex_behaviour, levels = c("msm", "msw", "wsm")))

m_sex <- glm(n_connections ~ Sex_behaviour + Read_depth, data = dat)
m_partner <- glm(n_connections ~ Sexual_partner + Read_depth, data = dat)
m_adj <- glm(n_connections ~ Sex_behaviour + Sexual_partner + Read_depth, data = dat)

pairs(emmeans(m_sex, ~ Sex_behaviour), adjust = "BH")
pairs(emmeans(m_partner, ~ Sexual_partner), adjust = "BH")
pairs(emmeans(m_adj, ~ Sex_behaviour), adjust = "BH")
summary(pairs(pairs(emmeans(m_adj, ~ Sexual_partner), adjust = "BH")))$p.value

###################
# C) GLM (only for men)
###################

dat_msm <- network_connectivity %>%
  filter(Sex_behaviour == "msm", !is.na(Anal_sex)) %>%
  mutate(Anal_sex = factor(Anal_sex, levels = 0:3))

m_anal <- glm(n_connections ~ Receptive_anal_sex + Read_depth, data = dat_msm)
m_anal_adj <- glm(n_connections ~ Sexual_partner + Receptive_anal_sex + Read_depth, data = dat_msm)

summary(pairs(emmeans(m_anal, ~ Receptive_anal_sex), adjust = "BH"))$p.value
summary(pairs(emmeans(m_anal_adj, ~ Receptive_anal_sex), adjust = "BH"))$p.value

###################
# Plots
###################

# Generate plots showing connectivity according to sex behaviour and n partners
pdf("8_STRAIN_TRANSMISSION/Plots/strain_sharing_network_connectivity.pdf", width = 3.2, height = 3.4)
ggplot(network_connectivity %>% filter(Sex_behaviour %in% c("msm", "msw", "wsm")),
       aes(x = Sex_behaviour, y = n_connections, fill = Sex_behaviour)) +
  geom_jitter(color = "lightgrey", alpha = 0.4, size = 1.8, width = 0.2) +
  geom_boxplot(width = 0.6, color = "black", linewidth = 0.4, outlier.shape = NA, alpha = 0.85) +
  scale_fill_manual(values = c("msm" = "#7A1F3D", "msw" = "#6E8CA6", "wsm" = "#8A9A5B")) +
  labs(x = "Sex behaviour", y = "Number of strain-sharing connections") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.75),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.position = "none")
dev.off()

# Barplot with the sexual behaviour distrbution of strain sharing connections
connection_types <- network_edges %>%
  mutate(
    Sex_behaviour = Sample_metadata$MSM_derived[match(Sample, Sample_metadata$SampleID)],
    Partner_behaviour = Sample_metadata$MSM_derived[match(Partner, Sample_metadata$SampleID)]
  ) %>%
  filter(
    !is.na(Sex_behaviour),
    !is.na(Partner_behaviour)
  )

connection_summary <- connection_types %>%
  count(Sex_behaviour, Partner_behaviour) %>%
  group_by(Sex_behaviour) %>%
  mutate(proportion = n / sum(n)) %>%
  ungroup()


pdf("8_STRAIN_TRANSMISSION/Plots/strain_sharing_network_connectivity_partners.pdf",
    width = 3.2, height = 3.4)
ggplot(network_connectivity %>% filter(!is.na(Sexual_partner)),
       aes(x = factor(Sexual_partner), y = n_connections, fill = factor(Sexual_partner))) +
  geom_jitter(color = "lightgrey", alpha = 0.4, size = 1.8, width = 0.2) +
  geom_boxplot(width = 0.6, color = "black", linewidth = 0.4, outlier.shape = NA, alpha = 0.85) +
  scale_fill_manual(values = c("0" = "#5DA5A4", "1" = "#D4A64A", "2" = "#B65F3D")) +
  labs(x = "Number of sexual partners", y = "Number of strain-sharing connections") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.75),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.position = "none")

dev.off()

pdf("8_STRAIN_TRANSMISSION/Plots/strain_sharing_network_connectivity_partners_by_MSM.pdf",
    width = 6.5, height = 3.4)
ggplot(network_connectivity %>% filter(!is.na(Sexual_partner),
                                       Sex_behaviour %in% c("msm", "msw", "wsm")),
       aes(x = factor(Sexual_partner), y = n_connections, fill = factor(Sexual_partner))) +
  geom_jitter(color = "lightgrey", alpha = 0.4, size = 1.8, width = 0.2) +
  geom_boxplot(width = 0.6, color = "black", linewidth = 0.4, outlier.shape = NA, alpha = 0.85) +
  scale_fill_manual(values = c("0" = "#5DA5A4", "1" = "#D4A64A", "2" = "#B65F3D")) +
  facet_wrap(~Sex_behaviour) +
  labs(x = "Number of sexual partners", y = "Number of strain-sharing connections") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.75),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    strip.text = element_text(size = 14),
    legend.position = "none"
  )
dev.off()

# Generate plot showing the connections among MSM according to receptive anal sex
pdf("8_STRAIN_TRANSMISSION/Plots/strain_sharing_network_connectivity_MSM_receptive_anal_sex.pdf",
    width = 2.9, height = 3.4)
ggplot(network_connectivity %>%
    filter(Sex_behaviour == "msm", !is.na(Receptive_anal_sex)),
  aes(x = Receptive_anal_sex, y = n_connections, fill = Receptive_anal_sex)) +
  geom_jitter(color = "lightgrey", alpha = 0.4, size = 1.8, width = 0.2) +
  geom_boxplot(width = 0.6, color = "black", linewidth = 0.4,
               outlier.shape = NA, alpha = 0.85) +
  scale_fill_manual(values = c("No" = "#7B6D8D", "Yes" = "#5B8068")) +
  labs(x = "Receptive anal sex", y = "Number of strain-sharing connections") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.75),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    legend.position = "none"
  )
dev.off()


pdf("8_STRAIN_TRANSMISSION/Plots/strain_sharing_network_connectivity_MSM_receptive_anal_sex_by_partners.pdf",
    width=6.5,height=3.4)
ggplot(network_connectivity %>%
         filter(Sex_behaviour == "msm", !is.na(Receptive_anal_sex), !is.na(Sexual_partner)),
       aes(Receptive_anal_sex, n_connections, fill = Receptive_anal_sex)) +
  geom_jitter(color = "lightgrey", alpha = 0.4, size = 1.8, width = 0.2) +
  geom_boxplot(width = 0.6, color = "black", linewidth = 0.4, outlier.shape = NA, alpha = 0.85) +
  scale_fill_manual(values = c("No" = "#7B6D8D", "Yes" = "#5B8068")) +
  facet_wrap(~factor(Sexual_partner, levels = c(0, 1, 2), labels = c("0", "1", "≥2"))) +
  labs(x = "Receptive anal sex", y = "Number of strain-sharing connections") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.75),
        axis.text = element_text(size = 14), axis.title = element_text(size = 16),
        strip.text = element_text(size = 14), legend.position = "none")
dev.off()

##************************************************************************
# 2. Generate a circular graph to visualize strain-sharing connections
#*************************************************************************

######################
# A. Sexual behaviour
#####################

g <- graph_from_data_frame(network_edges %>% distinct(Sample, Partner), directed = FALSE)

# Add metadata
V(g)$Sex_behaviour <- Sample_metadata$MSM_derived[match(V(g)$name, Sample_metadata$SampleID)]
V(g)$Sexual_partner <- Sample_metadata$Sexual_partner[match(V(g)$name, Sample_metadata$SampleID)]

# Remove isolates
g <- delete_vertices(g, V(g)[degree(g) == 0])

# Convert to tbl_graph + compute degree properly
g_tbl <- as_tbl_graph(g) %>%
  mutate(deg = centrality_degree()) %>%
  activate(nodes) %>%
  filter(deg >= 2) %>%
  filter(!is.na(Sex_behaviour))

# Add edge metadata
g_tbl <- g_tbl %>%
  activate(edges) %>%
  mutate(from_group = .N()$Sex_behaviour[from],
         to_group = .N()$Sex_behaviour[to],
         connection_type = case_when(
           from_group == "msm" & to_group == "msm" ~ "msm-msm",
           from_group == "msw" & to_group == "msw" ~ "msw-msw",
           from_group == "wsm" & to_group == "wsm" ~ "wsm-wsm",
           (from_group == "msw" & to_group == "wsm") | (from_group == "wsm" & to_group == "msw") ~ "msw-wsm",
           TRUE ~ "other"
         ))

# Order nodes and add circular layout
g_tbl <- g_tbl %>%
  activate(nodes) %>%
  arrange(factor(Sex_behaviour, levels = c("msm","msw","wsm")))

layout_circle <- create_layout(g_tbl, layout = "linear", circular = TRUE) %>%
  mutate(
    x_bar_start = x * 1.025,
    y_bar_start = y * 1.025,
    x_bar_end = x * (1.025 + deg / max(deg) * 0.22),
    y_bar_end = y * (1.025 + deg / max(deg) * 0.22)
  )

pdf("8_STRAIN_TRANSMISSION/Plots/strain_sharing_network_grouped_circular_degree.pdf", width = 8, height = 8)
ggraph(layout_circle) +
  geom_edge_arc(aes(colour = connection_type), alpha = 0.25, linewidth = 0.3) +
  geom_segment(data = layout_circle,
               aes(x = x_bar_start, y = y_bar_start, xend = x_bar_end, yend = y_bar_end),
               colour = "grey40", alpha = 0.7, linewidth = 0.45) +
  geom_node_point(aes(color = Sex_behaviour), size = 2, alpha = 0.85) +
  scale_color_manual(values = c(msm = "#7A1F3D", msw = "#6E8CA6", wsm = "#8A9A5B")) +
  scale_edge_colour_manual(values = c("msm-msm" = "#7A1F3D", "msw-msw" = "#6E8CA6",
                                      "wsm-wsm" = "#8A9A5B", "msw-wsm" = "#C98B4A",
                                      other = "grey85")) +
  coord_equal(clip = "off") +
  theme_void() +
  theme(legend.position = "right",
        plot.background = element_blank(),
        plot.margin = margin(10, 20, 10, 10))

dev.off()

######################
# B. Anal sex (MSM only)
#####################
msm_ids <- Sample_metadata$SampleID[Sample_metadata$MSM_derived == "msm"]

g_anal <- graph_from_data_frame(
  network_edges %>%
    filter(Sample %in% msm_ids, Partner %in% msm_ids) %>%
    distinct(Sample, Partner),
  directed = FALSE
)

V(g_anal)$Receptive_anal_sex <- Sample_metadata$Receptive_anal_sex[
  match(V(g_anal)$name, Sample_metadata$SampleID)
]
V(g_anal)$Sexual_partner <- Sample_metadata$Sexual_partner[
  match(V(g_anal)$name, Sample_metadata$SampleID)
]

g_anal_tbl <- as_tbl_graph(g_anal) %>%
  mutate(deg = centrality_degree()) %>%
  activate(nodes) %>%
  mutate(
    Receptive_anal_sex = factor(Receptive_anal_sex, levels = c("No", "Yes")),
    Sexual_partner = factor(Sexual_partner,
                            levels = c(0, 1, 2),
                            labels = c("0", "1", "≥2"))
  ) %>%
  filter(deg >= 2, !is.na(Receptive_anal_sex), !is.na(Sexual_partner)) %>%
  arrange(Receptive_anal_sex) %>%
  activate(edges) %>%
  mutate(
    from_group = .N()$Receptive_anal_sex[from],
    to_group = .N()$Receptive_anal_sex[to],
    connection_type = case_when(
      from_group == "No" & to_group == "No" ~ "No–No",
      from_group == "Yes" & to_group == "Yes" ~ "Yes–Yes",
      TRUE ~ "Yes–No"
    )
  )

layout_anal <- create_layout(g_anal_tbl, layout = "linear", circular = TRUE) %>%
  mutate(
    # 1. Receptive anal sex
    x_anal_start = x * 1.001, y_anal_start = y * 1.001,
    x_anal_end   = x * 1.031, y_anal_end   = y * 1.031,
    
    # 2. Sexual partners
    x_partner_start = x * 1.035, y_partner_start = y * 1.035,
    x_partner_end   = x * 1.065, y_partner_end   = y * 1.065,
    
    # 3. Total connections
    x_bar_start = x * 1.080, y_bar_start = y * 1.080,
    x_bar_end = x * (1.080 + deg / max(deg) * 0.22),
    y_bar_end = y * (1.080 + deg / max(deg) * 0.22)
    )

pdf("8_STRAIN_TRANSMISSION/Plots/strain_sharing_network_MSM_receptive_anal_partners_circular_degree.pdf",
  width = 8, height = 8)
ggraph(layout_anal) +
  geom_edge_arc(aes(colour = connection_type), alpha = 0.25, linewidth = 0.3) +
  
  # Receptive anal sex
  geom_segment(data = layout_anal, aes(x = x_anal_start, y = y_anal_start,
        xend = x_anal_end, yend = y_anal_end, colour = Receptive_anal_sex),
    linewidth = 1.5, lineend = "butt") +
  scale_color_manual(name = "Receptive anal sex", values = c("No" = "#7B6D8D", "Yes" = "#5B8068")) +
  ggnewscale::new_scale_color() +
  
  # Sexual partners
  geom_segment(data = layout_anal,aes(x = x_partner_start, y = y_partner_start,
        xend = x_partner_end, yend = y_partner_end, colour = Sexual_partner),
    linewidth = 1.5, lineend = "butt") +
  scale_color_manual(name = "Sexual partners",
    values = c("0" = "#5B8E8E", "1" = "#D19A4A", "≥2" = "#B65F3D")) +
  
  # Total strain-sharing connections
  geom_segment(data = layout_anal,
    aes(x = x_bar_start, y = y_bar_start, xend = x_bar_end, yend = y_bar_end),
    colour = "grey40", alpha = 0.7, linewidth = 0.45) +
  scale_edge_colour_manual(name = "Connection type",
                           values = c("No–No" = "#7B6D8D", "Yes–No" = "grey75",
                                      "Yes–Yes" = "#5B8068")) +
  coord_equal(clip = "off") +
  theme_void() +
  theme(
    legend.position = "right",
    plot.background = element_blank(),
    plot.margin = margin(10, 25, 10, 10)
  )
dev.off()

##************************************************************************
# 2. HIV phylogenetic clusters
#*************************************************************************

# Load HIV phylogeny information
linking_file <- read_excel("METADATA/HIV Phylogeny/HCRC_IDs.xlsx")
HIV_phylogeny <- read_excel("METADATA/HIV Phylogeny/HFGP2000_pol_matching_tree_Asier_ClusterPicker.xlsx", sheet = 2)
linking_file$new_ID <- sub("^HFGP", "", linking_file$HCRC_ID)

# Select only individuals included in our study
linking_file <- linking_file[linking_file$ID %in% Sample_metadata$SampleID,]

# Generate summary stats
n_clusters <- nrow(HIV_phylogeny) # Total number of clusters
n_10_individual_clusters <- sum(HIV_phylogeny$NumberOfTips > 5, na.rm = TRUE) # Clusters with >1 individual
HIV_phylogeny$GD <- as.numeric(HIV_phylogeny$GD)
HIV_phylogeny$NumberOfTips <- as.numeric(HIV_phylogeny$NumberOfTips)
mean_GD <- mean(HIV_phylogeny$GD, na.rm = TRUE) # Mean intra-cluster genetic distance
median_GD <- median(HIV_phylogeny$GD, na.rm = TRUE)
quantile_GD <- quantile(HIV_phylogeny$GD, probs = c(0.25, 0.75), na.rm = TRUE)
summary_cluster_sizes <- summary(HIV_phylogeny$NumberOfTips) #cluster size summary

###########################################################################
# 7.1. Check association with sexual behavior
###########################################################################

# Format HIV phylogeny clusters (one row per individual)
HIV_phylogeny_expanded <- HIV_phylogeny %>%
  mutate(Tip_vector = strsplit(gsub("\\[|\\]", "", TipNames), ",\\s*")) %>%
  tidyr::unnest(Tip_vector) %>%
  mutate(Tip_vector = trimws(Tip_vector)) %>%
  rename(new_ID = Tip_vector)

# Add sample IDs and sexual behaviour info
HIV_phylogeny_expanded <- HIV_phylogeny_expanded %>%
  left_join(
    linking_file %>% select(ID, new_ID),
    by = "new_ID"
  )

HIV_phylogeny_expanded <- HIV_phylogeny_expanded %>%
  left_join(
    Sample_metadata %>% select(SampleID, MSM_derived),
    by = c("ID" = "SampleID")
  )

# Summarize composition of each cluster
cluster_structure <- HIV_phylogeny_expanded %>%
  group_by(ClusterNumber) %>%
  summarise(
    msm = any(MSM_derived == "msm", na.rm = TRUE),
    msw = any(MSM_derived == "msw", na.rm = TRUE),
    wsm = any(MSM_derived == "wsm", na.rm = TRUE),
    n_individuals = n(),
    behaviours = paste(sort(unique(MSM_derived)), collapse = "-"),
    .groups = "drop"
  ) %>%
  mutate(
    msm_only = msm & !msw & !wsm
  )


# Add distance information
cluster_structure <- cluster_structure %>%
  left_join(
    HIV_phylogeny_expanded %>%
      mutate(ClusterNumber = as.numeric(ClusterNumber)) %>%
      select(ClusterNumber, GD) %>%
      distinct(),
    by = "ClusterNumber"
  ) %>%
  arrange(ClusterNumber)

# Count cluster types
table(cluster_structure$behaviours)

# Perform permutation test
num_permutations <- 10000
set.seed(122)

observed_stat <- sum(cluster_structure$msm & !cluster_structure$msw & !cluster_structure$wsm)
permuted_stats <- numeric(num_permutations)

for (i in 1:num_permutations) {
  permuted_labels <- sample(HIV_phylogeny_expanded$MSM_derived)
  tmp <- HIV_phylogeny_expanded
  tmp$permuted <- permuted_labels
  
  perm_cluster <- aggregate(permuted ~ ClusterNumber, data = tmp, FUN = function(x) {
    msm <- any(x == "msm")
    msw <- any(x == "msw")
    wsm <- any(x == "wsm")
    msm & !msw & !wsm
  })
  
  permuted_stats[i] <- sum(perm_cluster$permuted)
}

# p-value
p_value_msm_clusters <- mean(permuted_stats >= observed_stat, na.rm = TRUE)
p_value_msm_clusters

# Generate a histogram with the results
perm_df <- data.frame(value = permuted_stats)
pdf("8_STRAIN_TRANSMISSION/Plots/Histogram_HIV_phylogeny_MSMonly.pdf",width = 3.2, height = 3.4) 
ggplot(perm_df, aes(x = value)) +
  geom_histogram(bins = 30, fill = "#B0B0B0", color = "white") +
  geom_vline(xintercept = observed_stat,
             color = "#4B4FC5", size = 1.2, linetype = "dashed") +
  theme_bw() +
  labs(x = "Number of MSM-only clusters", y = "Permutation frequency") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.75),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14)
  )
dev.off()

# Generate also a barplot with the distribution of clusters - sex behavior
ggplot(cluster_structure %>% filter(!is.na(behaviours)),
       aes(x = reorder(behaviours, behaviours, function(x) -length(x)))) +
  geom_bar(fill = "#4CA6B1", width = 0.8) +
  theme_bw() +
  labs(x = "Cluster sexual behavior composition", y = "Number of clusters") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.75),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 12, hjust = 1),
    axis.title = element_text(size = 14)
  )


###########################################################################
# 7.2. Check if sharing events are more common within clusters
###########################################################################

# We restrict this analysis only to MSM individuals

# Add phylogenetic cluster info to inStrain results
inStrain_results <- inStrain_results %>%
  mutate(
    HIV_cluster1 = HIV_phylogeny_expanded$ClusterNumber[match(Sample1, HIV_phylogeny_expanded$ID)],
    HIV_cluster2 = HIV_phylogeny_expanded$ClusterNumber[match(Sample2, HIV_phylogeny_expanded$ID)],
  )

inStrain_results <- inStrain_results %>%
  mutate(
    same_HIV_cluster = HIV_cluster1 == HIV_cluster2
  )

# Collapse sharing results per pair for intra-cluster pairs
# Select only MSM
within_cluster_pairs <- inStrain_results %>%
  filter(same_HIV_cluster == TRUE & both_msm == TRUE) %>%
  group_by(Pair_ID) %>%
  summarise(
    Sharing = any(Sharing == "Yes"),
    .groups = "drop"
  )

# Collapse sharing results per pair for inter-cluster pairs
between_cluster_pairs <- inStrain_results %>%
  filter(same_HIV_cluster == FALSE  & both_msm == TRUE) %>%
  group_by(Pair_ID) %>%
  summarise(
    Sharing = any(Sharing == "Yes"),
    .groups = "drop"
  )

# Combine groups
sharing_summary_clusters <- bind_rows(
  within_cluster_pairs %>% mutate(Group = "Within cluster"),
  between_cluster_pairs %>% mutate(Group = "Between clusters")
)

# Build contingency table
contingency_table_clusters <- table(
  sharing_summary_clusters$Group,
  sharing_summary_clusters$Sharing
)

contingency_table_clusters

# Perform chi-square test
chi2_test_clusters <- chisq.test(contingency_table_clusters)
chi2_test_clusters

####################
# Generate the plot
####################

contingency_table_clusters_plot <- as.data.frame(as.table(contingency_table_clusters))
colnames(contingency_table_clusters_plot) <- c("Group", "Strain_sharing", "Count")

# Compute proportions
contingency_table_clusters_plot <- contingency_table_clusters_plot %>%
  group_by(Group) %>%
  mutate(Proportion = Count / sum(Count))

# Set plotting order
contingency_table_clusters_plot$Group <- factor(
  contingency_table_clusters_plot$Group,
  levels = c("Within cluster", "Between clusters")
)

# Generate stacked bar plot
pdf('8_STRAIN_TRANSMISSION/Plots/Within_vs_Between_HIV_cluster_sharing_events.pdf', width = 2.8, height = 3.6)
ggplot(contingency_table_clusters_plot, aes(x = Group, y = Proportion, fill = Strain_sharing)) +
  geom_bar(stat = "identity", width = 0.9) +
  geom_text(aes(label = paste0("n=", Count)), position = position_stack(vjust = 0.5), color = "black", size = 4.5) +
  scale_fill_manual(values = c("FALSE" = "#E0E0E0", "TRUE" = "#800020")) +
  labs(x = "Pair type", y = "Proportion of pairs", fill = "Strain sharing") +
  scale_y_continuous(labels = scales::percent_format(),limits = c(0, 1.15),
                     breaks = c(0, 0.25, 0.5, 0.75, 1), expand = expansion(add = c(0.05, 0.05))) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 0.75),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.position = "top",
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14)
  )
dev.off()

##************************************************************************
# Write results
##************************************************************************

write.table(cluster_structure, file = "8_STRAIN_TRANSMISSION//HIV_phylogeny_cluster_structure.txt",
  sep = "\t", row.names = FALSE, quote = FALSE)

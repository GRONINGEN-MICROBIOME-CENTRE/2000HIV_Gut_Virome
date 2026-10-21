################################################################################
##### 2000HIV: vOTU diversity/composition and association plots
### Author(s):Asier Fernández-Pato
### Last updated: 25th September, 2026
################################################################################

#****************
# Load modules
#****************
library(vegan)
library(dplyr)
library(ggplot2)
library(ggExtra)
library(tibble)
library(scales)
library(gggenes)

# Set working directory 
setwd("~/Desktop/PhD/Projects/HIV 2000/")

# Read abundance table and metadata tables
Abundance_table <- read.delim("ABUNDANCE TABLE/HIV_Abundance_Table_01052025.txt")
Abundance_table_host <- read.delim("5_HOST_PREDICTION/phage_abundance_by_host_genus_relative_abundance.tsv", row.names = 1)
Abundance_table_bact <- read.delim("ABUNDANCE TABLE/Bacterial_species_abundance_11052026.txt")
Abundance_table_bact_genus <- read.delim("ABUNDANCE TABLE/Bacterial_genus_abundance_11052026.txt")
Sample_metadata <- read.delim("METADATA/Sample_metadata_excluded_11052026.txt")
Virus_metadata <- read.delim("METADATA/HIV_Viral_Metadata_06052025.txt")

# Remove excluded participants
Abundance_table <- Abundance_table [, colnames(Abundance_table) %in% Sample_metadata$SampleID]

##************************************************************************
# 1. Alpha and beta diversity plots
#*************************************************************************

# First, generate a general plot with the number of sexual behaviour categories
pdf('2_DIVERSITY/Plots/MSM_derived_counts.pdf', width = 2.4, height = 3.0)
ggplot(Sample_metadata[!is.na(Sample_metadata$MSM_derived), ], aes(x = MSM_derived, fill = MSM_derived)) +
  geom_bar(width = 0.9, alpha = 0.9) +
  labs(x = "Sexual behaviour", y = "Number of samples") +
  scale_fill_manual(values = c("msm" = "#ec9060","msw" = "#a084ba","wsm" = "#619ac0")) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 15),
    axis.title.y = element_text(margin = margin(r = 10)),
    axis.title.x = element_text(margin = margin(t = 10)),
    legend.position = "none"
  )
dev.off()

##################################
# 1.1. Estimation of viral alpha diversity and richness in relation to sexual behaviour
##################################

# Estimate richness and alpha diversity (Shannon)
Sample_metadata$richness  <-specnumber(t(Abundance_table))                 
Sample_metadata$shannon <- vegan::diversity(t(Abundance_table), index = "shannon") 

# Convert to factor
Sample_metadata <- Sample_metadata %>%
  mutate(
    MSM_derived = factor(MSM_derived,
                         levels = c("msm", "msw", "wsm")),
    Sexual_orientation = factor(Sexual_orientation),
    Sexual_partner = factor(Sexual_partner,
                            levels = c(0, 1, 2),
                            labels = c("0", "1", "≥2"))
  )

# Plot the richness and diversity according to sex behaviour
pdf("2_DIVERSITY/Plots/Richness_Index_Sex_behaviour.pdf", width = 3.2, height = 3.4)
ggplot(Sample_metadata %>% filter(!is.na(MSM_derived)),
       aes(x = MSM_derived, y = richness, fill = MSM_derived)) +
  geom_jitter(color = "lightgrey", alpha = 0.4, size = 1.8, width = 0.2) +
  geom_boxplot(width = 0.6, color = "black", linewidth = 0.4, outlier.shape = NA, alpha = 0.85) +
  scale_fill_manual(values = c("msm" = "#7A1F3D", "msw" = "#6E8CA6", "wsm" = "#8A9A5B")) +
  labs(x = "Sex behaviour", y = "vOTU richness") +
  coord_cartesian(ylim = c(-20, 400)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.75),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.position = "none")
dev.off()


pdf('2_DIVERSITY/Plots/Shannon_Index_Sex_behaviour.pdf', width=2.7, height=3.2)
ggplot(Sample_metadata %>% filter(!is.na(MSM_derived)),
       aes(x = MSM_derived, y = shannon, fill = MSM_derived)) +
  geom_jitter(color = "lightgrey", alpha = 0.4, size = 1.8, width = 0.2) +
  geom_boxplot(width = 0.6, color = "black", linewidth = 0.4, outlier.shape = NA, alpha = 0.85) +
  scale_fill_manual(values = c("msm" = "#7A1F3D", "msw" = "#6E8CA6", "wsm" = "#8A9A5B")) +
  labs(x = "Sex behaviour", y = "Shannon index") +
  coord_cartesian(ylim = c(-0.5, 7)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.75),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.position = "none")
dev.off()

# Shannon diversity according to sexual orientation
pdf("2_DIVERSITY/Plots/Shannon_Index_Sex_orientation.pdf", width = 2.7, height = 3.2)
ggplot(Sample_metadata %>% filter(!is.na(Sexual_orientation)),
       aes(x = Sexual_orientation, y = shannon, fill = Sexual_orientation)) +
  geom_jitter(color = "lightgrey", alpha = 0.4, size = 1.8, width = 0.2) +
  geom_boxplot(width = 0.6, color = "black", linewidth = 0.4,
               outlier.shape = NA, alpha = 0.85) +
  scale_fill_manual(values = c("0" = "#54738C", "1" = "#9A7FA5", 
                               "2" = "#B56F72", "3" = "#8A9270")) +
  labs(x = "Sexual orientation", y = "Shannon index") +
  coord_cartesian(ylim = c(-0.5, 7)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.75),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.position = "none")
dev.off()

# Shannon diversity according to number of sexual partners
pdf("2_DIVERSITY/Plots/Shannon_Index_Sexual_partners.pdf", width = 2.7, height = 3.2)
ggplot(Sample_metadata %>% filter(!is.na(Sexual_partner)),
       aes(x = Sexual_partner, y = shannon, fill = Sexual_partner)) +
  geom_jitter(color = "lightgrey", alpha = 0.4, size = 1.8, width = 0.2) +
  geom_boxplot(width = 0.6, color = "black", linewidth = 0.4,
               outlier.shape = NA, alpha = 0.85) +
  scale_fill_manual(values = c("0" = "#5DA5A4", "1" = "#D4A64A", "≥2" = "#B65F3D")) +
  labs(x = "Number of sexual partners", y = "Shannon index") +
  coord_cartesian(ylim = c(-0.5, 7)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.75),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.position = "none")
dev.off()

##################################
# 1.2. Bacterial alpha diversity and sexual behaviour
##################################

# Estimate bacterial Shannon diversity
Sample_metadata$bacterial_shannon <- vegan::diversity(Abundance_table_bact,
  index = "shannon"
)

# Plot bacterial Shannon diversity according to sex behaviour
pdf("2_DIVERSITY/Plots/Bacterial_Shannon_Index_Sex_behaviour.pdf",
    width = 2.7, height = 3.2)
ggplot(Sample_metadata %>% filter(!is.na(MSM_derived)),
       aes(x = MSM_derived, y = bacterial_shannon, fill = MSM_derived)) +
  geom_jitter(color = "lightgrey", alpha = 0.4, size = 1.8, width = 0.2) +
  geom_boxplot(width = 0.6, color = "black", linewidth = 0.4,
               outlier.shape = NA, alpha = 0.85) +
  scale_fill_manual(values = c("msm" = "#7A1F3D","msw" = "#6E8CA6", "wsm" = "#8A9A5B")) +
  labs(x = "Sex behaviour", y = "Shannon index") +
  coord_cartesian(ylim = c(-0.5, 7)) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.75),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    legend.position = "none"
  )
dev.off()

##################################
# 1.3. Beta diversity analysis 
##################################

# Filter out vOTUs with prevalence < 5% (present in < 5% of individuals)
prevalence <- rowSums(Abundance_table > 0) / ncol(Abundance_table) * 100
Abundance_table_05prev <- Abundance_table[prevalence >= 5, ]
Abundance_table_05prev_nonzero <- Abundance_table_05prev[, colSums(Abundance_table_05prev) > 0 ]

# NMDS analysis
# Subset metadata to select phenotypes of interest
rownames(Sample_metadata) <- Sample_metadata$SampleID
NMDS_phenos <- Sample_metadata[Sample_metadata$SampleID %in% colnames(Abundance_table_05prev_nonzero),]
NMDS_phenos <- NMDS_phenos[, "MSM_derived", drop = FALSE]

# Run NMDS with vegan (metaMDS)
NMS <- metaMDS(t(Abundance_table_05prev_nonzero), distance = "bray", k=2)

# Perform fitting of metadata variables (how well variables correlate with NMDS axes)
en_all<- envfit(NMS, NMDS_phenos, permutations = 1000, na.rm = TRUE)

# Generate dataframe with NMDS site scores associated to phenotypes
data.scores_all <- as.data.frame(scores(NMS, "sites"))
data.scores_all <- merge(data.scores_all, NMDS_phenos, by='row.names')
row.names(data.scores_all) <- data.scores_all$Row.names
data.scores_all$Row.names <- NULL

# Calculate centroid and distances for each timepoint
data.scores_all <- data.scores_all %>%
  group_by(MSM_derived) %>%
  mutate(centroid_NMDS1 = mean(NMDS1),
         centroid_NMDS2 = mean(NMDS2),
         distance_from_centroid = sqrt((NMDS1 - centroid_NMDS1)^2 + (NMDS2 - centroid_NMDS2)^2))

# Use IQR method to identify and remove outliers within each timepoint
data.scores_all <- data.scores_all %>%
  group_by(MSM_derived) %>%
  mutate(Q1 = quantile(distance_from_centroid, 0.25),
         Q3 = quantile(distance_from_centroid, 0.75),
         IQR = Q3 - Q1,
         lower_bound = Q1 - 1.5 * IQR,
         upper_bound = Q3 + 1.5 * IQR) %>%
  filter(distance_from_centroid >= lower_bound & distance_from_centroid <= upper_bound) %>%
  ungroup()


# Generate NMDS plot
pdf('2_DIVERSITY/Plots/Viral_vOTUs_ALL_Bray_NMDS_5prev_Sex_behaviour.pdf', width=5.5, height=4)
NMDS_plot_all <- ggplot(data = data.scores_all[!is.na(data.scores_all$MSM_derived), ],
                        aes(x = NMDS1, y = NMDS2, color = MSM_derived)) + 
  geom_point(size = 1, alpha = 0.7) +  stat_ellipse(geom = "polygon", alpha = 0.0, 
               aes(group = MSM_derived, color = MSM_derived), 
               linetype = 2, linewidth = 1) +
  theme_bw() +
  scale_color_manual(values = c("msm" = "#7A1F3D", "msw" = "#6E8CA6", "wsm" = "#8A9A5B")) +
  theme(axis.text = element_text(size = 13),
        axis.title = element_text(size = 14), 
        panel.background = element_blank(), 
        panel.grid = element_blank(),  
        panel.border = element_rect(fill = NA, linewidth = 1.2, colour = "grey30"), 
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.title.x = element_text(margin = margin(t = 10)),
        legend.title = element_blank(), 
        legend.text = element_text(size = 14, colour = "grey30"),
        legend.position = "bottom") +
  guides(fill = "none") 
ggMarginal(NMDS_plot_all, type = "density", groupFill = TRUE, outlier.shape = NA)
dev.off()


##************************************************************************
# 2. Plots of Prevotella/Bacteroides ratio
#*************************************************************************

##############
# Viral ratio
##############

# Extract Prevotella and Bacteroides abundances
host_ratio <- data.frame(
  SampleID = colnames(Abundance_table_host),
  Prevotella = as.numeric(Abundance_table_host["Prevotella", ]),
  Bacteroides = as.numeric(Abundance_table_host["Bacteroides", ])
)

# Calculate log10(Prevotella/Bacteroides ratio)
pseudo <- 1e-6
host_ratio$log10_ratio <- log10((host_ratio$Prevotella + pseudo) /
                                  (host_ratio$Bacteroides + pseudo))

# Merge with metadata
host_ratio <- merge(host_ratio, Sample_metadata, by = "SampleID")

# Plot
pdf("2_DIVERSITY/Plots/Prevotella_Bacteroides_log10_ratio_Sex_behaviour.pdf",
    width = 2.9, height = 3.2)
ggplot(host_ratio %>% filter(!is.na(MSM_derived)), aes(x = MSM_derived, y = log10_ratio,
           fill = MSM_derived)) +
  geom_jitter(color = "lightgrey", alpha = 0.4, size = 1.8, width = 0.2) +
  geom_boxplot(width = 0.6, color = "black", linewidth = 0.4, outlier.shape = NA, alpha = 0.85) +
  scale_fill_manual(values = c("msm" = "#7A1F3D", "msw" = "#6E8CA6","wsm" = "#8A9A5B")) +
  labs(x = "Sex behaviour", y = expression(log[10]~"(Prevotella/Bacteroides ratio)")) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.75),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    legend.position = "none"
  )
dev.off()

##############
# Bacterial ratio
##############

# Calculate log10(Prevotella/Bacteroides ratio)
bacterial_ratio <- data.frame(
  SampleID = colnames(Abundance_table_bact_genus),
  Prevotella = as.numeric(Abundance_table_bact_genus["Prevotella", ]),
  Bacteroides = as.numeric(Abundance_table_bact_genus["Bacteroides", ])
) %>%
  mutate(
    log10_ratio = log10((Prevotella + 1e-6) / (Bacteroides + 1e-6))
  ) %>%
  left_join(Sample_metadata, by = "SampleID")

pdf("2_DIVERSITY/Plots/Bacterial_Prevotella_Bacteroides_log10_ratio_Sex_behaviour.pdf",
    width = 2.9, height = 3.2)
ggplot(filter(bacterial_ratio, !is.na(MSM_derived)),
       aes(MSM_derived, log10_ratio, fill = MSM_derived)) +
  geom_jitter(color = "lightgrey", alpha = 0.4, size = 1.8, width = 0.2) +
  geom_boxplot(width = 0.6, color = "black", linewidth = 0.4, outlier.shape = NA, alpha = 0.85) +
  scale_fill_manual(values = c(msm = "#7A1F3D", msw = "#6E8CA6", wsm = "#8A9A5B")) +
  labs(x = "Sex behaviour", y = expression(log[10]~"(Prevotella/Bacteroides ratio)")) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.75),
        axis.text = element_text(size = 14), axis.title = element_text(size = 16),
        legend.position = "none")
dev.off()


##************************************************************************
# 3. Association plots with clinical phenotypes
#*************************************************************************

######################################################
# 3.0. Prepare dataframe with viral groups
######################################################

# Convert immune cell counts to cells/µL
Sample_metadata$CD4_NADIR_cells <- Sample_metadata$CD4_NADIR * 1000
Sample_metadata$CD8_LATEST <- Sample_metadata$CD8_LATEST * 1000

# Master dataframe containing all host-group abundances + metadata
viral_df <- as.data.frame(t(Abundance_table_host)) %>%
  rownames_to_column("SampleID") %>%
  left_join(Sample_metadata, by = "SampleID")

# vOTU presence/absence dataframe
votu_df <- as.data.frame(t(Abundance_pa)) %>%
  rownames_to_column("SampleID") %>%
  left_join(Sample_metadata, by = "SampleID")


######################################################
# 3.1. HIV immune parameters
######################################################

# Prevotella abundance vs CD4 nadir
plot_df <- viral_df %>%
  filter(!is.na(Prevotella), !is.na(CD4_NADIR_cells))

pdf("Plots/Prevotella_CD4_NADIR.pdf", width = 3.2, height = 3.2)
ggplot(plot_df, aes(x = Prevotella, y = CD4_NADIR_cells)) +
  geom_point(size = 1.8, alpha = 0.6, colour = "#2A7F7B") +
  geom_smooth(method = "lm", colour = "black", linewidth = 0.8, se = TRUE) +
  labs(x = "Prevotella relative abundance", y = "CD4 nadir (cells/µL)") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16)
  )
dev.off()


# Alloprevotella presence vs CD4 nadir
plot_df <- viral_df %>%
  filter(!is.na(Alloprevotella), !is.na(CD4_NADIR_cells)) %>%
  mutate(
    presence = factor(
      ifelse(Alloprevotella > 0, "Present", "Absent"),
      levels = c("Absent", "Present")
    )
  )

pdf("Plots/Alloprevotella_CD4_nadir.pdf", width = 3, height = 3.4)
ggplot(plot_df, aes(x = presence, y = CD4_NADIR_cells, fill = presence)) +
  geom_jitter(color = "lightgrey", alpha = 0.4, size = 1.8, width = 0.15) +
  geom_boxplot(width = 0.6, color = "black", outlier.shape = NA, alpha = 0.85) +
  scale_fill_manual(values = c("Absent" = "#A6A6A6", "Present" = "#6D597A")) +
  labs(x = "Alloprevotella phages", y = "CD4 nadir (cells/µL)") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    legend.position = "none"
  )
dev.off()


# CAG-448 presence vs CD8
plot_df <- viral_df %>%
  filter(!is.na(`CAG-448`), !is.na(CD8_LATEST)) %>%
  mutate(
    presence = factor(
      ifelse(`CAG-448` > 0, "Present", "Absent"),
      levels = c("Absent", "Present")
    )
  )

pdf("Plots/CAG-448_CD8_LATEST.pdf", width = 2.7, height = 3.4)
ggplot(plot_df, aes(x = presence, y = CD8_LATEST, fill = presence)) +
  geom_jitter(color = "lightgrey", alpha = 0.4, size = 1.8, width = 0.15) +
  geom_boxplot(width = 0.6, color = "black", outlier.shape = NA, alpha = 0.85) +
  scale_fill_manual(values = c("Absent" = "#A6A6A6", "Present" = "#6D597A")) +
  labs(x = "CAG-448 phages", y = "CD8+ T cells (cells/µL)") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    legend.position = "none"
  )
dev.off()


# UBA4644 presence vs CD8
plot_df <- viral_df %>%
  filter(!is.na(UBA4644), !is.na(CD8_LATEST)) %>%
  mutate(
    presence = factor(
      ifelse(UBA4644 > 0, "Present", "Absent"),
      levels = c("Absent", "Present")
    )
  )

pdf("Plots/UBA4644_CD8_LATEST.pdf", width = 2.7, height = 3.4)
ggplot(plot_df, aes(x = presence, y = CD8_LATEST, fill = presence)) +
  geom_jitter(color = "lightgrey", alpha = 0.4, size = 1.8, width = 0.15) +
  geom_boxplot(width = 0.6, color = "black", outlier.shape = NA, alpha = 0.85) +
  scale_fill_manual(values = c("Absent" = "#A6A6A6", "Present" = "#6D597A")) +
  labs(x = "UBA4644 phages", y = "CD8+ T cells (cells/µL)") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    legend.position = "none"
  )
dev.off()

######################################################
# 3.2. Gut integrity markers
######################################################

# Prevotella abundance vs CD14
plot_df <- viral_df %>%
  filter(!is.na(Prevotella), !is.na(CD14_NPX))

pdf("Plots/Prevotella_CD14_NPX.pdf", width = 3.2, height = 3.2)
ggplot(plot_df, aes(x = Prevotella, y = CD14_NPX)) +
  geom_point(size = 1.8, alpha = 0.6, colour = "#2A7F7B") +
  geom_smooth(method = "lm", colour = "black", linewidth = 0.8, se = TRUE) +
  labs(x = "Prevotella-infecting phage abundance", y = "CD14") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16)
  )
dev.off()


# Agathobacter presence vs CD14
plot_df <- viral_df %>%
  filter(!is.na(Agathobacter), !is.na(CD14_NPX)) %>%
  mutate(
    presence = factor(
      ifelse(Agathobacter > 0, "Present", "Absent"),
      levels = c("Absent", "Present")
    )
  )

pdf("Plots/Agathobacter_CD14_NPX.pdf", width = 2.9, height = 3.4)
ggplot(plot_df, aes(x = presence, y = CD14_NPX, fill = presence)) +
  geom_jitter(color = "lightgrey", alpha = 0.4, size = 1.8, width = 0.15) +
  geom_boxplot(width = 0.6, color = "black", outlier.shape = NA, alpha = 0.85) +
  scale_fill_manual(values = c("Absent" = "#A6A6A6", "Present" = "#6D597A")) +
  labs(x = "Agathobacter-infecting phages", y = "CD14") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    legend.position = "none"
  )
dev.off()


# Acidaminococcus abundance vs I-FABP
plot_df <- viral_df %>%
  filter(!is.na(Acidaminococcus), !is.na(IFABP_NPX))

pdf("Plots/Acidaminococcus_IFABP_NPX.pdf", width = 3.2, height = 3.2)
ggplot(plot_df, aes(x = Acidaminococcus, y = IFABP_NPX)) +
  geom_point(size = 1.8, alpha = 0.6, colour = "#2A7F7B") +
  geom_smooth(method = "lm", colour = "black", linewidth = 0.8, se = TRUE) +
  labs(x = "Acidaminococcus-infecting phage abundance", y = "I-FABP") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16)
  )
dev.off()


# Spyradocola presence vs LBP
plot_df <- viral_df %>%
  filter(!is.na(Spyradocola), !is.na(LBP_NPX)) %>%
  mutate(
    presence = factor(
      ifelse(Spyradocola > 0, "Present", "Absent"),
      levels = c("Absent", "Present")
    )
  )

pdf("Plots/Spyradocola_LBP_NPX.pdf", width = 2.9, height = 3.4)
ggplot(plot_df, aes(x = presence, y = LBP_NPX, fill = presence)) +
  geom_jitter(color = "lightgrey", alpha = 0.4, size = 1.8, width = 0.15) +
  geom_boxplot(width = 0.6, color = "black", outlier.shape = NA, alpha = 0.85) +
  scale_fill_manual(values = c("Absent" = "#A6A6A6", "Present" = "#6D597A")) +
  labs(x = "Spyradocola-infecting phages", y = "LBP") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    legend.position = "none"
  )
dev.off()


######################################################
# 3.3. HIV infection and treatment history
######################################################

# CAG-245 presence vs HIV duration
plot_df <- viral_df %>%
  filter(!is.na(`CAG-245`), !is.na(HIV_DURATION)) %>%
  mutate(
    presence = factor(
      ifelse(`CAG-245` > 0, "Present", "Absent"),
      levels = c("Absent", "Present")
    )
  )

pdf("Plots/CAG-245_HIV_DURATION.pdf", width = 2.9, height = 3.4)
ggplot(plot_df, aes(x = presence, y = HIV_DURATION, fill = presence)) +
  geom_jitter(color = "lightgrey", alpha = 0.4, size = 1.8, width = 0.15) +
  geom_boxplot(width = 0.6, color = "black", outlier.shape = NA, alpha = 0.85) +
  scale_fill_manual(values = c("Absent" = "#A6A6A6", "Present" = "#6D597A")) +
  labs(x = "CAG-245-infecting phages", y = "HIV infection duration") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    legend.position = "none"
  )
dev.off()


# CAG-245 presence vs cART duration
plot_df <- viral_df %>%
  filter(!is.na(`CAG-245`), !is.na(CART_DURATION)) %>%
  mutate(
    presence = factor(
      ifelse(`CAG-245` > 0, "Present", "Absent"),
      levels = c("Absent", "Present")
    )
  )

pdf("Plots/CAG-245_cART_duration.pdf", width = 2.9, height = 3.4)
ggplot(plot_df, aes(x = presence, y = CART_DURATION, fill = presence)) +
  geom_jitter(color = "lightgrey", alpha = 0.4, size = 1.8, width = 0.15) +
  geom_boxplot(width = 0.6, color = "black", outlier.shape = NA, alpha = 0.85) +
  scale_fill_manual(values = c("Absent" = "#A6A6A6", "Present" = "#6D597A")) +
  labs(x = "CAG-245-infecting phages", y = "cART duration") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    legend.position = "none"
  )
dev.off()


# Blautia A presence vs HIV duration
plot_df <- viral_df %>%
  filter(!is.na(Blautia_A), !is.na(HIV_DURATION)) %>%
  mutate(
    presence = factor(
      ifelse(Blautia_A > 0, "Present", "Absent"),
      levels = c("Absent", "Present")
    )
  )

pdf("Plots/Blautia_A_HIV_DURATION.pdf", width = 2.9, height = 3.4)
ggplot(plot_df, aes(x = presence, y = HIV_DURATION, fill = presence)) +
  geom_jitter(color = "lightgrey", alpha = 0.4, size = 1.8, width = 0.15) +
  geom_boxplot(width = 0.6, color = "black", outlier.shape = NA, alpha = 0.85) +
  scale_fill_manual(values = c("Absent" = "#A6A6A6", "Present" = "#6D597A")) +
  labs(x = "Blautia A-infecting phages", y = "HIV infection duration") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    legend.position = "none"
  )
dev.off()


# HIV_43920 presence vs cART duration
plot_df <- votu_df %>%
  filter(!is.na(HIV_43920), !is.na(CART_DURATION))

pdf("Plots/HIV_43920_cART_duration.pdf", width = 2.9, height = 3.4)
ggplot(plot_df, aes(x = factor(HIV_43920), y = CART_DURATION,
                    fill = factor(HIV_43920))) +
  geom_jitter(color = "lightgrey", alpha = 0.4, size = 1.8, width = 0.15) +
  geom_boxplot(width = 0.6, color = "black", outlier.shape = NA, alpha = 0.85) +
  scale_fill_manual(values = c("0" = "#A6A6A6", "1" = "#6D597A"),
    labels = c("Absent", "Present")) +
  scale_x_discrete(labels = c("0" = "Absent", "1" = "Present")) +
  labs(x = "HIV_43920", y = "cART duration (years)") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    legend.position = "none"
  )
dev.off()


# HIV_43920 presence vs HIV duration
plot_df <- votu_df %>%
  filter(!is.na(HIV_43920), !is.na(HIV_DURATION))

pdf("Plots/HIV_43920_HIV_duration.pdf", width = 2.9, height = 3.4)
ggplot(plot_df, aes(x = factor(HIV_43920), y = HIV_DURATION,
                    fill = factor(HIV_43920))) +
  geom_jitter(color = "lightgrey", alpha = 0.4, size = 1.8, width = 0.15) +
  geom_boxplot(width = 0.6, color = "black", outlier.shape = NA, alpha = 0.85) +
  scale_fill_manual(
    values = c("0" = "#A6A6A6", "1" = "#6D597A"),
    labels = c("Absent", "Present")
  ) +
  scale_x_discrete(labels = c("0" = "Absent", "1" = "Present")) +
  labs(x = "HIV_43920", y = "HIV infection duration (years)") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    legend.position = "none"
  )
dev.off()


######################################################
# 3.4. Five-year clinical outcomes
######################################################

# Hungatella presence vs mortality
plot_df <- viral_df %>%
  filter(!is.na(Hungatella), !is.na(DEATH_SHM)) %>%
  mutate(
    presence = factor(
      ifelse(Hungatella > 0, "Present", "Absent"),
      levels = c("Absent", "Present")
    ),
    outcome = factor(
      DEATH_SHM,
      levels = c(0, 1),
      labels = c("No", "Yes")
    )
  )

prevalence <- plot_df %>%
  group_by(outcome) %>%
  summarise(prevalence = mean(presence == "Present") * 100)


pdf("Plots/Hungatella_DEATH_SHM.pdf", width = 2.9, height = 3.4)
ggplot(plot_df, aes(x = outcome, fill = presence)) +
  geom_bar(position = "fill", width = 0.65, color = "black") +
  geom_text(data = prevalence,
    aes(x = outcome, y = 1.05, label = paste0(round(prevalence, 1), "%")),
    inherit.aes = FALSE, size = 4) +
  scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.15))) +
  scale_fill_manual(values = c("Absent" = "#D9D9D9", "Present" = "#6D597A")) +
  labs(x = "All-cause mortality", y = "Hungatella prevalence", fill = NULL) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    legend.position = "bottom"
  )
dev.off()


# Anaerotardibacter presence vs incident malignancy
plot_df <- viral_df %>%
  filter(!is.na(Anaerotardibacter), !is.na(incident_MALIG_SHM)) %>%
  mutate(
    presence = factor(
      ifelse(Anaerotardibacter > 0, "Present", "Absent"),
      levels = c("Absent", "Present")
    ),
    outcome = factor(
      incident_MALIG_SHM,
      levels = c(0, 1),
      labels = c("No", "Yes")
    )
  )

prevalence <- plot_df %>%
  group_by(outcome) %>%
  summarise(prevalence = mean(presence == "Present") * 100)


pdf("Plots/Anaerotardibacter_incident_MALIG_SHM.pdf",
    width = 2.9, height = 3.4)
ggplot(plot_df, aes(x = outcome, fill = presence)) +
  geom_bar(position = "fill", width = 0.65, color = "black") +
  geom_text(data = prevalence,
    aes(x = outcome, y = 1.05, label = paste0(round(prevalence, 1), "%")),
    inherit.aes = FALSE, size = 4) +
  scale_y_continuous( labels = scales::percent, expand = expansion(mult = c(0, 0.15))) +
  scale_fill_manual(values = c("Absent" = "#D9D9D9", "Present" = "#6D597A")) +
  labs(x = "Incident malignancy", y = "Anaerotardibacter prevalence", fill = NULL) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    legend.position = "bottom"
  )
dev.off()


######################################################
# 3.5. HIV controller status
######################################################

# HIV_59378 presence vs EC_ART
plot_df <- votu_df %>%
  filter(!is.na(HIV_59378), !is.na(EC_ART))

prevalence <- plot_df %>%
  group_by(EC_ART) %>%
  summarise(prevalence = mean(HIV_59378 == 1) * 100)

pdf("Plots/HIV_59378_EC_proportion.pdf", width = 3.7, height = 3.4)
ggplot(plot_df, aes(EC_ART, fill = factor(HIV_59378))) +
  geom_bar(position = "fill", width = .65, color = "black") +
  geom_text(data = prevalence, aes(EC_ART, 1.05,
                                   label = paste0(round(prevalence, 1), "%")), inherit.aes = FALSE) +
  scale_y_continuous(labels = scales::percent,
                     expand = expansion(mult = c(0, .15))) +
  scale_fill_manual(values = c("0" = "#D9D9D9", "1" = "#6D597A"),
                    labels = c("Absent", "Present")) +
  labs(x = "HIV controller status", y = "HIV_59378 prevalence", fill = NULL) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.position = "bottom")
dev.off()

######################################################
# 3.6. COVID-19
######################################################

covid_df <- viral_df %>%
  filter(PANDEMIC_BEFOREAFTER == 1)

# Bacteroides abundance vs COVID group
plot_df <- covid_df %>%
  filter(!is.na(Bacteroides), !is.na(COVID_GROUP))

pdf("Plots/Bacteroides_COVID_GROUP.pdf", width = 3.7, height = 3.4)
ggplot(plot_df, aes(COVID_GROUP, Bacteroides, fill = COVID_GROUP)) +
  geom_jitter(color = "lightgrey", alpha = .4, size = 1.8, width = .2) +
  geom_boxplot(width = .6, color = "black", outlier.shape = NA, alpha = .85) +
  scale_fill_manual(values = c("#355070", "#6D597A", "#B56576", "#EAAC8B")) +
  labs(x = "COVID group", y = "Bacteroides relative abundance") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.position = "none")
dev.off()


# Mediterraneibacter prevalence vs COVID group
plot_df <- covid_df %>%
  filter(!is.na(Mediterraneibacter), !is.na(COVID_GROUP)) %>%
  mutate(presence = factor(ifelse(Mediterraneibacter > 0, "Present", "Absent"),
                           levels = c("Absent", "Present")))

prevalence <- plot_df %>%
  group_by(COVID_GROUP) %>%
  summarise(prevalence = mean(presence == "Present") * 100)

pdf("Plots/Mediterraneibacter_COVID_GROUP.pdf", width = 3.7, height = 3.4)
ggplot(plot_df, aes(COVID_GROUP, fill = presence)) +
  geom_bar(position = "fill", width = .65, color = "black") +
  geom_text(data = prevalence, aes(COVID_GROUP, 1.05,
                                   label = paste0(round(prevalence, 1), "%")), inherit.aes = FALSE) +
  scale_y_continuous(labels = scales::percent,
                     expand = expansion(mult = c(0, .15))) +
  scale_fill_manual(values = c("Absent" = "#D9D9D9", "Present" = "#6D597A")) +
  labs(x = "COVID group", y = "Mediterraneibacter prevalence", fill = NULL) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.position = "bottom")
dev.off()


##************************************************************************
# 4. Phage genome plot
#*************************************************************************

# After annotating associated genomes (Phold), we load results
annotation_vOTUs <- read.delim("GENOME_ANNOTATIONS/ANNOTATIONS/merged_per_cds_predictions.tsv", 
                               check.names = F)

# Keep only viruses of interest
phages <- annotation_vOTUs[annotation_vOTUs$contig_id %in% c("HIV_43920", "HIV_59378"),
                              c("contig_id", "cds_id", "start", "end", "strand",
                                "function", "product")]

# 2. Visualize vOTUs

# Standarize column names
votu_annot <- phages  %>%
  transmute(
    Genome    = contig_id,
    GeneID    = cds_id,
    Start     = as.numeric(start),
    Stop      = as.numeric(end),
    Direction = ifelse(strand == "+", TRUE, FALSE),
    Annot     = `function`,
    Product   = product
  )

# Estimate genome lengths
genome_lengths <- votu_annot %>%
  group_by(Genome) %>%
  summarise(GenomeLength = max(Stop, na.rm = TRUE)) %>%
  arrange(GenomeLength)

genome_lengths <- genome_lengths %>%
  arrange(desc(GenomeLength))

votu_annot <- votu_annot %>%
  mutate(
    Genome = factor(Genome, levels = genome_lengths$Genome)
  )

# Define bin annotations
categories <- c("connector", "DNA, RNA and nucleotide metabolism",
                "head and packaging", "integration and excision",
                "lysis", "moron, auxiliary metabolic gene and host takeover",
                "tail", "transcription regulation", "other")

votu_annot <- votu_annot %>%
  mutate(
    AnnotCat = ifelse(
      Annot %in% categories,
      Annot,
      "unknown"
    )
  )

votu_annot$AnnotCat <- factor(
  votu_annot$AnnotCat,
  levels = names(palette)
)

# Set a color palette
palette <- c(
  "connector" = "#2a9d8f",       
  "DNA, RNA and nucleotide metabolism" = "#264653", 
  "head and packaging" = "#e76f51", 
  "integration and excision" = "#f4a261", 
  "lysis" = "#a8dadc",         
  "moron, auxiliary metabolic gene and host takeover" = "#e9c46a", 
  "tail" = "#6a994e",           
  "transcription regulation" = "#457b9d", 
  "other" = "#bfc5c8",          
  "unknown" = "#f0efeb"         
)

# Generate the genome map plot
ggplot(votu_annot, aes(xmin = Start, xmax = Stop, y = Genome, forward = Direction,
                       color = AnnotCat, fill = AnnotCat)) +
  geom_gene_arrow() +
  scale_color_manual(values = palette) + 
  scale_x_continuous(labels = label_comma()) + 
  scale_fill_manual(values = palette) + 
  theme_bw() + 
  theme(strip.background = element_blank(),
        strip.text.x = element_blank(), 
        axis.title.y = element_blank(), 
        panel.border = element_blank(), 
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(), 
        axis.ticks.y = element_blank(), 
        panel.grid.major.y= element_blank(), 
        legend.text = element_text(size=10)
  )

ggsave("GENOME_ANNOTATIONS/Viral_genome_map.pdf", dpi = 300, w = 8, h = 4)


##************************************************************************
 # Write results
##************************************************************************
write.table(Sample_metadata, file = "METADATA/Sample_metadata_excluded_11052026.txt",
            quote = F, row.names = T, sep="\t")

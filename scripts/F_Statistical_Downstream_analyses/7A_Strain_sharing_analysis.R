################################################################################
##### 2000HIV: Strain sharing analysis - Data processing
### Author(s):Asier Fernández-Pato
### Last updated: 14th September, 2026
################################################################################

#****************
# Load modules
#****************
library(data.table) 
library(dplyr)       
library(readxl)  


# Set working directory 
setwd("~/Desktop/PhD/Projects/HIV 2000/")

# Read abundance table and metadata tables
Abundance_table <- read.delim("ABUNDANCE TABLE/HIV_Abundance_Table_01052025.txt")
Sample_metadata <- read.delim("METADATA/Sample_metadata_excluded_11052026.txt")
Viral_metadata <- read.delim("METADATA/HIV_Viral_Metadata_06052025.txt")
participants_removed <-  read.delim("METADATA/participants_removed.txt", col.names = "ID")

# Remove excluded participants
Abundance_table <- Abundance_table [, ! colnames(Abundance_table) %in% participants_removed$ID]

##************************************************************************
# 1. Process in inStrain results: Add metadata 
#*************************************************************************

# After running inStrain, we load the results
inStrain_results <- fread("8_STRAIN_TRANSMISSION/all_comparisons.tsv")

# Filter those comparisons with a minimum of 75% of the genome compared
inStrain_results <- inStrain_results[inStrain_results$percent_genome_compared >=0.75,]

# Format the results table
inStrain_results$name1 <- gsub("\\.sorted\\.bam$", "", inStrain_results$name1)
inStrain_results$name2 <- gsub("\\.sorted\\.bam$", "", inStrain_results$name2)
inStrain_results <- inStrain_results[, c("scaffold","name1","name2","popANI")]
colnames(inStrain_results) <- c("Virus_ID", "Sample1", "Sample2", "popANI")

# Change the sample IDs
Name_matching <-  readRDS("METADATA/name_trans.rds")
lookup <- Name_matching %>%
  select(sample_name, sample_name2)

inStrain_results$Sample1 <- lookup$sample_name[
  match(inStrain_results$Sample1, lookup$sample_name2)
]

inStrain_results$Sample2 <- lookup$sample_name[
  match(inStrain_results$Sample2, lookup$sample_name2)
]

# Add sample metadata
inStrain_results <- inStrain_results %>%
  mutate(
    Sample_ID1 = Sample_metadata$SampleID[match(Sample1, Sample_metadata$SampleID)],
    Sample_ID2 = Sample_metadata$SampleID[match(Sample2, Sample_metadata$SampleID)],
    Sex1 = Sample_metadata$SEX_BIRTH[match(Sample1, Sample_metadata$SampleID)],
    Sex2 = Sample_metadata$SEX_BIRTH[match(Sample2, Sample_metadata$SampleID)],
    Sex_behaviour1 = Sample_metadata$MSM_derived[match(Sample1, Sample_metadata$SampleID)],
    Sex_behaviour2 = Sample_metadata$MSM_derived[match(Sample2, Sample_metadata$SampleID)],
    Sex_orientation1 = Sample_metadata$Sexual_orientation[match(Sample1, Sample_metadata$SampleID)],
    Sex_orientation2 = Sample_metadata$Sexual_orientation[match(Sample2, Sample_metadata$SampleID)],
    Vaginal_sex1 = Sample_metadata$Vaginal_sexual_contact[match(Sample1, Sample_metadata$SampleID)],
    Vaginal_sex2 = Sample_metadata$Vaginal_sexual_contact[match(Sample2, Sample_metadata$SampleID)],
    Anal_sex1 = Sample_metadata$Anal_sexual_contact[match(Sample1, Sample_metadata$SampleID)],
    Anal_sex2 = Sample_metadata$Anal_sexual_contact[match(Sample2, Sample_metadata$SampleID)],
    Ethnicity1 = Sample_metadata$ETHNICITY[match(Sample1, Sample_metadata$SampleID)],
    Ethnicity2 = Sample_metadata$ETHNICITY[match(Sample2, Sample_metadata$SampleID)], 
    Center1 = Sample_metadata$CENTER[match(Sample1, Sample_metadata$SampleID)],
    Center2 = Sample_metadata$CENTER[match(Sample2, Sample_metadata$SampleID)],
    Sexual_partner1 = Sample_metadata$Sexual_partner[match(Sample1, Sample_metadata$SampleID)],
    Sexual_partner2 = Sample_metadata$Sexual_partner[match(Sample2, Sample_metadata$SampleID)],
  )

# Generate paired variables for comparisons
inStrain_results <- inStrain_results %>%
  mutate(
    Pair_ID = paste(pmin(Sample1, Sample2), pmax(Sample1, Sample2), sep = "_"),
    
    Sex_behaviour_pair = paste(pmin(Sex_behaviour1, Sex_behaviour2),
                               pmax(Sex_behaviour1, Sex_behaviour2), sep = "-"),
    
    Sex_orientation_pair = paste(pmin(Sex_orientation1, Sex_orientation2),
                                 pmax(Sex_orientation1, Sex_orientation2), sep = "-"),
    
    Sex_pair = paste(pmin(Sex1, Sex2), pmax(Sex1, Sex2), sep = "-"),
    
    Vaginal_sex_pair = paste(pmin(Vaginal_sex1, Vaginal_sex2),
                             pmax(Vaginal_sex1, Vaginal_sex2), sep = "-"),
    
    Anal_sex_pair = paste(pmin(Anal_sex1, Anal_sex2),
                          pmax(Anal_sex1, Anal_sex2), sep = "-"),
    
    Ethnicity_pair = paste(pmin(Ethnicity1, Ethnicity2),
                           pmax(Ethnicity1, Ethnicity2), sep = "-"),
    
    Center_pair = paste(pmin(Center1, Center2),
                        pmax(Center1, Center2), sep = "-"),
    
    Sexual_partner_pair = paste(pmin(Sexual_partner1, Sexual_partner2),
                                pmax(Sexual_partner1, Sexual_partner2), sep = "-")
  )

# Add extra variables for downstream analysis and comparisons
inStrain_results <- inStrain_results %>%
  mutate(
    same_sex_behaviour = Sex_behaviour1 == Sex_behaviour2,
    both_msm = Sex_behaviour1 == "msm" & Sex_behaviour2 == "msm",
    both_msw = Sex_behaviour1 == "msw" & Sex_behaviour2 == "msw",
    
    same_sex_orientation = Sex_orientation1 == Sex_orientation2,
    both_homo = Sex_orientation1 == 2 & Sex_orientation2 == 2,
    both_hetero = Sex_orientation1 == 0 & Sex_orientation2 == 0,
    
    same_vaginal_sex = Vaginal_sex1 == Vaginal_sex2,
    both_vaginal_sex = Vaginal_sex1 == 1 & Vaginal_sex2 == 1,
    
    same_anal_sex = Anal_sex1 == Anal_sex2,
    both_anal_sex = Anal_sex1 == 1 & Anal_sex2 == 1,
    
    same_ethnicity = Ethnicity1 == Ethnicity2,
    same_center = Center1 == Center2,
    same_n_partners = Sexual_partner1 == Sexual_partner2,
    more_than1_partner = Sexual_partner1 == "2" & Sexual_partner2 == "2"
  )

# Add information about potential transmission (based on inStrain viral cut-off)
inStrain_results$Sharing <- ifelse(inStrain_results$popANI > 0.99999 , "Yes", "No")

# Select vOTUs with strain-sharing
shared_vOTUs <- unique(inStrain_results$Virus_ID[inStrain_results$Sharing == "Yes"])

# Count number of total pairwise comparisons, number of vOTUs involved and number of MSM comparisons
nrow(inStrain_results) #321158
length(unique(inStrain_results$Virus_ID)) #5090
table(inStrain_results$Sex_behaviour_pair) #208421 MSM comparisons

# Remove excluded participants
inStrain_results <- inStrain_results[
  !(inStrain_results$Sample1 %in% participants_removed$ID |
      inStrain_results$Sample2 %in% participants_removed$ID),
]

# Remove HIV_91378 (potential contaminant based on strain sharing results)
inStrain_results <- inStrain_results[inStrain_results$Virus_ID != "HIV_91378", ]
  
# Count number of comparisons and vOTUs
nrow(inStrain_results) #273542
length(unique(inStrain_results$Virus_ID)) #4813
nrow(inStrain_results[inStrain_results$both_msm == TRUE,]) #200585
nrow(inStrain_results[inStrain_results$both_msw == TRUE,]) #2896

# Count potential sharing events
nrow(inStrain_results[inStrain_results$Sharing == "Yes",]) #16660

# Generate an inStrain results table separately for men only 
inStrain_results_men <- inStrain_results[inStrain_results$Sex_pair == "0-0",]


##************************************************************************
# Write results
##************************************************************************

write.table(inStrain_results, file = "8_STRAIN_TRANSMISSION/inStrain_results.txt",
            sep = "\t", row.names = F, quote = FALSE)
write.table(inStrain_results_men, file = "8_STRAIN_TRANSMISSION/inStrain_results_men.txt",
            sep = "\t", row.names = F, quote = FALSE)

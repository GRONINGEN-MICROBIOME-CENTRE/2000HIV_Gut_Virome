# Load R packages
library(stats)
library(Biostrings)
library(data.table)

# Get FASTA file of viral genomes from command line arguments
args <- commandArgs(trailingOnly = TRUE)
viral_sequences <- readDNAStringSet(args[1]) # FASTA file with the viral sequences 

# Get sequence IDs
viral_DB_origin <- data.frame(names(viral_sequences))

# Add DB information
colnames(viral_DB_origin)[1] <- "viral_seq"
viral_DB_origin$DB <- NA
viral_DB_origin[,"DB"] [grep("NODE", viral_DB_origin$viral_seq)] <- "HIV"  
viral_DB_origin[,"DB"] [grep("UHGV", viral_DB_origin$viral_seq)] <- "UHGV" 
viral_DB_origin[,"DB"] [grep("^NC.*\\.|^AC.*\\.", viral_DB_origin$viral_seq)] <- "RefSeq" #5,199

# Add new simplified viral IDs
viral_DB_origin$new_viral_ID <- ifelse(				    
    viral_DB_origin$DB == "UHGV",
    # Extract vOTU name from UHGV DB
    paste("UHGV", sub("UHGV-([0-9]+) .*", "\\1", viral_DB_origin$viral_seq), sep = "_"),
    # For remaining (HIV and RefSeq) sequences, generate a new name
    ave(viral_DB_origin$DB, 
        viral_DB_origin$DB, 
        FUN = function(x) paste0(x, "_", seq_along(x)))
)

# Replace the IDs in the FASTA file with the new simplified IDs
names(viral_sequences) <- viral_DB_origin$new_viral_ID

# Save final table with the X sequences used as input for STEP5 dereplication and their DB of origin (including 5 NEG-CONTROLS)
write.table(viral_DB_origin,"Dereplication_input_sequences_DB_origin.txt", sep = "\t", 
            row.names = F, col.names = F, quote = FALSE)

# Save the updated FASTA file
writeXStringSet(viral_sequences, "Dereplication_input_sequences_renamed.fa")

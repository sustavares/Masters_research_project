rm(list = ls())
graphics.off()

library(data.table)

# Import 

dt_cadd <- fread("../data/pathogenicity_scores/CADD_v1.7_hcm_dcm_cardiacG2P_all_transcript_CDS_splice.tsv.gz")

### Format ==========

# Subset columns and change column names 
colnames(dt_cadd) <- c("chromosome","pos","ref","alt", "CADD_raw","CADD_PHRED")

dt_cadd$chromosome <- paste0("chr",dt_cadd$chromosome)

# Convert pos to character 
dt_cadd$pos <- as.character(dt_cadd$pos)

dt_cadd <- dt_cadd[complete.cases(dt_cadd)]

# Export 
fwrite(dt_cadd, "../data/processed/CADD_v1.7_hcm_dcm_cardiacG2P_all_transcript_CDS_spliceprocessed.tsv.gz")

# Import processed CADD
dt_cadd <- fread("../data/processed/CADD_v1.7_hcm_dcm_cardiacG2P_all_transcript_CDS_spliceprocessed.tsv.gz")

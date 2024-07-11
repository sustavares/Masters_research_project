rm(list = ls())
graphics.off()

library(data.table)

install.packages("survminer")
install.packages("survival")

library(survival)
library(ggplot2)
library(survminer)

# Import
dt <- fread("../data/analysis_table.csv.gz")

survivorship_analysis_table <- fread("data/survivorship_analysis_table.csv")

survivorship_analysis_table_merged <- fread("data/survivorship_analysis_table_merged.czv.gz")

# Create VariantID column in dt
dt[, VariantID := paste(chromosome, pos, ref, alt, sep = "-")]

dt[, VariantID := sub("chr", "", VariantID)]

# Merge
survivorship_analysis_table_merged <- dt[survivorship_analysis_table, on = "VariantID"]

# Export

fwrite(survivorship_analysis_table_merged, "data/survivorship_analysis_table_merged.czv.gz")


# Clinical significance 
survivorship_analysis_table_merged$ClinVar_significance_binary <- NA
survivorship_analysis_table_merged$ClinVar_significance_binary[survivorship_analysis_table_merged$ClinVar_significance %in% c("Likely pathogenic",
                                                              "Pathogenic",
                                                              "Pathogenic/Likely pathogenic")] <- 1
survivorship_analysis_table_merged$ClinVar_significance_binary[survivorship_analysis_table_merged$ClinVar_significance %in% c("Likely benign",
                                                              "Benign",
                                                              "Benign/Likely benign")] <- 0

# DCM genes with evidence of missense pathogenicity
dcm_genes <- c("BAG3", "DES", "DSP", "LMNA", "MYH7", "PLN", "RBM20", "SCN5A", "TNNC1", "TNNT2", "TTN")

# HCM genes with evidence of missense pathogenicity
hcm_genes <- c("ACTC1", "MYBPC3", "MYH7", "MYL2", "MYL3", "PLN", "TNNI3", "TNNT2", "TPM1")

# HCM and DCM genes 
hcm_dcm_genes <- dcm_genes[dcm_genes %in% hcm_genes]


# Set score binary predictions based on recommended thresholds

survivorship_analysis_table_merged$AlphaMissense_binary <- NA
survivorship_analysis_table_merged$AlphaMissense_binary[survivorship_analysis_table_merged$AlphaMissense_class == "likely_pathogenic"] <- 1
survivorship_analysis_table_merged$AlphaMissense_binary[survivorship_analysis_table_merged$AlphaMissense_class == "likely_benign"] <- 0

survivorship_analysis_table_merged$CardioBoost_binary <- NA
survivorship_analysis_table_merged$CardioBoost_binary[survivorship_analysis_table_merged$CardioBoost_pathogenicity >= 0.1] <- 1
survivorship_analysis_table_merged$CardioBoost_binary[survivorship_analysis_table_merged$CardioBoost_pathogenicity <= 0.9] <- 0

survivorship_analysis_table_merged$REVEL_binary <- NA
survivorship_analysis_table_merged$REVEL_binary[survivorship_analysis_table_merged$REVEL >= 0.75] <- 1
survivorship_analysis_table_merged$REVEL_binary[survivorship_analysis_table_merged$REVEL <= 0.25] <- 0

survivorship_analysis_table_merged$CADD_binary <- NA
survivorship_analysis_table_merged$CADD_binary[survivorship_analysis_table_merged$CADD_PHRED >= 20] <- 1
survivorship_analysis_table_merged$CADD_binary[survivorship_analysis_table_merged$CADD_PHRED <= 10] <- 0

# Add a new column 'IsPathogenic' based on the binary columns
survivorship_analysis_table_merged[, IsPathogenic := ifelse(ClinVar_significance_binary == 1 | AlphaMissense_binary == 1 | CardioBoost_binary == 1 | REVEL_binary == 1 | CADD_binary == 1, "Yes", "No")]





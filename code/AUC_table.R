rm(list = ls())
graphics.off()

library(data.table)

# Import

AUC_table <- fread("../data/auc_table.csv.gz")

ClinVar_variants_table <- fread("../data/ClinVar_variants_table.csv.gz")

AUC_table_genes <- fread("../data/AUC_table_genes.csv.gz")

CV_variant_proportions <- fread("../data/cv_variant_proportions.csv.gz")

AUC_table_genes_merged <- fread("../data/auc_table_genes_merged.csv.gz")


# Create table of AUC values for AlphaMissense, CardioBoost and Revel

df <- data.frame(Genes = c("HCM genes", "DCM genes", "Total genes"), AlphaMissense_ROC_AUC = c(0.911352, 0.9233329, 0.9100733), AlphaMissense_PR_AUC = c(0.047468721, 0.4044371, 0.3941483), CardioBoost_ROC_AUC = c(0.9355648, 0.9255994, 0.9432624), CardioBoost_PR_AUC = c(0.03496824, 0.02510042, 0.03239188), REVEL_ROC_AUC = c(0.9103529, 0.9705611, 0.95422), REVEL_PR_AUC = c(0.04770143, 0.3982839, 0.3872899))

print(df)


AUC_table <- data.table(df)

# Export

fwrite(AUC_table, "../data/auc_table.csv.gz")


# Genes Table


# Create table of no. pathogenic and benign ClinVar variants for each gene 

df <- data.frame(Gene = c("ACTC1", "MYBPC3", "MYH7", "MYL2", "MYL3", "PLN", "TNNI3", "TNNT2", "TPM1", "BAG3", "DES", "DSP", "FLNC", "LMNA", "SCN5A", "TNNC1", "RBM20", "TTN", "Total"), Number_of_ClinVar_Pathogenic_Variants = c(13, 32, 307, 9, 3, 1, 43, 48, 35, 6, 31, 15, 15, 168, 88, 9, 11, 21, 855), Number_of_ClinVar_Benign_Variants = c(0, 27, 15, 0, 0, 0, 3, 3, 0, 17, 2, 154, 12, 1, 21, 0, 100, 818, 1283))

print(df)

ClinVar_variants_table <- data.table(df)

# Export 

fwrite(ClinVar_variants_table, "../data/ClinVar_variants_table.csv.gz")


# Create table for AUC values for individual genes for AlphaMissense, CardioBoost and REVEL

df <- data.frame(Gene = c("ACTC1", "MYBPC3", "MYH7", "MYL2", "MYL3", "PLN", "TNNI3", "TNNT2", "TPM1", "BAG3", "DES", "DSP", "FLNC", "LMNA", "SCN5A", "TNNC1", "RBM20", "TTN", "Total"), AlphaMissense_ROC_AUC = c(NA, 0.8229167, 0.8598263, NA, NA, NA, 0.9844961, 0.8958333, NA, 1, 1, 0.9350649, 0.8868852, 0.6666667, 0.9821429, NA, 0.9718182,0.8704254, 0.9100733), AlphaMissense_PR_AUC = c(NA, 0.3042065, 0.02528909, NA, NA, NA, 0.03352486, 0.03153205, NA, 0.5257407, 0.03192218, 0.7839659, 0.7554817, 0.004437889, 0.1033567, NA, 0.7522313, 0.9298989, 0.3941483), CardioBoost_ROC_AUC = c(NA, 0.9394531, 0.9164728, NA, NA, NA, 0.9883721, 0.8235294, NA, NA, 0.9, NA, NA, NA, 0.9967949, NA, NA, NA, 0.9432624), CardioBoost_PR_AUC = c(NA, 0.1936676, 0.02265677, NA, NA, NA, 0.02269315, 0.0340998, NA, NA, 0.01807361, NA, NA, NA, 0.115088, NA, NA, NA, 0.03239188), REVEL_ROC_AUC = c(NA, 0.9189815, 0.843325, NA, NA, NA, 0.8294574, 0.9791667, NA, 1, 1, 0.8329004, 0.8382514, 0.3095238, 0.9972944, NA, 0.9666048, 0.8177868, 0.95422), REVEL_PR_AUC = c(NA, 0.2827168, 0.02797986, NA, NA, NA, 0.03713616, 0.03061868, NA, 0.5257407, 0.03192218, 0.8256591, 0.7777076, 0.009493862, 0.1032248, NA, 0.7581775, 0.9333879, 0.3941483))

print(df)

AUC_table_genes <- data.table(df)

# Export 
fwrite(AUC_table_genes, "../data/auc_table_genes.csv.gz")

# Merge 
AUC_table_genes_merged <- unique(merge(AUC_table_genes, ClinVar_variants_table, by=c("Gene"), all.x=TRUE))

AUC_table_genes_merged <- unique(merge(AUC_table_genes_merged, CV_variant_proportions, by=c("Gene"), all.x=TRUE))

# Export
fwrite(AUC_table_genes_merged, "../data/auc_table_genes_merged.csv.gz")

# Proportions of pathogenic and benign ClinVar variants with values for each score

df <- data.frame(Gene = c("ACTC1", "MYBPC3", "MYH7", "MYL2", "MYL3", "PLN", "TNNI3", "TNNT2", "TPM1", "BAG3", "DES", "DSP", "FLNC", "LMNA", "SCN5A", "TNNC1", "RBM20", "TTN", "Total"), Proportion_Pathogenic_ClinVar_Variants_AlphaMissense = c(1, 1, 1, 0.888, 1, 1, 1, 1, 1, 1, 0.968, 1, 1, 1, 1, 1, 1, 1, 0.998), Proportion_Benign_ClinVar_Variants_AlphaMissense = c(0, 1, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 0.940, 0.962), Proportion_Pathogenic_ClinVar_Variants_CardioBoost = c(1, 1, 1, 0.888, 1, 1, 1, 0.354, 1, 0, 0.968, 0, 0, 1, 0.545, 0, 0, 0, 0.825), Proportion_Benign_ClinVar_Variants_CardioBoost = c(0, 0.593, 0.933, 0, 0, 0, 0.667, 0.333, 0, 0, 0.5, 0, 0, 0, 0.619, 0, 0, 0, 0.037), Proportion_Pathogenic_ClinVar_Variants_REVEL = c(1, 1, 1, 0.889, 1, 1, 1, 1, 1, 1, 0.968, 1, 1, 1, 1, 1, 1, 1, 0.998), Proportion_Benign_ClinVar_Variants_REVEL = c(0, 1, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0.98, 0.95, 0.967))

print(df)

CV_variant_proportions <- data.table(df)

# Export 
fwrite(CV_variant_proportions, "../data/cv_variant_proportions.csv.gz")





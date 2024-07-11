rm(list = ls())
graphics.off()

library(data.table)

library(PRROC)
library(pROC)

# Import ========

dt <- fread("../data/analysis_table.csv.gz")


# Generating pathogenic scores and benign scores for AlphaMissense

scores_pathogenic <- dt$AlphaMissense_pathogenicity[dt$ClinVar_ClinicalSignificance %in% c("Likely pathogenic", "Pathogenic", "Pathogenic/Likely pathogenic")]
scores_benign <- dt$AlphaMissense_pathogenicity[dt$ClinVar_ClinicalSignificance %in% c("Likely benign", "Benign", "Benign/Likely benign")]

scores_benign <- scores_benign[!is.na(scores_benign)]
scores_pathogenic <- scores_pathogenic[!is.na(scores_pathogenic)]

# Calculate PR-AUC for AlphaMissense - total genes

pr <- pr.curve(scores.class0 = scores_benign, scores.class1 = scores_pathogenic, curve = TRUE)
am_aupr <- pr$auc.integral


# Generating pathogenic scores and benign scores for CardioBoost 
scores_pathogenic <- dt$CardioBoost_pathogenicity[dt$ClinVar_ClinicalSignificance %in% c("Likely pathogenic", "Pathogenic", "Pathogenic/Likely pathogenic")]

scores_benign <- dt$CardioBoost_pathogenicity[dt$ClinVar_ClinicalSignificance %in% c("Likely benign", "Benign", "Benign/Likely benign")]

scores_benign <- scores_benign[!is.na(scores_benign)]
scores_pathogenic <- scores_pathogenic[!is.na(scores_pathogenic)]

# Calculate PR-AUC for CardioBoost - total genes
pr <- pr.curve(scores.class0 = scores_benign, scores.class1 = scores_pathogenic, curve = TRUE)
cb_aupr <- pr$auc.integral


# Generating pathogenic scores and benign scores for REVEL
scores_pathogenic <- dt$REVEL[dt$ClinVar_ClinicalSignificance %in% c("Likely pathogenic", "Pathogenic", "Pathogenic/Likely pathogenic")]

scores_benign <- dt$REVEL[dt$ClinVar_ClinicalSignificance %in% c("Likely benign", "Benign", "Benign/Likely benign")]

scores_benign <- scores_benign[!is.na(scores_benign)]
scores_pathogenic <- scores_pathogenic[!is.na(scores_pathogenic)]

# Calculate PR-AUC for REVEL - total genes
pr <- pr.curve(scores.class0 = scores_benign, scores.class1 = scores_pathogenic, curve = TRUE)
revel_aupr <- pr$auc.integral

# Subset for HCM genes
dt <- dt[gene_name %in% c("ACTC1", "MYBPC3", "MYH7", "MYL2", "MYL3", "PLN", "TNNI3", "TNNT2", "TPM1")]

# Calculate PR-AUC for Alphamissense - HCM genes

scores_pathogenic <- dt$AlphaMissense_pathogenicity[dt$ClinVar_ClinicalSignificance %in% c("Likely pathogenic", "Pathogenic", "Pathogenic/Likely pathogenic")]
scores_benign <- dt$AlphaMissense_pathogenicity[dt$ClinVar_ClinicalSignificance %in% c("Likely benign", "Benign", "Benign/Likely benign")]

scores_benign <- scores_benign[!is.na(scores_benign)]
scores_pathogenic <- scores_pathogenic[!is.na(scores_pathogenic)]

pr <- pr.curve(scores.class0 = scores_benign, scores.class1 = scores_pathogenic, curve = TRUE)
am_aupr_hcm <- pr$auc.integral


# Calculate PR-AUC for CardioBoost - HCM genes
scores_pathogenic <- dt$CardioBoost_pathogenicity[dt$ClinVar_ClinicalSignificance %in% c("Likely pathogenic", "Pathogenic", "Pathogenic/Likely pathogenic")]

scores_benign <- dt$CardioBoost_pathogenicity[dt$ClinVar_ClinicalSignificance %in% c("Likely benign", "Benign", "Benign/Likely benign")]

scores_benign <- scores_benign[!is.na(scores_benign)]
scores_pathogenic <- scores_pathogenic[!is.na(scores_pathogenic)]

pr <- pr.curve(scores.class0 = scores_benign, scores.class1 = scores_pathogenic, curve = TRUE)
cb_aupr_hcm <- pr$auc.integral

# Calculate PR-AUC for REVEL - HCM genes
scores_pathogenic <- dt$REVEL[dt$ClinVar_ClinicalSignificance %in% c("Likely pathogenic", "Pathogenic", "Pathogenic/Likely pathogenic")]

scores_benign <- dt$REVEL[dt$ClinVar_ClinicalSignificance %in% c("Likely benign", "Benign", "Benign/Likely benign")]

scores_benign <- scores_benign[!is.na(scores_benign)]
scores_pathogenic <- scores_pathogenic[!is.na(scores_pathogenic)]

pr <- pr.curve(scores.class0 = scores_benign, scores.class1 = scores_pathogenic, curve = TRUE)
revel_aupr_hcm <- pr$auc.integral

# Subset for DCM genes 

dt <- dt[gene_name %in% c("BAG3", "DES", "DSP", "LMNA", "MYH7", "PLN", "RBM20", "SCN5A", "TNNC1", "TNNT2", "TTN")]

# Calculate PR-AUC for AlphaMissense - DCM genes 

scores_pathogenic <- dt$AlphaMissense_pathogenicity[dt$ClinVar_ClinicalSignificance %in% c("Likely pathogenic", "Pathogenic", "Pathogenic/Likely pathogenic")]
scores_benign <- dt$AlphaMissense_pathogenicity[dt$ClinVar_ClinicalSignificance %in% c("Likely benign", "Benign", "Benign/Likely benign")]

scores_benign <- scores_benign[!is.na(scores_benign)]
scores_pathogenic <- scores_pathogenic[!is.na(scores_pathogenic)]

pr <- pr.curve(scores.class0 = scores_benign, scores.class1 = scores_pathogenic, curve = TRUE)
am_aupr_dcm <- pr$auc.integral

# Calculate PR-AUC for CardioBoost - DCM genes 

scores_pathogenic <- dt$CardioBoost_pathogenicity[dt$ClinVar_ClinicalSignificance %in% c("Likely pathogenic", "Pathogenic", "Pathogenic/Likely pathogenic")]

scores_benign <- dt$CardioBoost_pathogenicity[dt$ClinVar_ClinicalSignificance %in% c("Likely benign", "Benign", "Benign/Likely benign")]

scores_benign <- scores_benign[!is.na(scores_benign)]
scores_pathogenic <- scores_pathogenic[!is.na(scores_pathogenic)]

pr <- pr.curve(scores.class0 = scores_benign, scores.class1 = scores_pathogenic, curve = TRUE)
cb_aupr_dcm <- pr$auc.integral

# Calculate PR-AUC for Revel - DCM genes

scores_pathogenic <- dt$REVEL[dt$ClinVar_ClinicalSignificance %in% c("Likely pathogenic", "Pathogenic", "Pathogenic/Likely pathogenic")]

scores_benign <- dt$REVEL[dt$ClinVar_ClinicalSignificance %in% c("Likely benign", "Benign", "Benign/Likely benign")]

scores_benign <- scores_benign[!is.na(scores_benign)]
scores_pathogenic <- scores_pathogenic[!is.na(scores_pathogenic)]

pr <- pr.curve(scores.class0 = scores_benign, scores.class1 = scores_pathogenic, curve = TRUE)
revel_aupr_dcm <- pr$auc.integral

# Calculate PR-AUC for AM - individual genes 

# Subset for gene 
dt <- dt[gene_name %in% c("TTN")]

scores_pathogenic <- dt$AlphaMissense_pathogenicity[dt$ClinVar_ClinicalSignificance %in% c("Likely pathogenic", "Pathogenic", "Pathogenic/Likely pathogenic")]
scores_benign <- dt$AlphaMissense_pathogenicity[dt$ClinVar_ClinicalSignificance %in% c("Likely benign", "Benign", "Benign/Likely benign")]

scores_benign <- scores_benign[!is.na(scores_benign)]
scores_pathogenic <- scores_pathogenic[!is.na(scores_pathogenic)]

pr <- pr.curve(scores.class0 = scores_benign, scores.class1 = scores_pathogenic, curve = TRUE)

# Calculate PR-AUC for CB - individual genes 

# Subset for gene 
dt <- dt[gene_name %in% c("TTN")]

scores_pathogenic <- dt$CardioBoost_pathogenicity[dt$ClinVar_ClinicalSignificance %in% c("Likely pathogenic", "Pathogenic", "Pathogenic/Likely pathogenic")]

scores_benign <- dt$CardioBoost_pathogenicity[dt$ClinVar_ClinicalSignificance %in% c("Likely benign", "Benign", "Benign/Likely benign")]

scores_benign <- scores_benign[!is.na(scores_benign)]
scores_pathogenic <- scores_pathogenic[!is.na(scores_pathogenic)]

pr <- pr.curve(scores.class0 = scores_benign, scores.class1 = scores_pathogenic, curve = TRUE)

# Calculate PR-AUC for REVEL - individual genes 
dt <- fread("../data/analysis_table.csv.gz")


# Subset for gene 
dt <- dt[gene_name %in% c("TTN")]

scores_pathogenic <- dt$REVEL[dt$ClinVar_ClinicalSignificance %in% c("Likely pathogenic", "Pathogenic", "Pathogenic/Likely pathogenic")]

scores_benign <- dt$REVEL[dt$ClinVar_ClinicalSignificance %in% c("Likely benign", "Benign", "Benign/Likely benign")]

scores_benign <- scores_benign[!is.na(scores_benign)]
scores_pathogenic <- scores_pathogenic[!is.na(scores_pathogenic)]

pr <- pr.curve(scores.class0 = scores_benign, scores.class1 = scores_pathogenic, curve = TRUE)
revel_aupr_dcm <- pr$auc.integral


# Plot PR-AUC curves for Total genes
# Export
png(("../results/aupr_total_genes.png"), width = 12, height = 8, units = "in", res = 100)

# Create binaries for clinvar clinical significance 

dt$ClinVar_significance_binary[dt$ClinVar_ClinicalSignificance %in% c("Likely pathogenic", "Pathogenic", "Pathogenic/Likely pathogenic")] <- 1

dt$ClinVar_significance_binary[dt$ClinVar_ClinicalSignificance %in% c("Likely benign", "Benign", "Benign/Likely benign")] <- 0

# Remove NA values for pathogenicity scores

dt <- dt[!is.na(dt$AlphaMissense_pathogenicity), ]


pr_curve_logistic <- pr.curve(scores.class0 = dt$ClinVar_significance_binary, weights.class0 = dt$AlphaMissense_pathogenicity[dt$AlphaMissense_pathogenicity == 1], curve = TRUE)



# Plot PR-AUC curves - HCM Genes


# Plot PR-AUC curves - DCM Genes



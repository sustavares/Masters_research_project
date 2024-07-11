rm(list = ls())
graphics.off()

library(data.table)

# install.packages("pROC")
library(pROC)
library(PRROC)


# Import ========

dt <- fread("../data/analysis_table.csv.gz")

# Generating pathogenic scores and benign scores for AM
table(dt$ClinVar_ClinicalSignificance)
scores_pathogenic <- dt$AlphaMissense_pathogenicity[dt$ClinVar_ClinicalSignificance %in% c("Likely pathogenic", "Pathogenic", "Pathogenic/Likely pathogenic")]
scores_pathogenic <- dt$AlphaMissense_pathogenicity[dt$ClinVar_ClinicalSignificance %like% "athogenic"]

scores_benign <- dt$AlphaMissense_pathogenicity[dt$ClinVar_ClinicalSignificance %in% c("Likely benign", "Benign", "Benign/Likely benign")]

scores_benign <- scores_benign[!is.na(scores_benign)]
scores_pathogenic <- scores_pathogenic[!is.na(scores_pathogenic)]

# Calculate ROC-AUC for AM - total genes

roc <- roc.curve(scores.class0 = scores_pathogenic, scores.class1 = scores_benign, curve = TRUE)
am_auroc <- roc$auc


?roc.curve

# Generating pathogenic scores and benign scores for CB
scores_pathogenic <- dt$CardioBoost_pathogenicity[dt$ClinVar_ClinicalSignificance %in% c("Likely pathogenic", "Pathogenic", "Pathogenic/Likely pathogenic")]

scores_benign <- dt$CardioBoost_pathogenicity[dt$ClinVar_ClinicalSignificance %in% c("Likely benign", "Benign", "Benign/Likely benign")]

scores_benign <- scores_benign[!is.na(scores_benign)]
scores_pathogenic <- scores_pathogenic[!is.na(scores_pathogenic)]

# Calculate ROC-AUC for CB - total genes

roc <- roc.curve(scores.class0 = scores_pathogenic, scores.class1 = scores_benign, curve = TRUE)
cb_auroc <- roc$auc


# Generating pathogenic scores and benign scores for REVEL
scores_pathogenic <- dt$REVEL[dt$ClinVar_ClinicalSignificance %in% c("Likely pathogenic", "Pathogenic", "Pathogenic/Likely pathogenic")]

scores_benign <- dt$REVEL[dt$ClinVar_ClinicalSignificance %in% c("Likely benign", "Benign", "Benign/Likely benign")]

scores_benign <- scores_benign[!is.na(scores_benign)]
scores_pathogenic <- scores_pathogenic[!is.na(scores_pathogenic)]

# Calculate ROC-AUC for REVEL- total genes 

roc <- roc.curve(scores.class0 = scores_pathogenic, scores.class1 = scores_benign, curve = TRUE)
revel_auroc <- roc$auc

# Subset for HCM genes

dt <- dt[gene_name %in% c("ACTC1", "MYBPC3", "MYH7", "MYL2", "MYL3", "PLN", "TNNI3", "TNNT2", "TPM1")]

# Calculate ROC-AUC for AlphaMissense - HCM genes 

scores_pathogenic <- dt$AlphaMissense_pathogenicity[dt$ClinVar_ClinicalSignificance %in% c("Likely pathogenic", "Pathogenic", "Pathogenic/Likely pathogenic")]
scores_benign <- dt$AlphaMissense_pathogenicity[dt$ClinVar_ClinicalSignificance %in% c("Likely benign", "Benign", "Benign/Likely benign")]

scores_benign <- scores_benign[!is.na(scores_benign)]
scores_pathogenic <- scores_pathogenic[!is.na(scores_pathogenic)]

roc <- roc.curve(scores.class0 = scores_pathogenic, scores.class1 = scores_benign, curve = TRUE)
am_auroc_hcm <- roc$auc

# Calculate ROC-AUC for CardioBoost - HCM genes
scores_pathogenic <- dt$CardioBoost_pathogenicity[dt$ClinVar_ClinicalSignificance %in% c("Likely pathogenic", "Pathogenic", "Pathogenic/Likely pathogenic")]

scores_benign <- dt$CardioBoost_pathogenicity[dt$ClinVar_ClinicalSignificance %in% c("Likely benign", "Benign", "Benign/Likely benign")]

scores_benign <- scores_benign[!is.na(scores_benign)]
scores_pathogenic <- scores_pathogenic[!is.na(scores_pathogenic)]

roc <- roc.curve(scores.class0 = scores_pathogenic, scores.class1 = scores_benign, curve = TRUE)
cb_auroc_hcm <- roc$auc

# Calculate ROC-AUC for REVEL - HCM genes

scores_pathogenic <- dt$REVEL[dt$ClinVar_ClinicalSignificance %in% c("Likely pathogenic", "Pathogenic", "Pathogenic/Likely pathogenic")]

scores_benign <- dt$REVEL[dt$ClinVar_ClinicalSignificance %in% c("Likely benign", "Benign", "Benign/Likely benign")]

scores_benign <- scores_benign[!is.na(scores_benign)]
scores_pathogenic <- scores_pathogenic[!is.na(scores_pathogenic)]

roc <- roc.curve(scores.class0 = scores_pathogenic, scores.class1 = scores_benign, curve = TRUE)
revel_auroc_hcm <- roc$auc

# Subset for DCM genes 

dt <- dt[gene_name %in% c("BAG3", "DES", "DSP", "LMNA", "MYH7", "PLN", "RBM20", "SCN5A", "TNNC1", "TNNT2", "TTN")]

# Calculate ROC-AUC for AlphaMissense - DCM genes

scores_pathogenic <- dt$AlphaMissense_pathogenicity[dt$ClinVar_ClinicalSignificance %in% c("Likely pathogenic", "Pathogenic", "Pathogenic/Likely pathogenic")]
scores_benign <- dt$AlphaMissense_pathogenicity[dt$ClinVar_ClinicalSignificance %in% c("Likely benign", "Benign", "Benign/Likely benign")]

scores_benign <- scores_benign[!is.na(scores_benign)]
scores_pathogenic <- scores_pathogenic[!is.na(scores_pathogenic)]

roc <- roc.curve(scores.class0 = scores_pathogenic, scores.class1 = scores_benign, curve = TRUE)
am_auroc_dcm <- roc$auc

# Calculate ROC-AUC for CardioBoost - DCM genes
scores_pathogenic <- dt$CardioBoost_pathogenicity[dt$ClinVar_ClinicalSignificance %in% c("Likely pathogenic", "Pathogenic", "Pathogenic/Likely pathogenic")]

scores_benign <- dt$CardioBoost_pathogenicity[dt$ClinVar_ClinicalSignificance %in% c("Likely benign", "Benign", "Benign/Likely benign")]

scores_benign <- scores_benign[!is.na(scores_benign)]
scores_pathogenic <- scores_pathogenic[!is.na(scores_pathogenic)]

roc <- roc.curve(scores.class0 = scores_pathogenic, scores.class1 = scores_benign, curve = TRUE)
cb_auroc_dcm <- roc$auc


# Calculate ROC-AUC for REVEL - DCM genes

scores_pathogenic <- dt$REVEL[dt$ClinVar_ClinicalSignificance %in% c("Likely pathogenic", "Pathogenic", "Pathogenic/Likely pathogenic")]

scores_benign <- dt$REVEL[dt$ClinVar_ClinicalSignificance %in% c("Likely benign", "Benign", "Benign/Likely benign")]

scores_benign <- scores_benign[!is.na(scores_benign)]
scores_pathogenic <- scores_pathogenic[!is.na(scores_pathogenic)]

roc <- roc.curve(scores.class0 = scores_pathogenic, scores.class1 = scores_benign, curve = TRUE)
revel_auroc_dcm <- roc$auc

# Calculate ROC-AUC for AM for individual genes

# Subset for Individual genes

dt <- dt[gene_name %in% c("TTN")]

scores_pathogenic <- dt$AlphaMissense_pathogenicity[dt$ClinVar_ClinicalSignificance %in% c("Likely pathogenic", "Pathogenic", "Pathogenic/Likely pathogenic")]
scores_benign <- dt$AlphaMissense_pathogenicity[dt$ClinVar_ClinicalSignificance %in% c("Likely benign", "Benign", "Benign/Likely benign")]

scores_benign <- scores_benign[!is.na(scores_benign)]
scores_pathogenic <- scores_pathogenic[!is.na(scores_pathogenic)]

roc <- roc.curve(scores.class0 = scores_pathogenic, scores.class1 = scores_benign, curve = TRUE)
am_auroc_TTN <- roc$auc

# Calculate ROC-AUC for CB for individual genes
dt <- dt[gene_name %in% c("TTN")]

scores_pathogenic <- dt$CardioBoost_pathogenicity[dt$ClinVar_ClinicalSignificance %in% c("Likely pathogenic", "Pathogenic", "Pathogenic/Likely pathogenic")]

scores_benign <- dt$CardioBoost_pathogenicity[dt$ClinVar_ClinicalSignificance %in% c("Likely benign", "Benign", "Benign/Likely benign")]

scores_benign <- scores_benign[!is.na(scores_benign)]
scores_pathogenic <- scores_pathogenic[!is.na(scores_pathogenic)]

roc <- roc.curve(scores.class0 = scores_pathogenic, scores.class1 = scores_benign, curve = TRUE)

# Calculate ROC-AUC for REVEL for individual genes 

dt <- dt[gene_name %in% c("TTN")]
scores_pathogenic <- dt$REVEL[dt$ClinVar_ClinicalSignificance %in% c("Likely pathogenic", "Pathogenic", "Pathogenic/Likely pathogenic")]

scores_benign <- dt$REVEL[dt$ClinVar_ClinicalSignificance %in% c("Likely benign", "Benign", "Benign/Likely benign")]

scores_benign <- scores_benign[!is.na(scores_benign)]
scores_pathogenic <- scores_pathogenic[!is.na(scores_pathogenic)]

roc <- roc.curve(scores.class0 = scores_pathogenic, scores.class1 = scores_benign, curve = TRUE)


# Plot ROC curve - Total genes

# Export
png(("../results/auroc_total_genes.png"), width = 9, height = 4.5, units = "in", res = 100)

# Create binaries for clinvar clinical significance 

dt$ClinVar_significance_binary[dt$ClinVar_ClinicalSignificance %in% c("Likely pathogenic", "Pathogenic", "Pathogenic/Likely pathogenic")] <- 1

dt$ClinVar_significance_binary[dt$ClinVar_ClinicalSignificance %in% c("Likely benign", "Benign", "Benign/Likely benign")] <- 0

roc_logistic_am <- roc(dt$ClinVar_significance_binary, dt$AlphaMissense_pathogenicity)
auroc_logistic_am <- auc(roc_logistic_am)


roc_logistic_cb <- roc(dt$ClinVar_significance_binary, dt$CardioBoost_pathogenicity)
auroc_logistic_cb <- auc(roc_logistic_cb)


roc_logistic_revel <- roc(dt$ClinVar_significance_binary, dt$REVEL)
auroc_logistic_revel <- auc(roc_logistic_revel)

# Create legend text 

legend_text_auroc <- c(paste("AlphaMissense (AUROC:", format(auroc_logistic_am, digits = 2), ")"), paste("CardioBoost (AUROC:", format(auroc_logistic_cb, digits = 2),")"), paste("Revel (AUROC:", format(auroc_logistic_revel, digits = 2),")"))

# Plot roc curve 

plot(roc_logistic_am, main ="ROC Total Genes", col = "#89CFF0")

# Add CardioBoost and REVEL curves

lines(roc_logistic_cb, type = "l", col = "#009E73")
lines(roc_logistic_revel, type = "l", col = "#FF69B4")


legend("bottomright", legend = legend_text_auroc, col = c("#89CFF0", "#009E73", "#FF69B4"), lty = 1, lwd = 2, cex = 0.8, bty = "n")

dev.off()

# Plot ROC curves for HCM genes

# Subset for HCM genes 

dt <- dt[gene_name %in% c("ACTC1", "MYBPC3", "MYH7", "MYL2", "MYL3", "PLN", "TNNI3", "TNNT2", "TPM1")]

dt$ClinVar_significance_binary[dt$ClinVar_ClinicalSignificance %in% c("Likely pathogenic", "Pathogenic", "Pathogenic/Likely pathogenic")] <- 1

dt$ClinVar_significance_binary[dt$ClinVar_ClinicalSignificance %in% c("Likely benign", "Benign", "Benign/Likely benign")] <- 0

roc_logistic_am <- roc(dt$ClinVar_significance_binary, dt$AlphaMissense_pathogenicity)
auroc_logistic_am <- auc(roc_logistic_am)


roc_logistic_cb <- roc(dt$ClinVar_significance_binary, dt$CardioBoost_pathogenicity)
auroc_logistic_cb <- auc(roc_logistic_cb)


roc_logistic_revel <- roc(dt$ClinVar_significance_binary, dt$REVEL)
auroc_logistic_revel <- auc(roc_logistic_revel)

# Export
png(("../results/auroc_hcm_genes.png"), width = 9, height = 4.5, units = "in", res = 100)


# Create legend text 

legend_text_auroc <- c(paste("AlphaMissense (AUROC:", format(auroc_logistic_am, digits = 2), ")"), paste("CardioBoost (AUROC:", format(auroc_logistic_cb, digits = 2),")"), paste("Revel (AUROC:", format(auroc_logistic_revel, digits = 2),")"))

# Plot roc curve 

plot(roc_logistic_am, main ="ROC HCM Genes", col = "#89CFF0")

# Add CardioBoost and REVEL curves

lines(roc_logistic_cb, type = "l", col = "#009E73")
lines(roc_logistic_revel, type = "l", col = "#FF69B4")


legend("bottomright", legend = legend_text_auroc, col = c("#89CFF0", "#009E73", "#FF69B4"), lty = 1, lwd = 2, cex = 0.8, bty = "n")

dev.off()

# Plot ROC curves for DCM genes

# Subset for DCM genes

dt <- dt[gene_name %in% c("BAG3", "DES", "DSP", "LMNA", "MYH7", "PLN", "RBM20", "SCN5A", "TNNC1", "TNNT2", "TTN")]

# Export
png(("../results/auroc_dcm_genes.png"), width = 9, height = 4.5, units = "in", res = 100)


dt$ClinVar_significance_binary[dt$ClinVar_ClinicalSignificance %in% c("Likely pathogenic", "Pathogenic", "Pathogenic/Likely pathogenic")] <- 1

dt$ClinVar_significance_binary[dt$ClinVar_ClinicalSignificance %in% c("Likely benign", "Benign", "Benign/Likely benign")] <- 0

roc_logistic_am <- roc(dt$ClinVar_significance_binary, dt$AlphaMissense_pathogenicity)
auroc_logistic_am <- auc(roc_logistic_am)


roc_logistic_cb <- roc(dt$ClinVar_significance_binary, dt$CardioBoost_pathogenicity)
auroc_logistic_cb <- auc(roc_logistic_cb)


roc_logistic_revel <- roc(dt$ClinVar_significance_binary, dt$REVEL)
auroc_logistic_revel <- auc(roc_logistic_revel)

# Create legend text 

legend_text_auroc <- c(paste("AlphaMissense (AUROC:", format(auroc_logistic_am, digits = 2), ")"), paste("CardioBoost (AUROC:", format(auroc_logistic_cb, digits = 2),")"), paste("Revel (AUROC:", format(auroc_logistic_revel, digits = 2),")"))

# Plot roc curve 

plot(roc_logistic_am, main ="ROC DCM Genes", col = "#89CFF0")

# Add CardioBoost and REVEL curves

lines(roc_logistic_cb, type = "l", col = "#009E73")
lines(roc_logistic_revel, type = "l", col = "#FF69B4")


legend("bottomright", legend = legend_text_auroc, col = c("#89CFF0", "#009E73", "#FF69B4"), lty = 1, lwd = 2, cex = 0.8, bty = "n")

dev.off()





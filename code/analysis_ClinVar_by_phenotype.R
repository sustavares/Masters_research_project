# Table the number of pathogenic and benign variants associated with HCM or DCM per gene
# Calculate precision, recall, and specificity for each score on each gene
# Calculate AUROC and AUPRC for each score on each gene
# Table and plot results

rm(list = ls())
graphics.off()

library(data.table)
library(pROC)
library(PRROC)

### Import ==========

TTN_exons <- fread("../data/TTN_highPSI_GRCh38.bed")

dt <- fread("../data/analysis_table.csv.gz")

dt_model <- fread("/Users/susannaantunes/research project/cardio_path_scores/data/analysis_table.csv.gz")

### Set variables ==========

# DCM genes with evidence of missense pathogenicity
dcm_genes <- c("BAG3", "DES", "DSP", "LMNA", "MYH7", "PLN", "RBM20", "SCN5A", "TNNC1", "TNNT2", "TTN")

# HCM genes with evidence of missense pathogenicity
hcm_genes <- c("ACTC1", "MYBPC3", "MYH7", "MYL2", "MYL3", "PLN", "TNNI3", "TNNT2", "TPM1")


### Preprocessing ==========

# Code clinical significance 
dt$ClinVar_significance_binary <- NA
dt$ClinVar_significance_binary[dt$ClinVar_significance %in% c("Likely pathogenic",
                                                              "Pathogenic",
                                                              "Pathogenic/Likely pathogenic")] <- 1
dt$ClinVar_significance_binary[dt$ClinVar_significance %in% c("Likely benign",
                                                              "Benign",
                                                              "Benign/Likely benign")] <- 0

# Subset HCM and DCM genes with evidence of missnese pathogenicity
dt <- dt[consequence == "missense" & 
         ClinVar_significance_binary %in% c(0, 1) &
         gene_name %in% c(hcm_genes, dcm_genes), 
         c("gene_name", "chromosome", "pos", "ref", "alt", 
           "ClinVar_significance_binary", "ClinVar_phenotype",
           "AlphaMissense_pathogenicity", "AlphaMissense_class", 
           "CardioBoost_pathogenicity",      
           "REVEL")]

# Set score binary predictions based on recommended thresholds

dt$AlphaMissense_binary <- NA
dt$AlphaMissense_binary[dt$AlphaMissense_class == "likely_pathogenic"] <- 1
dt$AlphaMissense_binary[dt$AlphaMissense_class == "likely_benign"] <- 0

dt$CardioBoost_binary <- NA
dt$CardioBoost_binary[dt$CardioBoost_pathogenicity >= 0.9] <- 1
dt$CardioBoost_binary[dt$CardioBoost_pathogenicity <= 0.1] <- 0

dt$REVEL_binary <- NA
dt$REVEL_binary[dt$REVEL >= 0.75] <- 1
dt$REVEL_binary[dt$REVEL <= 0.25] <- 0

# Subset HCM genes and pathogenic variants
dt_hcm <- dt[ClinVar_phenotype == "HCM" & gene_name %in% hcm_genes]

# Subset DCM genes and pathogenic variants
dt_dcm <- dt[ClinVar_phenotype == "DCM" & gene_name %in% dcm_genes]


### Table pathogenic and benign ClinVar variants per gene ==========

# Table pathogenic and benign counts per HCM gene
dt_n_1 <- as.data.table(table(dt_hcm$gene_name[dt_hcm$ClinVar_significance_binary == 1]))
dt_n_2 <- as.data.table(table(dt_hcm$gene_name[dt_hcm$ClinVar_significance_binary == 0]))
colnames(dt_n_1) <- c("gene", "n_pathogenic")
colnames(dt_n_2) <- c("gene", "n_benign")
dt_counts_hcm <- merge(dt_n_1, dt_n_2, all = TRUE)
dt_counts_hcm[is.na(dt_counts_hcm)] <- 0
dt_counts_hcm$disease <- "HCM"

# Table pathogenic and benign counts per DCM gene
dt_n_1 <- as.data.table(table(dt_dcm$gene_name[dt_dcm$ClinVar_significance_binary == 1]))
dt_n_2 <- as.data.table(table(dt_dcm$gene_name[dt_dcm$ClinVar_significance_binary == 0]))
colnames(dt_n_1) <- c("gene", "n_pathogenic")
colnames(dt_n_2) <- c("gene", "n_benign")
dt_counts_dcm <- merge(dt_n_1, dt_n_2, all = TRUE)
dt_counts_dcm[is.na(dt_counts_dcm)] <- 0
dt_counts_dcm$disease <- "DCM"
dt_counts <- rbindlist(list(dt_counts_hcm,
                            dt_counts_dcm,
                                  data.table(gene = c("total", "total", "total"),
                                             disease = c("HCM", "DCM", "HCM and DCM"),
                                             n_pathogenic = c(sum(dt_counts_hcm$n_pathogenic, na.rm = TRUE),
                                                              sum(dt_counts_dcm$n_pathogenic, na.rm = TRUE),
                                                              sum(dt_counts_hcm$n_pathogenic, dt_counts_dcm$n_pathogenic, na.rm = TRUE)),
                                             n_benign = c(sum(dt_counts_hcm$n_benign, na.rm = TRUE),
                                                          sum(dt_counts_dcm$n_benign, na.rm = TRUE),
                                                          sum(dt_counts_hcm$n_benign, dt_counts_dcm$n_benign, na.rm = TRUE)))
                            ), fill = TRUE)



### Evaluate performance ==========

# score <- "AlphaMissense_binary"
# Function that calculates precision, recall, and specificity
eval_metrics <- function(dt, score){

  # Create score vectors and remove NAs
  score_vec <- dt[, get(score)]
  cv_vec <- dt$ClinVar_significance_binary
  cv_vec <- cv_vec[!is.na(score_vec)]
  score_vec <- score_vec[!is.na(score_vec)]

  # Calculate true positives, false positives, true negatives, and false negatives
  TP <- sum(cv_vec)
  FP <- sum(score_vec == 1 & cv_vec == 0)
  TN <- sum(cv_vec == 0)
  FN <- sum(score_vec == 0 & cv_vec == 1)

  # Calculate precision, recall, and specificity
  precision <- TP / (TP + FP)
  recall <- TP / (TP + FN)
  specificity <- TN / (TN + FP)
  
  return(c(TP, FP, TN, FN, precision, recall, specificity))
}

# Calculate precision, recall, and specificity for each score on each HCM gene
for (i in seq_along(unique(hcm_genes))){
  gene <- unique(hcm_genes)[i]
  dt_sub <- dt_hcm[gene_name == gene]

  for (score in c("AlphaMissense_binary", "CardioBoost_binary", "REVEL_binary")){
    metrics <- eval_metrics(dt_sub, score)

    if (i == 1 & score == "AlphaMissense_binary"){
      dt_metrics_hcm <- data.table(gene = gene,
                               score = score,
                               TP = metrics[1],
                               FP = metrics[2],
                               TN = metrics[3],
                               FN = metrics[4],
                               precision = metrics[5],
                               recall = metrics[6],
                               specificity = metrics[7])
    } else {
      dt_metrics_hcm <- rbindlist(list(dt_metrics_hcm,
                                  data.table(gene = gene,
                                             score = score,
                                             TP = metrics[1],
                                             FP = metrics[2],
                                             TN = metrics[3],
                                             FN = metrics[4],
                                             precision = metrics[5],
                                             recall = metrics[6],
                                             specificity = metrics[7])),
                              fill = TRUE)
    }
  }
}
dt_metrics_hcm$disease <- "HCM"

# Calculate precision, recall, and specificity for each score on each DCM gene
for (i in seq_along(unique(dcm_genes))){
  gene <- unique(dcm_genes)[i]
  dt_sub <- dt_dcm[gene_name == gene]
  
  for (score in c("AlphaMissense_binary", "CardioBoost_binary", "REVEL_binary")){
    metrics <- eval_metrics(dt_sub, score)
    
    if (i == 1 & score == "AlphaMissense_binary"){
      dt_metrics_dcm <- data.table(gene = gene,
                                   score = score,
                                   TP = metrics[1],
                                   FP = metrics[2],
                                   TN = metrics[3],
                                   FN = metrics[4],
                                   precision = metrics[5],
                                   recall = metrics[6],
                                   specificity = metrics[7])
    } else {
      dt_metrics_dcm <- rbindlist(list(dt_metrics_dcm,
                                       data.table(gene = gene,
                                                  score = score,
                                                  TP = metrics[1],
                                                  FP = metrics[2],
                                                  TN = metrics[3],
                                                  FN = metrics[4],
                                                  precision = metrics[5],
                                                  recall = metrics[6],
                                                  specificity = metrics[7])),
                                  fill = TRUE)
    }
  }
}
dt_metrics_dcm$disease <- "DCM"

# Calculate precision, recall, and specificity for each score by all HCM genes, all DCM genes, and all HCM and DCM genes
gene_sets <- list(dt_hcm, dt_dcm, rbind(dt_hcm, dt_dcm))
disease_groups <- c("HCM", "DCM", "HCM and DCM")
for (i in seq_along(gene_sets)){
  dt_sub <- gene_sets[[i]]
  
  for (score in c("AlphaMissense_binary", "CardioBoost_binary", "REVEL_binary")){
    metrics <- eval_metrics(dt_sub, score)
    
    if (i == 1 & score == "AlphaMissense_binary"){
      dt_metrics_all <- data.table(disease = disease_groups[i],
                               gene = "total",
                               score = score,
                               TP = metrics[1],
                               FP = metrics[2],
                               TN = metrics[3],
                               FN = metrics[4],
                               precision = metrics[5],
                               recall = metrics[6],
                               specificity = metrics[7])
    } else {
      dt_metrics_all <- rbindlist(list(dt_metrics_all,
                                   data.table(disease = disease_groups[i],
                                              gene = "total",
                                              score = score,
                                              TP = metrics[1],
                                              FP = metrics[2],
                                              TN = metrics[3],
                                              FN = metrics[4],
                                              precision = metrics[5],
                                              recall = metrics[6],
                                              specificity = metrics[7])),
                              fill = TRUE)
    }
  }
}
dt_metrics <- rbindlist(list(dt_metrics_hcm, dt_metrics_dcm, dt_metrics_all), fill = TRUE)
dt_metrics$score <- gsub("_binary", "", dt_metrics$score)

### Calculate AUROC and AUPRC ==========

# Function that calculates AUROC and AUPRC
# score <- "AlphaMissense_pathogenicity"
eval_auc <- function(dt, score) {
  
  # Create score vectors and remove NAs
  # score_vec <- dt[, get(score), with = FALSE][[1]]  
  score_vec <- dt[, .SD[[1]], .SDcols = score] # Access the first column directly
  cv_vec <- dt$ClinVar_significance_binary
  cv_vec <- cv_vec[!is.na(score_vec)]
  score_vec <- score_vec[!is.na(score_vec)]
  
  # Check if 'cv_vec' has exactly two levels
  if (length(unique(cv_vec)) != 2) {
    return(list(data.table(score = score, auroc = NA, prauc = NA),
                roc_curve = NA,
                pr_curve = NA))
  }
  
  # Calculate AUROC
  roc_curve <- roc(cv_vec, score_vec)
  au_roc <- auc(roc_curve)
  
  # Calculate AUPRC
  pr_curve <- pr.curve(scores.class0 = score_vec, weights.class0 = cv_vec, curve = TRUE)
  pr_auc <- pr_curve$auc.integral
  
  # Return structured output
  return(list(data.table(score = score, auroc = au_roc, prauc = pr_auc),
              roc_curve = roc_curve,
              pr_curve = pr_curve))
}

# Initialize a list to store AUROC and AUPRC tables and curves
results_list <- list()

# Calculate AUROC and AUPRC for each score by all HCM genes, all DCM genes, and all HCM and DCM genes
gene_sets <- list(dt_hcm, dt_dcm, rbind(dt_hcm, dt_dcm))
disease_groups <- c("HCM", "DCM", "HCM and DCM")

for (i in seq_along(gene_sets)){
  dt_sub <- gene_sets[[i]]
  
  for (score in c("AlphaMissense_pathogenicity", "CardioBoost_pathogenicity", "REVEL")){
    auc_metrics <- eval_auc(dt_sub, score)
    
    # Store metrics in the list with a unique identifier
    results_list[[paste(disease_groups[i], score, "metrics")]] <- auc_metrics[[1]]
    results_list[[paste(disease_groups[i], score, "roc_curve")]] <- auc_metrics[[2]]
    results_list[[paste(disease_groups[i], score, "pr_curve")]] <- auc_metrics[[3]]
    print(score)
  }
  print(i)
}

# Filter to get only every third element for auc values
auc_list <- results_list[seq(1, length(results_list), by = 3)]
# Convert the second column of all data tables to character
auc_list <- lapply(auc_list, function(dt) {
  dt[, 2 := as.character(get(names(dt)[2]))]
  return(dt)
})
dt_auc <- rbindlist(auc_list, use.names = TRUE, fill = TRUE)
dt_auc$prauc <- as.numeric(dt_auc$prauc)
dt_auc$auroc <- as.numeric(dt_auc$auroc)
dt_auc$disease <- rep(c("HCM", "DCM", "HCM and DCM"), each = 3)
dt_auc$score <- gsub("_pathogenicity", "", dt_auc$score)


### Plot AUROC and AUPRC ==========

auroc_list <- results_list[seq(2, length(results_list), by = 3)]
prauc_list <- results_list[seq(3, length(results_list), by = 3)]

# Create legend text
legend_text_auroc_hcm <- c(
  paste(dt_auc$score[1], " (AUROC:", format(dt_auc$auroc[1], digits = 2), ")"),
  paste(dt_auc$score[2], " (AUROC:", format(dt_auc$auroc[2], digits = 2), ")"),
  paste(dt_auc$score[3], " (AUROC:", format(dt_auc$auroc[3], digits = 2), ")")
)
legend_text_auroc_dcm <- c(
  paste(dt_auc$score[4], " (AUROC:", format(dt_auc$auroc[4], digits = 2), ")"),
  paste(dt_auc$score[5], " (AUROC:", format(dt_auc$auroc[5], digits = 2), ")"),
  paste(dt_auc$score[6], " (AUROC:", format(dt_auc$auroc[6], digits = 2), ")")
)
legend_text_auroc_all <- c(
  paste(dt_auc$score[7], " (AUROC:", format(dt_auc$auroc[7], digits = 2), ")"),
  paste(dt_auc$score[8], " (AUROC:", format(dt_auc$auroc[8], digits = 2), ")"),
  paste(dt_auc$score[9], " (AUROC:", format(dt_auc$auroc[9], digits = 2), ")")
)
legend_text_auprc_hcm <- c(
  paste(dt_auc$score[1], " (AUPRC:", format(dt_auc$prauc[1], digits = 2), ")"),
  paste(dt_auc$score[2], " (AUPRC:", format(dt_auc$prauc[2], digits = 2), ")"),
  paste(dt_auc$score[3], " (AUPRC:", format(dt_auc$prauc[3], digits = 2), ")")
)
legend_text_auprc_dcm <- c(
  paste(dt_auc$score[4], " (AUPRC:", format(dt_auc$prauc[4], digits = 2), ")"),
  paste(dt_auc$score[5], " (AUPRC:", format(dt_auc$prauc[5], digits = 2), ")"),
  paste(dt_auc$score[6], " (AUPRC:", format(dt_auc$prauc[6], digits = 2), ")")
)
legend_text_auprc_all <- c(
  paste(dt_auc$score[7], " (AUPRC:", format(dt_auc$prauc[7], digits = 2), ")"),
  paste(dt_auc$score[8], " (AUPRC:", format(dt_auc$prauc[8], digits = 2), ")"),
  paste(dt_auc$score[9], " (AUPRC:", format(dt_auc$prauc[9], digits = 2), ")")
)


# Specify the file to save the plot
png("../results/Figure_auroc_auprc_by_phenotype.png", width = 12, height = 8, units = "in", res = 100)

# Set global text size parameters
par(mfrow = c(2, 3),  # 1 row, 2 columns
    cex = 1,         # Increase text size for content
    # cex.main = 1.5,    # Increase main title text size
    cex.axis = 1,    # Increase axis text size
    cex.lab = 1)     # Increase axis label text size

## Plot AUROC HCM
plot(auroc_list[[1]], col = "#E69F00", lwd = 2, main = "A) ROC Curves (HCM)", percent = TRUE)
try(lines(auroc_list[[2]], col = "#56B4E9", lwd = 2, type = "l"))
try(lines(auroc_list[[3]], col = "#009E73", lwd = 2, type = "l"))
legend("bottomright", legend = legend_text_auroc_hcm,
       col = c("#E69F00", "#56B4E9", "#009E73"), lty = 1, lwd = 2, cex = 0.8, bty = "n")

## Plot AUROC DCM
plot(auroc_list[[4]], col = "#E69F00", lwd = 2, main = "B) ROC Curves (DCM)", percent = TRUE)
try(lines(auroc_list[[5]], col = "#56B4E9", lwd = 2, type = "l"))
try(lines(auroc_list[[6]], col = "#009E73", lwd = 2, type = "l"))
legend("bottomright", legend = legend_text_auroc_dcm,
       col = c("#E69F00", "#56B4E9", "#009E73"), lty = 1, lwd = 2, cex = 0.8, bty = "n")

## Plot AUROC All
plot(auroc_list[[7]], col = "#E69F00", lwd = 2, main = "C) ROC Curves (HCM and DCM)", percent = TRUE)
try(lines(auroc_list[[8]], col = "#56B4E9", lwd = 2, type = "l"))
try(lines(auroc_list[[9]], col = "#009E73", lwd = 2, type = "l"))
legend("bottomright", legend = legend_text_auroc_all,
       col = c("#E69F00", "#56B4E9", "#009E73"), lty = 1, lwd = 2, cex = 0.8, bty = "n")

## Plot AUPRC HCM
plot(prauc_list[[1]]$curve, col = "#E69F00", lwd = 2, type = "l", xlim = c(0, 1), ylim = c(0, 1),
     xlab = "Recall", ylab = "Precision",
     main = "C) P-R Curves (HCM)")
try(lines(prauc_list[[2]]$curve, col = "#56B4E9", lwd = 2, type = "l"))
try(lines(prauc_list[[3]]$curve, col = "#009E73", lwd = 2, type = "l"))
legend("bottomleft", legend = legend_text_auprc_hcm,
       col = c("#E69F00", "#56B4E9", "#009E73"), lty = 1, lwd = 2, cex = 0.8, bty = "n")

## Plot AUPRC DCM
plot(prauc_list[[4]]$curve, col = "#E69F00", lwd = 2, type = "l", xlim = c(0, 1), ylim = c(0, 1),
     xlab = "Recall", ylab = "Precision",
     main = "D) P-R Curves (DCM)")
try(lines(prauc_list[[5]]$curve, col = "#56B4E9", lwd = 2, type = "l"))
try(lines(prauc_list[[6]]$curve, col = "#009E73", lwd = 2, type = "l"))
legend("bottomleft", legend = legend_text_auprc_dcm,
       col = c("#E69F00", "#56B4E9", "#009E73"), lty = 1, lwd = 2, cex = 0.8, bty = "n")

## Plot AUPRC All
plot(prauc_list[[7]]$curve, col = "#E69F00", lwd = 2, type = "l", xlim = c(0, 1), ylim = c(0, 1),
     xlab = "Recall", ylab = "Precision",
     main = "E) P-R Curves (HCM and DCM)")
try(lines(prauc_list[[8]]$curve, col = "#56B4E9", lwd = 2, type = "l"))
try(lines(prauc_list[[9]]$curve, col = "#009E73", lwd = 2, type = "l"))
legend("bottomleft", legend = legend_text_auprc_all,
       col = c("#E69F00", "#56B4E9", "#009E73"), lty = 1, lwd = 2, cex = 0.8, bty = "n")

# Close the device to save the file
dev.off()

### Export =========

fwrite(dt_counts, "../results/Table_counts_allPB.csv")
fwrite(dt_metrics, "../results/Table_evaluation_metrics_allPB.csv")
fwrite(dt_auc, "../results/Table_auc_allPB.csv")


# Subset for TTN

# create vec of all positions in exons of interest
all_positions <- unlist(mapply(seq, from = TTN_exons$V2, to = TTN_exons$V3))

# Subset ttn 
dt_ttn <- dt[gene_name %in% c("TTN")]
dt_nottn <- dt[gene_name != "TTN"]

# Subset variants in exons of interest
dt_ttn <- dt_ttn[pos %in% all_positions]

# Put back together
dt <- rbind(dt_nottn, dt_ttn)

# Exclude TTN exons 
dt <- subset(dt, gene_name != 'TTN')

# Subset for only TTN gene
dt <- subset(dt, gene_name == 'TTN')



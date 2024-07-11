# Table counts of pathogenic and benign variants for each score (each cohort and total)
## Exclude variants with MAF > 0.001 (UKB control)
## Exclude variants classified as pathogenic or benign by ClinVar.
# Calculate OR and 95% CIs in cases vs controls. 

rm(list = ls())
graphics.off()

library(data.table)
library(epitools)
library(ggplot2)
library(ggpubr)

### Import ==========

dt <- fread("../data/analysis_table.csv.gz")

### Preprocessing ==========

# Calculate All cohort totals (VFX and UKB), treating NAs as 0s
dt$ALL_HCM_AC <- replace(dt$VFX_ALL_HCM_AC, is.na(dt$VFX_ALL_HCM_AC), 0) + 
                  replace(dt$UKB_hcm_all_AC, is.na(dt$UKB_hcm_all_AC), 0) 
dt$ALL_HCM_AN <- replace(dt$VFX_ALL_HCM_AN, is.na(dt$VFX_ALL_HCM_AN), 0) + 
                  replace(dt$UKB_hcm_all_AN, is.na(dt$UKB_hcm_all_AN), 0)

dt$ALL_DCM_AC <- replace(dt$VFX_ALL_DCM_AC, is.na(dt$VFX_ALL_DCM_AC), 0) + 
                  replace(dt$UKB_dcm_all_AC, is.na(dt$UKB_dcm_all_AC), 0) +
                  replace(dt$goDCM_EUR_AC, is.na(dt$goDCM_EUR_AC), 0)
dt$ALL_DCM_AN <- replace(dt$VFX_ALL_DCM_AN, is.na(dt$VFX_ALL_DCM_AN), 0) + 
                  replace(dt$UKB_dcm_all_AN, is.na(dt$UKB_dcm_all_AN), 0) +
                  replace(dt$goDCM_EUR_AN, is.na(dt$goDCM_EUR_AN), 0)


dt$ALL_HCM_DCM_AC <- dt$ALL_HCM_AC + dt$ALL_DCM_AC
dt$ALL_HCM_DCM_AN <- dt$ALL_HCM_AN + dt$ALL_DCM_AN

# Set score binary predictions based on recommended thresholds

dt$AlphaMissense_binary <- "VUS"
dt$AlphaMissense_binary[dt$AlphaMissense_class == "likely_pathogenic"] <- "P"
dt$AlphaMissense_binary[dt$AlphaMissense_class == "likely_benign"] <- "B"

dt$CardioBoost_binary <- "VUS"
dt$CardioBoost_binary[dt$CardioBoost_pathogenicity >= 0.9] <- "P"
dt$CardioBoost_binary[dt$CardioBoost_pathogenicity <= 0.1] <- "B"

dt$REVEL_binary <- "VUS"
dt$REVEL_binary[dt$REVEL >= 0.75] <- "P"
dt$REVEL_binary[dt$REVEL <= 0.75] <- "B"

dt$CADD_binary <- "VUS"
dt$CADD_binary[dt$CADD_PHRED >= 20] <- "P"
dt$CADD_binary[dt$CADD_PHRED <= 10] <- "B"


### QC ==========

# Exclude variants with MAF > 0.0001 (gnomAD)
dt$gnomAD_AF <- as.numeric(dt$gnomAD_AC) / as.numeric(dt$gnomAD_AN)
dt$gnomAD_AF[is.na(dt$gnomAD_AF)] <- 0
dt <- dt[gnomAD_AF <= 0.0001,]

# Exclude variants classified as pathogenic or benign by ClinVar.
dt <- dt[!ClinVar_significance %in% c("Benign",
                                     "Benign/Likely benign",
                                     "Likely benign",
                                     "Likely pathogenic",
                                     "Pathogenic Pathogenic/Likely pathogenic")]

# Subset synonymous and missense variants
dt_syn <- dt[consequence == "synonymous"]
dt <- dt[consequence == "missense"]

### Set variables ==========

v_cohort_hcm <- c("UKB_hcm_all",
                  "VFX_BLF_HCM",
                  "VFX_EGY_HCM", 
                  "VFX_GDX_HCM", 
                  "VFX_LMM_HCM", 
                  "VFX_OMG_HCM", 
                  "VFX_RBH_HCM",
                  "VFX_SNG_HCM", 
                  "ALL_HCM")

v_cohort_dcm <- c("UKB_dcm_all", 
                  "goDCM_EUR", 
                  "VFX_EGY_DCM",
                  "VFX_LMM_DCM",
                  "VFX_OMG_DCM",
                  "VFX_RBH_DCM",
                  "VFX_SNG_DCM",
                  "ALL_DCM")

# vec of scores
v_scores <- c("AlphaMissense", "CardioBoost", "REVEL", "CADD")

# DCM genes with evidence of missense pathogenicity
dcm_genes <- c("BAG3", "DES", "DSP", "LMNA", "MYH7", "PLN", "RBM20", "SCN5A", "TNNC1", "TNNT2", "TTN")

# HCM genes with evidence of missense pathogenicity
hcm_genes <- c("ACTC1", "MYBPC3", "MYH7", "MYL2", "MYL3", "PLN", "TNNI3", "TNNT2", "TPM1")


### Functions ==========

# Function to manually compute OR, CI, Fisher and Chi-Square p-values for each cohort
compute_or_manual <- function(AC_case, AN_case, AC_control, AN_control) {
  
  A <- as.numeric(AC_case)
  B <- as.numeric(AN_case - AC_case)
  C <- as.numeric(AC_control)
  D <- as.numeric(AN_control - AC_control)
  
  # Calculate Odds Ratio (OR)
  OR <- (A * D) / (B * C)
  
  # Standard Error (SE) of log(OR)
  SE_log_OR <- sqrt(1 / A + 1 / B + 1 / C + 1 / D)
  
  # Confidence Interval of log(OR)
  log_OR <- log(OR)
  Z <- 1.96  # For 95% confidence interval
  log_lower <- log_OR - Z * SE_log_OR
  log_upper <- log_OR + Z * SE_log_OR
  
  # Convert log CI back to the original scale
  CI_lower <- exp(log_lower)
  CI_upper <- exp(log_upper)
  
  # Initialize p-values
  fisher_p_value <- NA
  chi_square_p_value <- NA
  
  # Perform statistical tests
  contingency_table <- matrix(c(A, B, C, D), nrow = 2, byrow = TRUE)
  
  # # Fisher's Exact Test
  # try({
  #   fisher_result <- fisher.test(contingency_table)
  #   fisher_p_value <- fisher_result$p.value
  # }, silent = TRUE)
  # 
  # # Chi-Square Test
  # try({
  #   chi_square_result <- chisq.test(contingency_table)
  #   chi_square_p_value <- chi_square_result$p.value
  # }, silent = TRUE)
  
  data.table(
    OR = OR,
    Lower_CI = CI_lower,
    Upper_CI = CI_upper
    # fisher_p_value = fisher_p_value,
    # chi_square_p_value = chi_square_p_value
  )
}

test_enrichment <- function(dt_sub, v_cohort){
  
  # Subset by score and pathogenecity prediction
  dt_sub <- dt_sub[,c("chromosome", "pos", "ref", "alt", paste0(v_cohort, "_AC"), paste0(v_cohort, "_AN")), with = FALSE]
  
  # Custom function to count non-zero allele counts
  count_nonzero <- function(x) {
    x <- x[!is.na(x)]  # Remove NA values
    return(length(x[x > 0]))  # Count elements that are greater than zero
  }
  
  dt_mod <- data.table(cohort = v_cohort,
                       n_cohort = apply(dt[, paste0(v_cohort, "_AN"), with = FALSE], 2, max, na.rm = TRUE) / 2,
                       AC_unique = apply(dt_sub[, paste0(v_cohort, "_AC"), with = FALSE], 2, count_nonzero),
                       AN_unique = apply(dt_sub[, paste0(v_cohort, "_AN"), with = FALSE], 2, length),
                       AC = apply(dt_sub[, paste0(v_cohort, "_AC"), with = FALSE], 2, sum, na.rm = TRUE),
                       AN = apply(dt_sub[, paste0(v_cohort, "_AN"), with = FALSE], 2, sum, na.rm = TRUE)
  )
  
  # # Apply the function to each cohort
  # dt_result <- cbind(dt_mod, 
  #                    dt_mod[, compute_or_manual(AC_case = AC, 
  #                                               AN_case = AN,
  #                                               AC_control = dt_mod$AC[dt_mod$cohort == "gnomAD"],
  #                                               AN_control = dt_mod$AN[dt_mod$cohort == "gnomAD"])]
  # )
  
  # Apply the function to each cohort (treat AN as n_cohort so not to inflate confidence intervals)
  dt_result <- cbind(dt_mod, 
                     dt_mod[, compute_or_manual(AC_case = AC, 
                                                AN_case = n_cohort,
                                                AC_control = dt_mod$AC[dt_mod$cohort == "gnomAD"],
                                                AN_control = dt_mod$n_cohort[dt_mod$cohort == "gnomAD"])]
  )
  
  return(dt_result)
}

### Run analysis by gene ==========

list_dcm <- list()
for(i in seq_along(dcm_genes)){
  dt_sub_gene <- dt[gene_name == dcm_genes[i]]
  list_pathogenic <- list() 
  for(j in seq_along(v_scores)){
    list_pathogenic[[j]] <- test_enrichment(dt_sub = dt_sub_gene[get(paste0(v_scores[j], "_binary")) == "P", 
                                                                  with = TRUE], 
                                             v_cohort = c("gnomAD", v_cohort_dcm))
    list_pathogenic[[j]]$score <- v_scores[j]
    list_pathogenic[[j]]$prediction <- "pathogenic"
    
  }
  list_benign <- list() 
  for(j in seq_along(v_scores)){
    list_benign[[j]] <- test_enrichment(dt_sub = dt_sub_gene[get(paste0(v_scores[j], "_binary")) == "B", 
                                                                  with = TRUE], 
                                             v_cohort = c("gnomAD", v_cohort_dcm))
    list_benign[[j]]$score <- v_scores[j]
    list_benign[[j]]$prediction <- "benign"
    
  }
  list_all <- list() 
  for(j in seq_along(v_scores)){
    list_all[[j]] <- test_enrichment(dt_sub = dt_sub_gene[get(paste0(v_scores[j], "_binary")) == "VUS", 
                                                          with = TRUE], 
                                         v_cohort = c("gnomAD", v_cohort_dcm))
    list_all[[j]]$score <- v_scores[j]
    list_all[[j]]$prediction <- "VUS"
    
  }
  list_dcm[[i]] <- rbind(rbindlist(list_pathogenic), rbindlist(list_benign), rbindlist(list_all))
  list_dcm[[i]]$gene <- dcm_genes[i]
  print(dcm_genes[i])
}
list_hcm <- list()
for(i in seq_along(hcm_genes)){
  dt_sub_gene <- dt[gene_name == hcm_genes[i]]
  list_pathogenic <- list() 
  for(j in seq_along(v_scores)){
    list_pathogenic[[j]] <- test_enrichment(dt_sub = dt_sub_gene[get(paste0(v_scores[j], "_binary")) == "P", 
                                                                  with = TRUE], 
                                             v_cohort = c("gnomAD", v_cohort_hcm))
    list_pathogenic[[j]]$score <- v_scores[j]
    list_pathogenic[[j]]$prediction <- "pathogenic"
    
  }
  list_benign <- list() 
  for(j in seq_along(v_scores)){
    list_benign[[j]] <- test_enrichment(dt_sub = dt_sub_gene[get(paste0(v_scores[j], "_binary")) == "B", 
                                                              with = TRUE], 
                                         v_cohort = c("gnomAD", v_cohort_hcm))
    list_benign[[j]]$score <- v_scores[j]
    list_benign[[j]]$prediction <- "benign"
    
  }
  list_all <- list() 
  for(j in seq_along(v_scores)){
    list_all[[j]] <- test_enrichment(dt_sub = dt_sub_gene[get(paste0(v_scores[j], "_binary")) == "VUS", 
                                                          with = TRUE], 
                                      v_cohort = c("gnomAD", v_cohort_hcm))
    list_all[[j]]$score <- v_scores[j]
    list_all[[j]]$prediction <- "VUS"
    
  }
  list_hcm[[i]] <- rbind(rbindlist(list_pathogenic), rbindlist(list_benign), rbindlist(list_all))
  list_hcm[[i]]$gene <- hcm_genes[i]
  print(hcm_genes[i])
}
dt_gene <- unique(rbind(rbindlist(list_dcm), rbindlist(list_hcm)))
dt_gene$n_cohort[dt_gene$n_cohort == -Inf] <- 0
dt_gene$disease <- NA
dt_gene$disease[dt_gene$gene %in% hcm_genes] <- "HCM"
dt_gene$disease[dt_gene$gene %in% dcm_genes] <- "DCM"
dt_gene$disease[dt_gene$gene %in% hcm_genes & dt_gene$gene %in% dcm_genes] <- "HCM and DCM"

### Run analysis by disease and total ==========

# HCM genes
dt_sub_hcm <- dt[gene_name %in% hcm_genes]
list_pathogenic <- list() 
for(i in seq_along(v_scores)){
  list_pathogenic[[i]] <- test_enrichment(dt_sub =  dt_sub_hcm[get(paste0(v_scores[i], "_binary")) == "P", 
                                                                with = TRUE],
                                           v_cohort = c("gnomAD", v_cohort_hcm))
  list_pathogenic[[i]]$score <- v_scores[i]
  list_pathogenic[[i]]$prediction <- "pathogenic"
}
list_benign <- list() 
for(i in seq_along(v_scores)){
  list_benign[[i]] <- test_enrichment(dt_sub =  dt_sub_hcm[get(paste0(v_scores[i], "_binary")) == "B", 
                                                            with = TRUE],
                                       v_cohort = c("gnomAD", v_cohort_hcm))
  list_benign[[i]]$score <- v_scores[i]
  list_benign[[i]]$prediction <- "benign"
}
list_all <- list() 
for(i in seq_along(v_scores)){
  list_all[[i]] <- test_enrichment(dt_sub =  dt_sub_hcm[get(paste0(v_scores[i], "_binary")) == "VUS", 
                                                        with = TRUE],
                                       v_cohort = c("gnomAD", v_cohort_hcm))
  list_all[[i]]$score <- v_scores[i]
  list_all[[i]]$prediction <- "VUS"
}
dt_hcm <- rbind(rbindlist(list_pathogenic), rbindlist(list_benign), rbindlist(list_all))
dt_hcm$disease <- "HCM"

# DCM genes
dt_sub_dcm <- dt[gene_name %in% dcm_genes]
list_pathogenic <- list() 
for(i in seq_along(v_scores)){
  list_pathogenic[[i]] <- test_enrichment(dt_sub =  dt_sub_dcm[get(paste0(v_scores[i], "_binary")) == "P", 
                                                                with = TRUE],
                                           v_cohort = c("gnomAD", v_cohort_dcm))
  list_pathogenic[[i]]$score <- v_scores[i]
  list_pathogenic[[i]]$prediction <- "pathogenic"
}
list_benign <- list() 
for(i in seq_along(v_scores)){
  list_benign[[i]] <- test_enrichment(dt_sub =  dt_sub_dcm[get(paste0(v_scores[i], "_binary")) == "B", 
                                                            with = TRUE],
                                       v_cohort = c("gnomAD", v_cohort_dcm))
  list_benign[[i]]$score <- v_scores[i]
  list_benign[[i]]$prediction <- "benign"
}
list_all <- list() 
for(i in seq_along(v_scores)){
  list_all[[i]] <- test_enrichment(dt_sub =  dt_sub_dcm[get(paste0(v_scores[i], "_binary")) == "VUS", 
                                                        with = TRUE],
                                    v_cohort = c("gnomAD", v_cohort_dcm))
  list_all[[i]]$score <- v_scores[i]
  list_all[[i]]$prediction <- "VUS"
}
dt_dcm <- rbind(rbindlist(list_pathogenic), rbindlist(list_benign), rbindlist(list_all))
dt_dcm$disease <- "DCM"

# HCM and DCM genes
dt_sub <- dt[gene_name %in% c(hcm_genes, dcm_genes)]
list_pathogenic <- list() 
for(i in seq_along(v_scores)){
  list_pathogenic[[i]] <- test_enrichment(dt_sub = dt_sub[get(paste0(v_scores[i], "_binary")) == "P", 
                                                                with = TRUE],
                                           v_cohort = c("gnomAD", "ALL_HCM_DCM"))
  list_pathogenic[[i]]$score <- v_scores[i]
  list_pathogenic[[i]]$prediction <- "pathogenic"
}
list_benign <- list() 
for(i in seq_along(v_scores)){
  list_benign[[i]] <- test_enrichment(dt_sub =  dt_sub[get(paste0(v_scores[i], "_binary")) == "B", 
                                                            with = TRUE],
                                       v_cohort = c("gnomAD", "ALL_HCM_DCM"))
  list_benign[[i]]$score <- v_scores[i]
  list_benign[[i]]$prediction <- "benign"
}
list_all <- list() 
for(i in seq_along(v_scores)){
  list_all[[i]] <- test_enrichment(dt_sub =  dt_sub[get(paste0(v_scores[i], "_binary")) == "VUS", 
                                                    with = TRUE],
                                    v_cohort = c("gnomAD", "ALL_HCM_DCM"))
  list_all[[i]]$score <- v_scores[i]
  list_all[[i]]$prediction <- "VUS"
}
dt_hcm_dcm <- rbind(rbindlist(list_pathogenic), rbindlist(list_benign), rbindlist(list_all))
dt_hcm_dcm$disease <- "HCM and DCM"

dt_totals <- rbindlist(list(dt_hcm, dt_dcm, dt_hcm_dcm))
dt_totals <- dt_totals[cohort != "gnomAD"]

### Calculate OR for synonymous variants ==========

# HCM genes
dt_hcm_syn <- test_enrichment(dt_sub =  dt_syn[gene_name %in% hcm_genes],
                              v_cohort = c("gnomAD", v_cohort_hcm))
dt_hcm_syn$consequence <- "synonymous"
dt_hcm_syn$disease <- "HCM"

# DCM genes
dt_dcm_syn <- test_enrichment(dt_sub =  dt_syn[gene_name %in% dcm_genes],
                              v_cohort = c("gnomAD", v_cohort_dcm))
dt_dcm_syn$consequence <- "synonymous"
dt_dcm_syn$disease <- "DCM"

dt_hcm_dcm_syn <- rbindlist(list(dt_hcm_syn, dt_dcm_syn))

### Export ==========

dt_gene_out <- rbind(dt_gene[(disease == "HCM" | disease == "HCM and DCM") & cohort == "ALL_HCM"],
                     dt_gene[(disease == "DCM" | disease == "HCM and DCM") & cohort == "ALL_DCM"])

fwrite(dt_gene_out, "../results/OR_by_gene.csv")
fwrite(dt_totals, "../results/OR_by_cohort.csv")

########

# x <- dt_gene[(disease == "DCM" | disease == "HCM and DCM") & cohort %in% c("goDCM_EUR", "UKB_dcm_all")]
# 
# dt_dcm_m <- dt[gene_name %in% dcm_genes,c("gene_name", "chromosome", "pos", "ref", "alt",
#                                                                       "goDCM_EUR_AC", "goDCM_EUR_AN", "VFX_ALL_DCM_AC","VFX_ALL_DCM_AN",
#                                                                      "UKB_dcm_all_AC", "UKB_dcm_all_AN", "gnomAD_AC","gnomAD_AN")]
# 
# max(dt_dcm_m$goDCM_EUR_AC, na.rm = TRUE)
# max(dt_dcm_m$gnomAD_AC, na.rm = TRUE)
# 
# sum(dt_dcm_m$goDCM_EUR_AC, na.rm = TRUE) / sum(dt_dcm_m$goDCM_EUR_AN, na.rm = TRUE)
# sum(dt_dcm_m$gnomAD_AC, na.rm = TRUE) / sum(dt_dcm_m$gnomAD_AN, na.rm = TRUE)
# sum(dt_dcm_m$UKB_dcm_all_AC, na.rm = TRUE) / sum(dt_dcm_m$UKB_dcm_all_AN, na.rm = TRUE)
# sum(dt_dcm_m$VFX_ALL_DCM_AC, na.rm = TRUE) / sum(dt_dcm_m$VFX_ALL_DCM_AN, na.rm = TRUE)
# 
# dt_dcm_s <- dt_syn[gene_name %in% dcm_genes,c("gene_name", "chromosome", "pos", "ref", "alt",
#                                           "goDCM_EUR_AC", "goDCM_EUR_AN", "VFX_ALL_DCM_AC","VFX_ALL_DCM_AN",
#                                           "UKB_dcm_all_AC", "UKB_dcm_all_AN", "gnomAD_AC","gnomAD_AN")]
# 
# max(dt_dcm_s$goDCM_EUR_AC, na.rm = TRUE)
# max(dt_dcm_s$gnomAD_AC, na.rm = TRUE)
# 
# sum(dt_dcm_s$goDCM_EUR_AC, na.rm = TRUE) / sum(dt_dcm_s$goDCM_EUR_AN, na.rm = TRUE)
# sum(dt_dcm_s$gnomAD_AC, na.rm = TRUE) / sum(dt_dcm_s$gnomAD_AN, na.rm = TRUE)
# sum(dt_dcm_s$UKB_dcm_all_AC, na.rm = TRUE) / sum(dt_dcm_s$UKB_dcm_all_AN, na.rm = TRUE)
# sum(dt_dcm_s$VFX_ALL_DCM_AC, na.rm = TRUE) / sum(dt_dcm_s$VFX_ALL_DCM_AN, na.rm = TRUE)


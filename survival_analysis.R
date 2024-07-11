rm(list=ls())
graphics.off()

library(data.table)
library(ggpubr)
library(survival)
library(survminer)
library(ggplot2)
library(gridExtra)
library(forestplot)
library(RColorBrewer)

### Functions ==========

## Function that fits a Cox proportional hazards model and returns the coefficients, HRs, CIs, and p-values
# Input:  
# dt_mod: time, event, x, covariates 
# x: predictor of interest column name 
# vec_x: order of x factors (first factor will be model reference).  
# covars: covariate column names 
# Output: 
# table of HRs and P-values 

# x = "genetic_status"
# vec_x = vec_x_vneg
# covars = c("Sex")

cox_fit <- function(dt_mod, x, vec_x, covars){
  
  # QC
  dt_mod <- dt_mod[!is.na(event) & !is.na(time) & time >= 0,]
  
  # Set reference level
  dt_mod <- dt_mod[get(x) %in% vec_x]
  dt_mod[, x := factor(get(x), levels = vec_x)]
  
  # Create the survival object
  surv_obj <- Surv(time = dt_mod$time, event = dt_mod$event)
  
  # Create formula for cox model 
  formula <- as.formula(paste("surv_obj ~", paste(c("x", covars), collapse = " + ")))
  
  # Fit Cox proportional hazards model 
  cox_model <- coxph(formula, 
                     data = dt_mod, 
                     control = coxph.control(iter.max = 1000))
  
  # Extract the summary of the fitted Cox model
  summary_cox <- summary(cox_model)
  summary_cox
  
  # Extract the coefficients, standard errors, and p-values
  dt_cox_coef <- as.data.table(summary_cox$coefficients, keep.rownames = "covariate")
  
  # Calculate hazard ratios and confidence intervals
  dt_cox_coef[, `:=`(
    covariate_string = covariate,
    HR = exp(coef),
    CI95_lower = exp(coef - 1.96 * `se(coef)`),
    CI95_upper = exp(coef + 1.96 * `se(coef)`),
    p_val = `Pr(>|z|)`
  )]
  dt_cox_coef$`Pr(>|z|)` <- NULL
  dt_cox_coef$covariate <- NULL
  dt_cox_coef
  
  # Split the covariates and the factors into seperate columns
  # create lookup table 
  list_covar <- list()  
  list_factor <- list()  
  list_reference <- list()
  for(i in seq_along(covars)){
    if(is.factor(dt_mod[[covars[i]]]) == TRUE){
      list_covar[[i]] <- rep(covars[i], length(unique(dt_mod[[covars[i]]]))-1)
      list_factor[[i]] <- levels(dt_mod[[covars[i]]])[2:length(levels(dt_mod[[covars[i]]]))]
      list_reference[[i]] <- levels(dt_mod[[covars[i]]])[1]
    }
    if(is.factor(dt_mod[[covars[i]]]) == FALSE){
      list_covar[[i]] <- covars[i]
      list_factor[[i]] <- NA
      list_reference[[i]] <- NA
    }
  }
  dt_lookup <- data.table(covariate = c(rep("x", length(levels(dt_mod$x))-1), unlist(list_covar)),
                          reference = c(rep(levels(dt_mod[["x"]])[1], length(levels(dt_mod$x))-1), unlist(list_reference)),
                          factor = c(levels(dt_mod[["x"]])[2:length(levels(dt_mod[["x"]]))], unlist(list_factor))
  )
  dt_lookup$covariate_string <- NA
  dt_lookup$covariate_string[is.na(dt_lookup$factor)] <- dt_lookup$covariate[is.na(dt_lookup$factor)]
  dt_lookup$covariate_string[!is.na(dt_lookup$factor)] <- paste0(dt_lookup$covariate[!is.na(dt_lookup$factor)], dt_lookup$factor[!is.na(dt_lookup$factor)])
  dt_cox_coef <- merge(dt_lookup, dt_cox_coef, by = "covariate_string", all.y = TRUE)
  dt_cox_coef$covariate_string <- NULL
  dt_cox_coef$covariate[dt_cox_coef$covariate == "x"] <- x
  
  # Add n factor and n reference
  list_n <- list()
  for(i in seq_along(covars)){ list_n[[i]] <- as.data.table(table(dt_mod[[covars[i]]])) }
  dt_n <- rbindlist(list_n)
  dt_n <- rbindlist(list(dt_n, as.data.table(table(dt_mod$x))))
  colnames(dt_n) <- c("factor", "n_factor")
  dt_cox_coef <- merge(dt_cox_coef, dt_n, by = "factor", all.x = TRUE)
  colnames(dt_n) <- c("reference", "n_reference")
  dt_cox_coef <- merge(dt_cox_coef, dt_n, by = "reference", all.x = TRUE)
  
  return(dt_cox_coef)
}

## Function that fits a Kaplan-Meier model and returns the plot
# Input:  
# dt_mod: time, event, x, covariates 
# x: predictor of interest column name 
# vec_x: order of x factors (first factor will be model reference).  
# covars: covariate column names 
# Output: 
# table of HRs and P-values 

# dt_mod = dt_mod
# x = "Sarc_status"
# vec_x = c("sarc(-)", "sarc(+)")
# x_label = "Age"
# y_label = "Proportion free of endpoint"
# plot_title = dt_outcomes_birth$outcome[i]
# p_val_text = paste0("p = ", format(dt_cox$p_val[dt_cox$factor == "sarc(+)"], digits = 3))
# vec_colors = c(
#   "sarc(+)" = "grey40",
#   "sarc(-)" = "grey80"
# )

km_fit <- function(dt_mod, x, vec_x, x_label, y_label, plot_title, p_val_text, vec_colors){
  
  # QC
  dt_mod <- dt_mod[!is.na(event) & !is.na(time),]
  
  # Set reference level
  dt_mod <- dt_mod[get(x) %in% vec_x]
  dt_mod[, x := factor(get(x), levels = vec_x)]
  
  # Create the survival object
  surv_obj <- Surv(time = dt_mod$time, event = dt_mod$event)
  
  # Fit km model
  # km_fit <- survfit(formula, data = dt_mod)
  km_fit <- survfit(Surv(time = dt_mod$time, event = dt_mod$event) ~ x, data = dt_mod)
  summary(km_fit)
  
  # Remove "Group=" from strata
  names(km_fit$strata) <- gsub(paste0("x", "="), "", names(km_fit$strata))
  
  # Plot Kaplan-Meier curves
  km_plot <- ggsurvplot(
    km_fit,
    data = dt_mod,
    pval = p_val_text,
    pval.method = FALSE,
    conf.int = FALSE,
    risk.table = TRUE,
    palette = vec_colors,
    title = plot_title,
    xlab = x_label, 
    ylab = y_label
  )
  km_plot
  
  # Adjust the legend size and remove the legend title
  km_plot$plot <- km_plot$plot + 
    theme(
      legend.text = element_text(size = 10),  # Adjust the size as needed
      legend.title = element_blank(), # Remove the legend title
      legend.position = "right",
    )
  km_plot
  
  return(km_plot)
}

## Function that plots a forest plot of hazard ratios
# Input:
# dt_forest = dt_cox_diagnosis

# dt_forest <- dt_plot[order(n_factor)]
# dt_forest[, factor := factor(factor, levels = rev(dt_forest$factor))]
# table_text <- cbind(c("Genotype", as.character(dt_forest$factor)),
#                     c("n Genotype", dt_forest$n_factor),
#                     c("Hazard Ratio", round(dt_forest$HR, 2)),
#                     c("P-value", format(dt_forest$p_val_fdr, digits = 2)))
# plot_title <- dt_forest$y_var[1]

plot_HR <- function(dt_forest, table_text, plot_title){
  
  # Plot forest plot
  forest_plot <- forestplot(
    labeltext = table_text,
    mean = c(NA, dt_forest$HR),
    lower = c(NA, dt_forest$CI95_lower),
    upper = c(NA, dt_forest$CI95_upper),
    is.summary = c(TRUE, rep(FALSE, nrow(dt_forest))),
    zero = 1,
    xlog = TRUE,
    col = forestplot::fpColors(box = "blue", lines = "black", zero = "red"),
    xlab = "Hazard Ratio (log scale)",
    ci.vertices = TRUE,
    ci.vertices.height = 0.1,
    title = plot_title
  )
  forest_plot
  
  return(forest_plot)
}


### Import ==========

dt <- fread("../data/analysis_table.csv.gz")
dt_share <- fread("../data/SHaRe_analysis_table.csv")


### Set variables ==========

# HCM genes with evidence of missense pathogenicity
hcm_genes <- c("ACTC1", "MYBPC3", "MYH7", "MYL2", "MYL3", "PLN", "TNNI3", "TNNT2", "TPM1")

# Outcomes
dt_outcomes_birth <- data.table(col_name = c("Diagnosis", "Arrhythmia_AFib", "Composite_HF", "Composite_VArrhythmia", "Composite_Overall"),
                                outcome = c("Diagnosis Age", "Artial Fibrilation", "Heart Failure", "Ventricular Arrhythmia", "Overall Composite"))

dt_outcomes_pd <- data.table(col_name = c("Arrhythmia_AFib", "Composite_HF", "Composite_VArrhythmia", "Composite_Overall"),
                             outcome = c("Artial Fibrilation", "Heart Failure", "Ventricular Arrhythmia", "Overall Composite"))


### Prepocessing ==========

# Set score binary predictions based on recommended thresholds
dt$AlphaMissense_binary <- NA
dt$AlphaMissense_binary[dt$AlphaMissense_class == "likely_pathogenic"] <- 1
dt$AlphaMissense_binary[dt$AlphaMissense_class == "likely_benign"] <- 0

dt$CardioBoost_binary <- NA
dt$CardioBoost_binary[dt$CardioBoost_pathogenicity >= 0.1] <- 1
dt$CardioBoost_binary[dt$CardioBoost_pathogenicity <= 0.9] <- 0

dt$REVEL_binary <- NA
dt$REVEL_binary[dt$REVEL >= 0.75] <- 1
dt$REVEL_binary[dt$REVEL <= 0.25] <- 0

dt$CADD_binary <- NA
dt$CADD_binary[dt$CADD_PHRED >= 20] <- 1
dt$CADD_binary[dt$CADD_PHRED <= 10] <- 0

# Merge
dt$VariantID <- paste0(gsub("chr", "", dt$chromosome), "-", dt$pos, "-", dt$ref, "-", dt$alt)
dt <- merge(dt_share, dt[, c("VariantID",
                             "AlphaMissense_binary",
                             "CardioBoost_binary",  
                             "REVEL_binary",
                             "CADD_binary")], by="VariantID", all.x=TRUE)

# Filter for missense 
dt <- dt[VEP_consequence1 == "missense_variant" | Sarc_status == "sarc(-)", ]

# Add age at last record to age_ columns with NA
dt$age_Arrhythmia_AFib[is.na(dt$age_Arrhythmia_AFib)] <- dt$LastEncounter_Age[is.na(dt$age_Arrhythmia_AFib)]
dt$age_Composite_HF[is.na(dt$age_Composite_HF)] <- dt$LastEncounter_Age[is.na(dt$age_Composite_HF)]
dt$age_Composite_VArrhythmia[is.na(dt$age_Composite_VArrhythmia)] <- dt$LastEncounter_Age[is.na(dt$age_Composite_VArrhythmia)]
dt$age_Composite_Overall[is.na(dt$age_Composite_Overall)] <- dt$LastEncounter_Age[is.na(dt$age_Composite_Overall)]

# Set reference level 
dt <- dt[Sex %in% c("Female", "Male")]
dt[, Sex := factor(Sex, levels = c("Female", "Male"))]

# Format diagnosis age outcome
dt$event_Diagnosis <- 1
dt$event_Diagnosis[is.na(dt$Primary_Diagnosis_Age)] <- NA
dt$age_Diagnosis <- dt$Primary_Diagnosis_Age

# Subset VUS and format factors
dt_vus <- dt[Sarc_status == "sarc(VUS)", ]
vec_AM <- rep(NA, nrow(dt_vus))
vec_AM[dt_vus$AlphaMissense_binary == 1] <- "AlphaMissense(P)"
vec_AM[dt_vus$AlphaMissense_binary == 0] <- "AlphaMissense(B)"

vec_CB <- rep(NA, nrow(dt_vus))
vec_CB[dt_vus$CardioBoost_binary == 1] <- "CardioBoost(P)"
vec_CB[dt_vus$CardioBoost_binary == 0] <- "CardioBoost(B)"

vec_RV <- rep(NA, nrow(dt_vus))
vec_RV[dt_vus$REVEL_binary == 1] <- "REVEL(P)"
vec_RV[dt_vus$REVEL_binary == 0] <- "REVEL(B)"

vec_CD <- rep(NA, nrow(dt_vus))
vec_CD[dt_vus$CADD_binary == 1] <- "CADD(P)"
vec_CD[dt_vus$CADD_binary == 0] <- "CADD(B)"

dt_vus <- rbind(dt_vus, dt_vus, dt_vus, dt_vus)
dt_vus$x_var <- c(vec_AM, vec_CB, vec_RV, vec_CD)

dt$x_var <- dt$Sarc_status
dt <- rbind(dt, dt_vus)

### Set colours for variant status ==========

# Step 1: Define base colors for known strata
base_colors <- c(
  "sarc(+)" = "grey40",
  "sarc(-)" = "grey80",
  "sarc(VUS)" = "red"
)

# Generate a large enough color palette
color_palette <- rev(colorRampPalette(brewer.pal(8, "Paired"))(8))

# Create a named vector for variant colors
colors_named <- setNames(color_palette, c("AlphaMissense(P)", "AlphaMissense(B)", "CardioBoost(P)", "CardioBoost(B)", "REVEL(P)", "REVEL(B)", "CADD(P)", "CADD(B)"))

# Step 3: Combine the base colors and variant colors
strata_colors <- c(base_colors, colors_named)


### Model fitting ==========

# Fit Cox model and record outcomes: 
# sarc(+) and sarc(VUS) vs sarc(-)
# AlphaMissense(P) vs AlphaMissense(B)
# CardioBoost(P) vs CardioBoost(B)
# REVEL(P) vs REVEL(B)
# CADD(P) vs CADD(B)

# Plot KM:
# sarc(+) sarc(VUS) sarc(-) AlphaMissense(P) vs AlphaMissense(B)
# sarc(+) sarc(VUS) sarc(-) CardioBoost(P) vs CardioBoost(B)
# sarc(+) sarc(VUS) sarc(-) REVEL(P) vs REVEL(B)
# sarc(+) sarc(VUS) sarc(-) CADD(P) vs CADD(B)

# Run for Sarc(+) vs Sarc(-)

dt_mod <- dt[x_var %in% c("sarc(-)", "sarc(+)", "sarc(VUS)")]
dt_mod <- dt_mod[,c(paste0(c("event_", "age_"), "Diagnosis"), "Sex", "x_var"), with = FALSE]
colnames(dt_mod) <- c("event", "time", "Sex", "x_var")
dt_mod <- dt_mod[time != 0 & !is.na(time)]
dt_cox_sarc <- cox_fit(dt_mod = dt_mod, 
                       x = "x_var",
                       vec_x = c("sarc(-)", "sarc(VUS)", "sarc(+)"),
                       covars = c("Sex"))
dt_cox_sarc
dt_cox_sarc$y_var <- "Diagnosis Age"


# Run for each score on VUS

vec_scores <- c("AlphaMissense", "CardioBoost", "REVEL", "CADD")
list_cox <- list()
list_km <- list()
for (i in seq_along(vec_scores)){
  
  # Fit cox model
  dt_mod <- dt[x_var %in% c(paste0(vec_scores[i], "(B)"), 
                            paste0(vec_scores[i], "(P)"))]
  dt_mod <- dt_mod[,c(paste0(c("event_", "age_"), "Diagnosis"), "Sex", "x_var"), with = FALSE]
  colnames(dt_mod) <- c("event", "time", "Sex", "x_var")
  dt_mod <- dt_mod[time != 0 & !is.na(time)]
  dt_cox <- cox_fit(dt_mod = dt_mod, 
                    x = "x_var",
                    vec_x = c(paste0(vec_scores[i], "(B)"), 
                              paste0(vec_scores[i], "(P)")),
                    covars = c("Sex"))
  dt_cox
  dt_cox$y_var <- "Diagnosis Age"
  
  list_cox[[i]] <- dt_cox
  
  # Fit km model
  dt_mod <- dt[x_var %in% c("sarc(-)", "sarc(VUS)", "sarc(+)", paste0(vec_scores[i], "(B)"), 
                            paste0(vec_scores[i], "(P)"))]
  dt_mod <- dt_mod[,c(paste0(c("event_", "age_"), "Diagnosis"), "Sex", "x_var"), with = FALSE]
  colnames(dt_mod) <- c("event", "time", "Sex", "x_var")
  dt_mod <- dt_mod[time != 0 & !is.na(time)]
  km_fit <- survfit(Surv(time = dt_mod$time, event = dt_mod$event) ~ x_var, data = dt_mod)
  summary(km_fit)
  
  # Remove "Group=" from strata
  names(km_fit$strata) <- gsub(paste0("x_var", "="), "", names(km_fit$strata))
  
  # Plot Kaplan-Meier curves
  km_plot <- ggsurvplot(
    km_fit,
    data = dt_mod,
    pval = "",
    pval.method = FALSE,
    conf.int = FALSE,
    risk.table = TRUE,
    palette = strata_colors,
    title = "Diagnosis Age",
    xlab = "Age", 
    ylab = "Proportion free of endpoint"
  )
  km_plot
  
  # Adjust the legend size and remove the legend title
  km_plot$plot <- km_plot$plot + 
    theme(
      legend.text = element_text(size = 10),  # Adjust the size as needed
      legend.title = element_blank(), # Remove the legend title
      legend.position = "right",
    )
  km_plot
  
  list_km[[i]] <- km_plot$plot
}

dt_cox <- rbindlist(list_cox)
dt_cox <- rbind(dt_cox_sarc, dt_cox)
dt_cox <- dt_cox[covariate != "Sex"]

# Plot forest plot

dt_forest <- dt_cox
dt_forest[, factor := factor(factor, levels = (dt_forest$factor))]
table_text <- cbind(c("Group 1", as.character(dt_forest$reference)),
                    c("Group 2", as.character(dt_forest$factor)),
                    c("n Group 1", dt_forest$n_reference),
                    c("n Group 2", dt_forest$n_factor),
                    c("Hazard Ratio", round(dt_forest$HR, 2)),
                    c("P-value", format(dt_forest$p_val, digits = 2)))
plot_title <- "Diagnosis Age"
plot_hr <- plot_HR(dt_forest = dt_forest, table_text = table_text, plot_title = plot_title)

jpeg(paste0("../results/Figure_diagnosis_age_HR_from_birth.jpg"), 
     width = 14, height = 4, units = "in", res = 100)
print(plot_hr)
dev.off()

# Plot KM
p_km_panel <- ggarrange(list_km[[1]], list_km[[2]], list_km[[3]], list_km[[4]], ncol = 2, nrow = 2)
ggsave(p_km_panel, filename = "../results/Figure_KM_diagnosis_age_from_birth.jpg", width = 9, height = 7, units = "in")


### Run for composite overall 

# Run for Sarc(+) vs Sarc(-)

dt_mod <- dt[x_var %in% c("sarc(-)", "sarc(+)", "sarc(VUS)")]
dt_mod <- dt_mod[,c(paste0(c("event_", "age_"), "Composite_Overall"), "Sex", "x_var"), with = FALSE]
colnames(dt_mod) <- c("event", "time", "Sex", "x_var")
dt_mod <- dt_mod[time != 0 & !is.na(time)]
dt_cox_sarc <- cox_fit(dt_mod = dt_mod, 
                       x = "x_var",
                       vec_x = c("sarc(-)", "sarc(VUS)", "sarc(+)"),
                       covars = c("Sex"))
dt_cox_sarc
dt_cox_sarc$y_var <- "Composite Overall"


# Run for each score on VUS

vec_scores <- c("AlphaMissense", "CardioBoost", "REVEL", "CADD")
list_cox <- list()
list_km <- list()
for (i in seq_along(vec_scores)){
  
  # Fit cox model
  dt_mod <- dt[x_var %in% c(paste0(vec_scores[i], "(B)"), 
                            paste0(vec_scores[i], "(P)"))]
  dt_mod <- dt_mod[,c(paste0(c("event_", "age_"), "Composite_Overall"), "Sex", "x_var"), with = FALSE]
  colnames(dt_mod) <- c("event", "time", "Sex", "x_var")
  dt_mod <- dt_mod[time != 0 & !is.na(time)]
  dt_cox <- cox_fit(dt_mod = dt_mod, 
                    x = "x_var",
                    vec_x = c(paste0(vec_scores[i], "(B)"), 
                              paste0(vec_scores[i], "(P)")),
                    covars = c("Sex"))
  dt_cox
  dt_cox$y_var <- "Composite Overall"
  
  list_cox[[i]] <- dt_cox
  
  # Fit km model
  dt_mod <- dt[x_var %in% c("sarc(-)", "sarc(VUS)", "sarc(+)", paste0(vec_scores[i], "(B)"), 
                            paste0(vec_scores[i], "(P)"))]
  dt_mod <- dt_mod[,c(paste0(c("event_", "age_"), "Composite_Overall"), "Sex", "x_var"), with = FALSE]
  colnames(dt_mod) <- c("event", "time", "Sex", "x_var")
  dt_mod <- dt_mod[time != 0 & !is.na(time)]
  km_fit <- survfit(Surv(time = dt_mod$time, event = dt_mod$event) ~ x_var, data = dt_mod)
  summary(km_fit)
  
  # Remove "Group=" from strata
  names(km_fit$strata) <- gsub(paste0("x_var", "="), "", names(km_fit$strata))
  
  # Plot Kaplan-Meier curves
  km_plot <- ggsurvplot(
    km_fit,
    data = dt_mod,
    pval = "",
    pval.method = FALSE,
    conf.int = FALSE,
    risk.table = TRUE,
    palette = strata_colors,
    title = "Composite Overall",
    xlab = "Age", 
    ylab = "Proportion free of endpoint"
  )
  km_plot
  
  # Adjust the legend size and remove the legend title
  km_plot$plot <- km_plot$plot + 
    theme(
      legend.text = element_text(size = 10),  # Adjust the size as needed
      legend.title = element_blank(), # Remove the legend title
      legend.position = "right",
    )
  km_plot
  
  list_km[[i]] <- km_plot$plot
}

dt_cox <- rbindlist(list_cox)
dt_cox <- rbind(dt_cox_sarc, dt_cox)
dt_cox <- dt_cox[covariate != "Sex"]

# Plot forest plot

dt_forest <- dt_cox
dt_forest[, factor := factor(factor, levels = (dt_forest$factor))]
table_text <- cbind(c("Group 1", as.character(dt_forest$reference)),
                    c("Group 2", as.character(dt_forest$factor)),
                    c("n Group 1", dt_forest$n_reference),
                    c("n Group 2", dt_forest$n_factor),
                    c("Hazard Ratio", round(dt_forest$HR, 2)),
                    c("P-value", format(dt_forest$p_val, digits = 2)))
plot_title <- "Composite Overall"
plot_hr <- plot_HR(dt_forest = dt_forest, table_text = table_text, plot_title = plot_title)

jpeg(paste0("../results/Figure_composite_overall_HR_from_birth.jpg"), 
     width = 14, height = 4, units = "in", res = 100)
print(plot_hr)
dev.off()

# Plot KM
p_km_panel <- ggarrange(list_km[[1]], list_km[[2]], list_km[[3]], list_km[[4]], ncol = 2, nrow = 2)
ggsave(p_km_panel, filename = "../results/Figure_KM_composite_overall_from_birth.jpg", width = 12, height = 9, units = "in")


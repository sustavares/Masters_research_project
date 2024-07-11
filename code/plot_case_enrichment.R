# Plot pathogenic and benign OR for each score in each cohort (and combined cohorts). 
# Figure 1 -- Three panel plot (HCM, DCM, HCM and DCM), Score on y axis, and OR on X axis. Colored by pathogenic or benign prediction
# Figure 2 -- Same plot as Figure 1 but stacked panels for each cohort (combined genes)
# Figure 3 -- Same plot as Figure 1 but stacked panels for each gene (combined cohort).

rm(list = ls())
graphics.off()


library(data.table)
library(epitools)
library(ggplot2)
library(ggpubr)

### Import ==========

dt_cohort <- fread("../results/OR_by_cohort.csv")
dt_gene <- fread("../results/OR_by_gene.csv")

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

### QC ==========

dt_cohort <- dt_cohort[AC >= 3]
dt_gene <- dt_gene[AC >= 3]

### Format ==========

dt_cohort$logOR <- log(dt_cohort$OR)
dt_cohort$logOR_lower_CI <- log(dt_cohort$Lower_CI)
dt_cohort$logOR_upper_CI <- log(dt_cohort$Upper_CI)

dt_gene$logOR <- log(dt_gene$OR)
dt_gene$logOR_lower_CI <- log(dt_gene$Lower_CI)
dt_gene$logOR_upper_CI <- log(dt_gene$Upper_CI)

# Vectors defining the mapping
v1 <- c("UKB_hcm_all", "VFX_BLF_HCM", "VFX_EGY_HCM", "VFX_GDX_HCM", "VFX_LMM_HCM", "VFX_OMG_HCM",
        "VFX_RBH_HCM", "VFX_SNG_HCM", "ALL_HCM", "UKB_dcm_all", "goDCM_EUR", "VFX_EGY_DCM",
        "VFX_LMM_DCM", "VFX_OMG_DCM", "VFX_RBH_DCM", "VFX_SNG_DCM", "ALL_DCM", "ALL_HCM_DCM")
v2 <- c("UKB", "BLF", "EGY", "GDX", "LMM", "OMG",
        "RBH", "SNG", "ALL", "UKB", "goDCM", "EGY",
        "LMM", "OMG", "RBH", "SNG", "ALL", "ALL")

# Create a named vector for mapping
labels_map <- setNames(v2, v1)

# Map the cohorts to new labels
dt_cohort[, label := labels_map[cohort]]

### Plot totals ==========

dt_plot_hcm <- dt_cohort[disease == "HCM" & cohort == "ALL_HCM"]
dt_plot_dcm <- dt_cohort[disease == "DCM" & cohort == "ALL_DCM"]
dt_plot_hdcm <- dt_cohort[disease == "HCM and DCM" & cohort == "ALL_HCM_DCM"]

custom_palette <- c("#1b9e77", "#d95f02", "black")
x_max <- max(c(dt_plot_dcm$Upper_CI, dt_plot_hcm$Upper_CI, dt_plot_hdcm$Upper_CI)) + 0.1
x_min <- min(c(dt_plot_dcm$Lower_CI, dt_plot_hcm$Lower_CI, dt_plot_hdcm$Lower_CI)) / 2

p_hcm <- ggplot(dt_plot_hcm, aes(x = OR, y = score, color = prediction)) +
  geom_point(size = 3, position = position_dodge(width = 0.5)) +
  geom_errorbarh(aes(xmin = Lower_CI, xmax = Upper_CI), height = 0.2, position = position_dodge(width = 0.5)) +
  # scale_x_continuous(trans = "log10") +  # Log-transform the x-axis
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +  # Add a dashed line where OR = 1
  scale_color_manual(values = custom_palette) +  # Apply the custom color palette
  scale_x_continuous(limits = c(0, x_max)) +
  labs(
    title = "A) HCM",
    x = "Odds Ratio",
    y = "",
    color = ""
  ) +
  theme_minimal() +
  theme(
    legend.position = "top",
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12, angle = 0),
    axis.line.x = element_line()
  )
p_hcm

p_dcm <- ggplot(dt_plot_dcm, aes(x = OR, y = score, color = prediction)) +
  geom_point(size = 3, position = position_dodge(width = 0.5)) +
  geom_errorbarh(aes(xmin = Lower_CI, xmax = Upper_CI), height = 0.2, position = position_dodge(width = 0.5)) +
  # scale_x_continuous(trans = "log10") +  # Log-transform the x-axis
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +  # Add a dashed line where OR = 1
  scale_color_manual(values = custom_palette) +  # Apply the custom color palette
  scale_x_continuous(limits = c(0, x_max)) +
  labs(
    title = "B) DCM",
    x = "Odds Ratio",
    y = "",
    color = ""
  ) +
  theme_minimal() +
  theme(
    legend.position = "top",
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12, angle = 0),
    axis.line.x = element_line()
  )
p_dcm

p_hcm_dcm <- ggplot(dt_plot_hdcm, aes(x = OR, y = score, color = prediction)) +
  geom_point(size = 3, position = position_dodge(width = 0.5)) +
  geom_errorbarh(aes(xmin = Lower_CI, xmax = Upper_CI), height = 0.2, position = position_dodge(width = 0.5)) +
  # scale_x_continuous(trans = "log10") +  # Log-transform the x-axis
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +  # Add a dashed line where OR = 1
  scale_color_manual(values = custom_palette) +  # Apply the custom color palette
  scale_x_continuous(limits = c(0, x_max)) +
  labs(
    title = "C) HCM and DCM",
    x = "Odds Ratio",
    y = "",
    color = ""
  ) +
  theme_minimal() +
  theme(
    legend.position = "top",
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12, angle = 0),
    axis.line.x = element_line()
  )
p_hcm_dcm

# Arrange plots
p_totals <- ggarrange(p_hcm, p_dcm, p_hcm_dcm,
                        ncol = 3, nrow = 1)
p_totals
ggsave("../results/Figure_case_enrichment_totals.png", p_totals, height = 4, width = 12, bg = "white")

### Plot by cohort for combined genes ==========

dt_plot_hcm <- dt_cohort[disease == "HCM" & cohort %in% v_cohort_hcm[v_cohort_hcm != "ALL_HCM"]]
dt_plot_dcm <- dt_cohort[disease == "DCM" & cohort %in% v_cohort_dcm[v_cohort_dcm != "ALL_DCM"]]

x_max <- max(c(dt_plot_dcm$Upper_CI, dt_plot_hcm$Upper_CI)) + 0.1
x_min <- min(c(dt_plot_dcm$Lower_CI, dt_plot_hcm$Lower_CI)) / 2
# x_min <- 0.1

list_hcm <- list()
for(i in seq_along(unique(dt_plot_hcm$cohort))){
  dt_plot <- dt_plot_hcm[cohort == unique(dt_plot_hcm$cohort)[i]]
  custom_palette <- c("#1b9e77", "#d95f02", "black")
  p_out <- ggplot(dt_plot, aes(x = OR, y = score, color = prediction)) +
    geom_point(size = 3, position = position_dodge(width = 0.5)) +
    geom_errorbarh(aes(xmin = Lower_CI, xmax = Upper_CI), height = 0.2, position = position_dodge(width = 0.5)) +
    scale_x_continuous(trans = "log10", limits = c(x_min, x_max)) +  # Log-transform the x-axis
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +  # Add a dashed line where OR = 1
    scale_color_manual(values = custom_palette) +  # Apply the custom color palette
    # scale_x_continuous(limits = c(0, x_max)) +
    labs(
      title = paste0(unique(dt_plot_hcm$label)[i], " (n = ", dt_plot$n_cohort[1], ")"),
      x = "log(Odds Ratio)",
      y = "",
      color = ""
    ) +
    theme_minimal() +
    theme(
      legend.position = "top",
      axis.text.y = element_text(size = 12),
      axis.text.x = element_text(size = 12, angle = 0),
      axis.line.x = element_line()
    )
  p_out
  list_hcm[[i]] <- p_out
}

list_dcm <- list()
for(i in seq_along(unique(dt_plot_dcm$cohort))){
  dt_plot <- dt_plot_dcm[cohort == unique(dt_plot_dcm$cohort)[i]]
  custom_palette <- c("#1b9e77", "#d95f02", "black")
  p_out <- ggplot(dt_plot, aes(x = OR, y = score, color = prediction)) +
    geom_point(size = 3, position = position_dodge(width = 0.5)) +
    geom_errorbarh(aes(xmin = Lower_CI, xmax = Upper_CI), height = 0.2, position = position_dodge(width = 0.5)) +
    # scale_x_continuous(trans = "log10") +  # Log-transform the x-axis
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +  # Add a dashed line where OR = 1
    scale_color_manual(values = custom_palette) +  # Apply the custom color palette
    scale_x_continuous(trans = "log10", limits = c(x_min, x_max)) +
    labs(
      title = paste0(unique(dt_plot_dcm$label)[i], " (n = ", dt_plot$n_cohort[1], ")"),
      x = "log(Odds Ratio)",
      y = "",
      color = ""
    ) +
    theme_minimal() +
    theme(
      legend.position = "top",
      axis.text.y = element_text(size = 12),
      axis.text.x = element_text(size = 12, angle = 0),
      axis.line.x = element_line()
    )
  p_out
  list_dcm[[i]] <- p_out
}

# Arrange and annotate
p_cohort_hcm <- ggarrange(plotlist = list_hcm, ncol = 2, nrow = ceiling(length(list_hcm)/2))
p_cohort_hcm
p_cohort_dcm <- ggarrange(plotlist = list_dcm, ncol = 2, nrow = ceiling(length(list_dcm)/2))
p_cohort_dcm

p_cohort_hcm <- annotate_figure(p_cohort_hcm, top = text_grob("HCM cohorts", size = 14, face = "bold"))
p_cohort_dcm <- annotate_figure(p_cohort_dcm, top = text_grob("DCM cohorts", size = 14, face = "bold"))

ggsave("../results/Figure_case_enrichment_by_cohort_HCM.png", p_cohort_hcm, height = 3*ceiling(length(list_hcm)/2), width = 10, bg = "white")
ggsave("../results/Figure_case_enrichment_by_cohort_DCM.png", p_cohort_dcm, height = 3*ceiling(length(list_dcm)/2), width = 10, bg = "white")

### Plot by gene ==========

dt_plot_hcm <- dt_gene[cohort == "ALL_HCM" & gene %in% hcm_genes]
dt_plot_dcm <- dt_gene[cohort == "ALL_DCM" & gene %in% dcm_genes]

x_max <- max(c(dt_plot_dcm$Upper_CI, dt_plot_hcm$Upper_CI)) + 0.1
x_min <- min(c(dt_plot_dcm$Lower_CI, dt_plot_hcm$Lower_CI)) / 2

list_hcm <- list()
for(i in seq_along(unique(dt_plot_hcm$gene))){
  dt_plot <- dt_plot_hcm[gene == unique(dt_plot_hcm$gene)[i]]
  custom_palette <- c("#1b9e77", "#d95f02", "black")
  p_out <- ggplot(dt_plot, aes(x = OR, y = score, color = prediction)) +
    geom_point(size = 3, position = position_dodge(width = 0.5)) +
    geom_errorbarh(aes(xmin = Lower_CI, xmax = Upper_CI), height = 0.2, position = position_dodge(width = 0.5)) +
    scale_x_continuous(trans = "log10", limits = c(x_min, x_max)) +  # Log-transform the x-axis
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +  # Add a dashed line where OR = 1
    scale_color_manual(values = custom_palette) +  # Apply the custom color palette
    # scale_x_continuous(limits = c(0, x_max)) +
    labs(
      title = paste0(unique(dt_plot_hcm$gene)[i], " (HCM)"),
      x = "log(Odds Ratio)",
      y = "",
      color = ""
    ) +
    theme_minimal() +
    theme(
      legend.position = "top",
      axis.text.y = element_text(size = 12),
      axis.text.x = element_text(size = 12, angle = 0),
      axis.line.x = element_line()
    )
  p_out
  list_hcm[[i]] <- p_out
}

list_dcm <- list()
for(i in seq_along(unique(dt_plot_dcm$gene))){
  dt_plot <- dt_plot_dcm[gene == unique(dt_plot_dcm$gene)[i]]
  custom_palette <- c("#1b9e77", "#d95f02", "black")
  p_out <- ggplot(dt_plot, aes(x = OR, y = score, color = prediction)) +
    geom_point(size = 3, position = position_dodge(width = 0.5)) +
    geom_errorbarh(aes(xmin = Lower_CI, xmax = Upper_CI), height = 0.2, position = position_dodge(width = 0.5)) +
    scale_x_continuous(trans = "log10", limits = c(x_min, x_max)) +  # Log-transform the x-axis
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +  # Add a dashed line where OR = 1
    scale_color_manual(values = custom_palette) +  # Apply the custom color palette
    # scale_x_continuous(limits = c(0, x_max)) +
    labs(
      title = paste0(unique(dt_plot_dcm$gene)[i], " (DCM)"),
      x = "log(Odds Ratio)",
      y = "",
      color = ""
    ) +
    theme_minimal() +
    theme(
      legend.position = "top",
      axis.text.y = element_text(size = 12),
      axis.text.x = element_text(size = 12, angle = 0),
      axis.line.x = element_line()
    )
  p_out
  list_dcm[[i]] <- p_out
}

# Arrange and annotate
p_gene_hcm <- ggarrange(plotlist = list_hcm, ncol = 2, nrow = ceiling(length(list_hcm)/2))
p_gene_hcm
p_gene_dcm <- ggarrange(plotlist = list_dcm, ncol = 2, nrow = ceiling(length(list_dcm)/2))
p_gene_dcm

p_gene_hcm <- annotate_figure(p_gene_hcm, top = text_grob("HCM Genes", size = 14, face = "bold"))
p_gene_dcm <- annotate_figure(p_gene_dcm, top = text_grob("DCM Genes", size = 14, face = "bold"))

ggsave("../results/Figure_case_enrichment_by_gene_HCM.png", p_gene_hcm, height = 3*ceiling(length(list_hcm)/2), width = 10, bg = "white")
ggsave("../results/Figure_case_enrichment_by_gene_DCM.png", p_gene_dcm, height = 3*ceiling(length(list_dcm)/2), width = 10, bg = "white")


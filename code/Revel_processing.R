rm(list = ls())
graphics.off()

library(data.table)

# Import REVEL data

dt_rv <- fread("../data/pathogenicity_scores/REVEL_cardiacG2P_GRCh38.csv.gz")


# Subset for desired columns

dt_rv <- dt_rv[,c("chromosome", "pos", "ref", "alt","REVEL")]

# Rename and reorder columns 

colnames(dt_rv) <- c("chromosome", "pos", "ref", "alt", "REVEL")

dt_rv <- dt_rv[, c("chromosome", "pos", "ref", "alt", "REVEL")]

dt_rv <- dt_rv[complete.cases(dt_rv)]

# Output Revel 

fwrite(dt_rv, "../data/processed/REVEL_cardiacG2P_GRCh38processed.csv.gz")



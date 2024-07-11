rm(list = ls())

library(data.table)

# Import AlphaMissense data

dt_am <- fread("./data/pathogenicity_scores/AlphaMissense_cardiacG2P_GRCh38.csv.gz")


dt_cv <- fread("../data/clinvar_data/ClinVar_HCM_DCM_variant_summary.csv.gz")

# Subset for desired columns

dt_am <- dt_am[,c("chromosome", "pos", "ref", "alt", "am_pathogenicity", "am_class", "external_gene_name", "ensembl_transcript_id")]

dt_am <- dt_am[,c("Chromosome", "Pos", "Ref", "Alt", "Gene", "Ensembl_transcript_id", "AlphaMissense_score")]

# Rename and reorder columns 

colnames(dt_am) <- c("Chromosome", "Pos", "Ref", "Alt", "AlphaMissense_score", "AlphaMissense_class", "Gene", "Ensembl_transcript_id")

dt_am <- dt_am[, c("Chromosome", "Pos", "Ref", "Alt", "Gene", "Ensembl_transcript_id", "AlphaMissense_score", "AlphaMissense_class")]


# Output
fwrite(dt_am, "./data/pathogenicity_scores/AlphaMissense_cardiacG2P_GRCh38processed.csv.gz")




               
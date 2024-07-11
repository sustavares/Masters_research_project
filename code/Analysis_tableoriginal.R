rm(list = ls())
graphics.off()

library(data.table)

### Import ==========
dt_cv <- fread("./data/clinvar_data/variant_summary_140324.txt.gz")
dt_cb <- fread("./data/pathogenicity_scores/CardioBoost_all_predictions.csv.gz")
dt_am <- fread("./data/pathogenicity_scores/AlphaMissense_cardiacG2P_GRCh38.csv.gz")
dt_rv <- fread("./data/pathogenicity_scores/REVEL_cardiacG2P_GRCh38.csv.gz")
dt_rv <- fread("./data/pathogenicity_scores/REVEL_cardiacG2P_GRCh38.csv.gz", na.strings = c("NA", "NaN", "", "."))
dt_mut <- fread("./data/mutations/cardiacG2P_hcm_dcm_canonical_CDS_point_mutations_GRCh38.csv.gz")

analysis_table <- fread("./data/analysis/analysis_table.csv")

analysis_table1 <- fread("./data/analysis/analysis_table1.csv")

### Format ========== 
colnames(dt_cv)
colnames(dt_cb)

dt_cv$Assembly <- NULL
colnames(dt_cv) <- c("Gene", "ClinicalSignificance_ClinVar", "ReviewStatus_ClinVar", "PhenotypeList_ClinVar", "OriginSimple_ClinVar", "Type_ClinVar", "PhenotypeCM_ClinVar", "Chromosome", "Pos", "Ref", "Alt")

colnames(dt_cb) <- c("Chromosome", "Pos", "Ref", "Alt", "Gene", "HGVSC_CardioBoost", "Pathogenicity_CardioBoost", "Classification_CardioBoost")

colnames(dt_am) <- c("Chromosome", "Pos", "Ref", "Alt", "Pathogenicity_AlphaMissense", "Classification_AlphaMissense")


dt_am$genome <- NULL 

dt_am$uniprot_id <- NULL

dt_am$transcript_id <- NULL

dt_am$protein_variant <- NULL

dt_am$ensembl_transcript_id <- NULL

dt_am$ensembl_transcript_version <- NULL


# Rename columns to merge 
colnames(dt_cb)[colnames(dt_cb) == "chrom"] <- "Chromosome"

colnames(dt_cv)[colnames(dt_cv) == "PositionVCF"] <- "pos"

colnames(dt_cv)[colnames(dt_cv) == "ReferenceAlleleVCF"] <- "ref"

colnames(dt_cv)[colnames(dt_cv) == "AlternateAlleleVCF"] <- "alt"

colnames(dt_cv)[colnames(dt_cv) == "GeneSymbol"] <- "Gene"

colnames(dt_cb)[colnames(dt_cb) == "gene"] <- "Gene"

colnames(dt_cv)[colnames(dt_cv) == "ClinicalSignificance"] <- "classification"

str(dt_cv)
str(dt_cb)

# Change str
dt_cv$Chromosome <- as.character(dt_cv$Chromosome)

# Reorder columns 
analysis_table2 <- analysis_table2[,c("Chromosome", "Pos", "Ref", "Alt", "Gene", "Ensembl_transcript_id", "Mutation_consequence", "REVEL_score", "AlphaMissense_score", "CardioBoost_score")]

# Remove rows with non-numeric values in pos column for REVEL
dt_rv <- dt_rv[complete.cases(dt_rv)]

duplicated(dt_rv, by = c("Chromosome", "Pos", "Ref", "Alt"), fromLast = TRUE)

# Remove duplicates
dt_rv <- unique(dt_rv)

analysis_table <- unique(analysis_table)

# Merge data tables

analysis_table <- merge(dt_rv, dt_mut, by = c("Chromosome", "Pos", "Ref", "Alt", "Gene", "Ensembl_transcript_id"), all = TRUE)

analysis_table1 <- merge(analysis_table, dt_am, by = c("Chromosome", "Pos", "Ref", "Alt", "Gene", "Ensembl_transcript_id"), all = TRUE)

test <- merge(analysis_table1, dt_cb, by=c("Pos_cds","Ref","Alt", "Gene", "Ensembl_transcript_id"), all.x=TRUE, allow.cartesian = TRUE)

analysis_table1 <- unique(analysis_table1)

### Output
fwrite(analysis_table, "./data/analysis/analysis_table.csv")

fwrite(analysis_table1, "./data/analysis/analysis_table1.csv")

fwrite(analysis_table2),



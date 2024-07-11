### Script that merges pathogenicity scores and ACs for downstream analysis
# Be aware some score (eg AlphaMissense and REVEL) are transcript specific and can have multiple scores per variant

rm(list = ls())
graphics.off()

library(data.table)

### Import ==========

# All canonical CDS point mutations
dt_mu <- fread("../data/cardiacG2P_hcm_dcm_canonical_CDS_point_mutations_GRCh38.csv.gz", na.strings = c("NA", "NaN", "", "."))

# CardioBoost predictions
dt_cb <- fread("../data/CardioBoost_all_predictions.csv.gz", na.strings = c("NA", "NaN", "", "."))

# REVEL predictions
dt_rv <- fread("../data/REVEL_cardiacG2P_GRCh38.csv.gz", na.strings = c("NA", "NaN", "", "."))

# AlphaMissense predictions
dt_am <- fread("../data/AlphaMissense_cardiacG2P_GRCh38.csv.gz", na.strings = c("NA", "NaN", "", "."))

# CADD predictions
dt_cad <- fread("../data/CADD_v1.7_hcm_dcm_cardiacG2P_all_transcript_CDS_splice.tsv.gz", na.strings = c("NA", "NaN", "", "."))

# ClinVar
dt_cv <- fread("../../ClinVar/data/ClinVar_variant_summary_cardiacG2P_hg38_QCed.txt.gz", na.strings = c("NA", "NaN", "", "."))

# AC AN
dt_a <- fread("../../cardiacG2P_variant_counts/results/cardiacG2P_variant_counts_simple.csv")

### Format ==========

# CLinVar
dt_cv <- dt_cv[,c("Chromosome","PositionVCF","ReferenceAlleleVCF","AlternateAlleleVCF","ClinicalSignificance","PhenotypeCM")]
colnames(dt_cv) <- c("chromosome","pos","ref","alt","ClinVar_significance","ClinVar_phenotype")
dt_cv$chromosome <- paste0("chr",dt_cv$chromosome)

# AlphaMissense
dt_am <- dt_am[,c("chromosome","pos","ref","alt", "am_pathogenicity","am_class")]
colnames(dt_am) <- c("chromosome","pos","ref","alt", "AlphaMissense_pathogenicity","AlphaMissense_class")

# CardioBoost
colnames(dt_cb) <- c("gene_name","ensembl_transcript_id","pos_cds","ref","alt","CardioBoost_pathogenicity","CardioBoost_isTrain")

# CADD
colnames(dt_cad) <- c("chromosome", "pos", "ref", "alt", "CADD_raw", "CADD_PHRED")
dt_cad$chromosome <- paste0("chr",dt_cad$chromosome)

# REVEL
dt_rv <- unique(dt_rv[,c("chromosome","pos","ref","alt", "REVEL")])
dt_rv$chromosome <- paste0("chr",dt_rv$chromosome)
dt_rv <- dt_rv[complete.cases(dt_rv)] # There are some rows with non-numeric values in the pos column for REVEL -- remove them
dt_rv <- dt_rv[!duplicated(dt_rv[, c("chromosome", "pos", "ref", "alt")])] # exclude duplicates 

# Cohort AC AN

# Split ID column by - delimeter into chromosome, pos, ref, alt
dt_a <- dt_a[,c("chromosome","pos","ref","alt") := tstrsplit(ID, "-", fixed = TRUE)]
dt_a$pos <- as.integer(dt_a$pos)
dt_a$ID <- NULL
dt_a$na_count <- NULL
dt_a$na_percent <- NULL

# Merge
dt <- unique(merge(dt_mu, dt_cv, by=c("chromosome", "pos","ref","alt"), all.x=TRUE))
dt <- unique(merge(dt, dt_am, by=c("chromosome", "pos","ref","alt"), all.x=TRUE))
dt <- unique(merge(dt, dt_cb, by=c("ensembl_transcript_id","gene_name","pos_cds","ref","alt"), all.x=TRUE, allow.cartesian = TRUE))
dt <- unique(merge(dt, dt_rv, by=c("chromosome","pos","ref","alt"), all.x=TRUE))
dt <- unique(merge(dt, dt_cad, by=c("chromosome","pos","ref","alt"), all.x=TRUE))
dt <- unique(merge(dt, dt_a, by=c("chromosome","pos","ref","alt"), all.x=TRUE))

# Check for duplicates
dt_duplicates <- dt[duplicated(dt, by = c("chromosome", "pos", "ref", "alt")) | duplicated(dt, by = c("chromosome", "pos", "ref", "alt"), fromLast = TRUE)]
dt_duplicates <- dt_duplicates[order(chromosome, pos, ref, alt)]

# Remove duplicates
dt_out <- dt[!dt_duplicates, on = c("chromosome", "pos", "ref", "alt")]
# (nrow(dt_out) + nrow(dt_duplicates)) == nrow(dt)

### Export ==========

fwrite(dt_out, "../data/analysis_table.csv.gz")



###

# x <- dt[!is.na(dt$ClinVar_significance)]
# x <- dt_out[gene_name == "TTN"]



rm(list = ls())
graphics.off()

library(data.table)

### Import ==========

# Import clinvar data
dt_cv <- fread("../data/clinvar_data/variant_summary_140324.txt.gz")

# List of gene names 
dcm_genes <- c("BAG3","DES","DSP","FLNC" ,"LMNA","MYH7","PLN","RBM20","SCN5A","TNNC1","TNNT2","TTN")
hcm_genes <- c("ACTC1","MYBPC3","MYH7","MYL2","MYL3","PLN","TNNI3","TNNT2","TPM1")

### Filtering ==========

# Subset gene names 
dt_cv <- dt_cv[GeneSymbol %in% c(hcm_genes, dcm_genes)]

# filters by review status
table(dt_cv$ReviewStatus)
dt_cv <- dt_cv[ReviewStatus %in% c("criteria provided, multiple submitters, no conflicts", "criteria provided, single submitter", "reviewed by expert panel")]

# filters by genome build 
dt_cv <- dt_cv[Assembly %in% c("GRCh38")]

# filters for germline variants
dt_cv <- dt_cv[OriginSimple %in% c("germline")]

# filters for single nucleotide variants
dt_cv <- dt_cv[Type %in% c("single nucleotide variant")]

# filters for clinical significance of variants
dt_cv <- dt_cv[ClinicalSignificance %in% c("Benign", "Pathogenic", "Pathogenic/Likely pathogenic", "Likely pathogenic", "Benign/Likely benign", "Likely benign")]

# Add coulumn for phenotype
dt_cv$PhenotypeCM <- NA
dt_cv$PhenotypeCM[dt_cv$PhenotypeList %like% "Hypertrophic cardiomyopathy" |
                    dt_cv$PhenotypeList %like% "hypertrophic cardiomyopathy"] <- "HCM"
dt_cv$PhenotypeCM[dt_cv$PhenotypeList %like% "Dilated cardiomyopathy" |
                    dt_cv$PhenotypeList %like% "dilated cardiomyopathy"] <- "DCM"
ind_dcm <- which(dt_cv$PhenotypeList %like% "Dilated cardiomyopathy" |
                   dt_cv$PhenotypeList %like% "dilated cardiomyopathy")
ind_hcm <- which(dt_cv$PhenotypeList %like% "Hypertrophic cardiomyopathy" |
                   dt_cv$PhenotypeList %like% "hypertrophic cardiomyopathy")
if(length(ind_hcm[ind_hcm %in% ind_dcm]) != length(ind_dcm[ind_dcm %in% ind_hcm])){ print("Error!")}
dt_cv$PhenotypeCM[ind_hcm[ind_hcm %in% ind_dcm]] <- "HCM/DCM"

# Subset desired columns
dt_cv <- dt_cv[,c("GeneSymbol", "ClinicalSignificance", "ReviewStatus", "Assembly", "PhenotypeList", "OriginSimple", "Type", "PhenotypeCM", "Chromosome", "PositionVCF", "ReferenceAlleleVCF", "AlternateAlleleVCF")]

dt_cv <- dt_cv[,c("GeneSymbol", "ClinicalSignificance", "PhenotypeCM", "Chromosome", "PositionVCF", "ReferenceAlleleVCF", "AlternateAlleleVCF")]

# Rename and Reorder columns 

colnames(dt_cv) <- c("gene_name", "ClinicalSignificance", "PhenotypeCM", "chromosome", "pos", "ref", "alt")

dt_cv <- dt_cv[, c("Chromosome", "Pos", "Ref", "Alt", "Gene", "ClinicalSignificance", "ReviewStatus", "Assembly", "PhenotypeList", "OriginSimple", "Type", "PhenotypeCM")]

### Output
fwrite(dt_cv, "./data/clinvar_data/variant_summary_140324.txt.gz", compress = "gzip")







# https://bioconductor.org/packages/release/bioc/vignettes/biomaRt/inst/doc/accessing_ensembl.R

rm(list = ls())
graphics.off()

# libraries
library(biomaRt)
library(data.table)


# variables ==========

dcm_genes <- c("BAG3","DES","DSP","FLNC" ,"LMNA","MYH7","PLN","RBM20","SCN5A","TNNC1","TNNT2","TTN")
hcm_genes <- c("ACTC1","MYBPC3","MYH7","MYL2","MYL3","PLN","TNNI3","TNNT2","TPM1")
gene_names <- c(hcm_genes, dcm_genes)

# Use the ENSEMBL Mart and select the human dataset for hg38
hsap_mart_38 <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

# dt_attributes <- listAttributes(mart = hsap_mart_38)
# table(dt_attributes$page)
# dt_attributes_features <- dt_attributes[dt_attributes$page == "feature_page",]
# dt_attributes_seq <- dt_attributes[dt_attributes$page == "sequences",]
# dt_filters <- listFilters(mart = hsap_mart_38)

# Fetch the canonical transcript IDs for genes
dt_canonical <- as.data.table(
  getBM(attributes = c("external_gene_name", 
                       "ensembl_transcript_id", 
                       "uniprotswissprot", 
                       "transcript_is_canonical"),
  filters = c("external_gene_name"),
  values = list(gene_names),
  mart = hsap_mart_38)
)

# Fetch the CDS coordinates
dt_exon <- as.data.table(
  getBM(attributes = c("external_gene_name", "ensembl_transcript_id", "ensembl_exon_id",
                       "chromosome_name", "strand", "rank",
                       "genomic_coding_start", "genomic_coding_end",
                       "cds_start", "cds_end"),
        filters = "ensembl_transcript_id",
        values = dt_canonical$ensembl_transcript_id,
        mart = hsap_mart_38)
)

# Fetch the CDS sequence
dt_seq <- as.data.table(
  getBM(attributes = c("external_gene_name", "ensembl_transcript_id",
                       "chromosome_name", "coding"),
        filters = "ensembl_transcript_id",
        values = dt_canonical$ensembl_transcript_id,
        mart = hsap_mart_38)
)

# Flag canonical transcripts
dt_exon <- merge(dt_exon, dt_canonical[,c("ensembl_transcript_id", "transcript_is_canonical")], by = "ensembl_transcript_id", all.x =  TRUE)
dt_seq <- merge(dt_seq, dt_canonical[,c("ensembl_transcript_id", "transcript_is_canonical")], by = "ensembl_transcript_id", all.x =  TRUE)

# Flag HCM and DCM genes 
dt_canonical$HCM_gene[dt_canonical$external_gene_name %in% hcm_genes] <- 1
dt_canonical$DCM_gene[dt_canonical$external_gene_name %in% dcm_genes] <- 1

dt_exon$HCM_gene[dt_exon$external_gene_name %in% hcm_genes] <- 1
dt_exon$DCM_gene[dt_exon$external_gene_name %in% dcm_genes] <- 1

dt_seq$HCM_gene[dt_seq$external_gene_name %in% hcm_genes] <- 1
dt_seq$DCM_gene[dt_seq$external_gene_name %in% dcm_genes] <- 1


### EXPORT ==========

# export canonical info
fwrite(dt_canonical, 
       "../data/biomart_cardiacG2P_transcript_ids_GRCh38.tsv", 
       sep = "\t")

# export exon info
fwrite(dt_exon, 
       "../data/biomart_cardiacG2P_exon_coordinates_GRCh38.tsv", 
       sep = "\t")

# export canonical IDs and chromosome
fwrite(dt_seq, 
       "../data/biomart_cardiacG2P_CDS_seq_GRCh38.tsv", 
       sep = "\t")




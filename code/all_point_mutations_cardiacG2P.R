rm(list = ls())
graphics.off()

library(data.table)
library(stringi)
library(stringr)

### Import ==========
dt_exon <- fread("../data/biomart_cardiacG2P_exon_coordinates_GRCh38.tsv")
dt_seq <- fread("../data/biomart_cardiacG2P_CDS_seq_GRCh38.tsv")
dt_aa <- fread("../data/AA_mutation_table.csv.gz")


### QC ==========

dt_removed <- data.frame()
dt_seq$cds_length <- nchar(dt_seq$coding)


# remove seq not divisible by 3
ncod <- dt_seq$cds_length/3
rm.id <- which(ncod%%1!=0)
if (length(rm.id) != 0){
  dt_removed <- rbind(dt_removed, dt_seq[rm.id,])
  dt_seq <- dt_seq[-rm.id,]
}

# remove seq with characters other than ATCG 
rm.id <- grep('[^ATGC]', dt_seq$coding)
if (length(rm.id) != 0){
  dt_removed <- rbind(dt_removed, dt_seq[rm.id,])
  dt_seq <- dt_seq[-rm.id,]
}

# remove seq that do not start with start codon
c1 <- substring(dt_seq$coding, 1, 3)
rm.id <- which(c1 != "ATG")
if (length(rm.id) != 0){
  dt_removed <- rbind(dt_removed, dt_seq[rm.id,])
  dt_seq <- dt_seq[-rm.id,]
}

# remove seq that do not end with stop codon
cn <- substring(dt_seq$coding, dt_seq$cds_length-2, dt_seq$cds_length)
rm.id <- which(cn != "TAG" & cn != "TAA" & cn != "TGA")
if (length(rm.id) != 0){
  dt_removed <- rbind(dt_removed, dt_seq[rm.id,])
  dt_seq <- dt_seq[-rm.id,]
}

rm(c1, cn, ncod, rm.id)


### Format ==========

# AA table
dt_aa <- dt_aa[,c("ref_codon", "alt_codon", "ref_aa", "alt_aa", "consequence")]

# Remove unneeded columns
dt_exon <- dt_exon[, c("ensembl_transcript_id", "genomic_coding_start", "genomic_coding_end", "strand", "rank", "cds_start", "cds_end")]

# Function to find the differing nucleotide in codons
find_mutation <- function(ref_codon, alt_codon) {
  for (i in 1:nchar(ref_codon)) {
    if (substr(ref_codon, i, i) != substr(alt_codon, i, i)) {
      return(list(ref = substr(ref_codon, i, i), alt = substr(alt_codon, i, i)))
    }
  }
  return(list(ref = NA, alt = NA)) # In case there's no difference
}

# Apply the function to each row and create new columns
dt_aa[, c("ref_coding", "alt_coding") := find_mutation(ref_codon, alt_codon), by = 1:nrow(dt_aa)]

# Function to find the position of change in the codon
find_change_position <- function(ref_codon, alt_codon) {
  for (i in 1:nchar(ref_codon)) {
    if (substr(ref_codon, i, i) != substr(alt_codon, i, i)) {
      return(i)
    }
  }
  return(NA)  # Return NA if no change is found
}

# Apply the function to each row and create a new column
dt_aa[, pos_codon := mapply(find_change_position, ref_codon, alt_codon)]

# Convert each gene to long format
list_out <- list()
for (j in seq_along(dt_seq$ensembl_transcript_id)){
  
transcript_id <- dt_seq$ensembl_transcript_id[j]
dt_exon_sub <- dt_exon[ensembl_transcript_id == transcript_id]
dt_seq_sub <- dt_seq[ensembl_transcript_id == transcript_id]

dt_exon_sub <- dt_exon_sub[complete.cases(dt_exon_sub),]
dt_exon_sub <- dt_exon_sub[order(rank)]
coding_seq <- unlist(str_split(dt_seq_sub$coding, ""))

complement <- function(coding_seq){
  ref <- coding_seq
  ref <- gsub("A", "A1", ref)
  ref <- gsub("C", "C1", ref)
  ref <- gsub("G", "C", ref)
  ref <- gsub("T", "A", ref)
  ref <- gsub("A1", "T", ref)
  ref <- gsub("C1", "G", ref)
  return(ref)
}

if(dt_exon_sub$strand[1] == -1){ 
  pos_seq <- unlist(lapply(1:nrow(dt_exon_sub), function(i) {
    seq(from = dt_exon_sub$genomic_coding_end[i], 
        to = dt_exon_sub$genomic_coding_start[i])
  }))
  ref <- complement(coding_seq)
}

if(dt_exon_sub$strand[1] == 1){ 
  pos_seq <- unlist(sequences_list <- lapply(1:nrow(dt_exon_sub), function(i) {
    seq(from = dt_exon_sub$genomic_coding_start[i], 
        to = dt_exon_sub$genomic_coding_end[i])
  }))
  ref <- coding_seq
}

three.split <- function(sequence){
  vec <- strsplit(sequence, "")[[1]]
  out <- paste0(vec[c(TRUE, FALSE, FALSE)], vec[c(FALSE, TRUE, FALSE)], vec[c(FALSE, FALSE, TRUE)])
  return(out)
}
ref_codon <- rep(three.split(paste0(coding_seq, collapse = "")), each = 3)

dt_long <- data.table(pos = pos_seq,
                      ref = ref, 
                      ref_coding = coding_seq,
                      pos_cds = 1:length(coding_seq),
                      ref_codon = ref_codon,
                      pos_codon = rep(1:3, length(coding_seq)/3))

dt_long2 <- unique(dt_long[dt_aa, on = c("ref_coding", "ref_codon", "pos_codon"), allow.cartesian = TRUE])
dt_long2 <- dt_long2[complete.cases(dt_long2),]
dt_long2 <- dt_long2[order(pos_cds),]

if(nrow(dt_long)*(3) != nrow(dt_long2)){
  print(paste0(transcript_id, " not merging correctly."))
  }

dt_long2$chromosome <- paste0("chr", dt_seq_sub$chromosome_name)
dt_long2$gene_name <- dt_seq_sub$external_gene_name
dt_long2$ensembl_transcript_id <- dt_seq_sub$ensembl_transcript_id
dt_long2$strand <- dt_exon_sub$strand[1]
if(dt_exon_sub$strand[1] == -1){dt_long2$alt <- complement(dt_long2$alt_coding)}
if(dt_exon_sub$strand[1] == 1){dt_long2$alt <- dt_long2$alt_coding}

list_out[[j]] <- dt_long2[,c("chromosome", "pos", "ref", "alt", "gene_name", "ensembl_transcript_id", "strand", "pos_cds", "ref_coding", "alt_coding", "ref_codon", "alt_codon", "ref_aa", "alt_aa", "consequence")]
print(paste0(transcript_id, " done! ", j))
}

dt_out <- rbindlist(list_out)

# adjust consequences
dt_out$consequence[dt_out$pos_cds %in% 1:3] <- "start lost"
dt_out$consequence[dt_out$ref_aa == "STOP" & dt_out$alt_aa != "STOP"] <- "stop lost"
table(dt_out$consequence)

# Subset canonical
dt_can <- dt_out[ensembl_transcript_id %in% dt_seq$ensembl_transcript_id[dt_seq$transcript_is_canonical == 1]]


### Export ==========

fwrite(dt_out, "../data/cardiacG2P_hcm_dcm_CDS_point_mutations_GRCh38.csv.gz", compress = "gzip")
fwrite(dt_can, "../data/cardiacG2P_hcm_dcm_canonical_CDS_point_mutations_GRCh38.csv.gz", compress = "gzip")

# ### Export VCFs ==========
# 
# # create vcf
# dt_vcf <- dt_out[,c("chromosome", "pos", "ref", "alt")]
# dt_vcf$id <- dt_out$gene_name
# dt_vcf$qual <- "."
# dt_vcf$filter <- "."
# dt_vcf$info <- "."
# dt_vcf$format <- "."
# dt_vcf <- dt_vcf[,c("chromosome", "pos", "id", "ref", "alt", "qual", "filter", "info", "format")]
# colnames(dt_vcf) <- c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT")
# dt_vcf <- dt_vcf[order(`#CHROM`, POS, ALT)]
# 
# dcm_genes <- c("BAG3","DES","DSP","FLNC" ,"LMNA","MYH7","PLN","RBM20","SCN5A","TNNC1","TNNT2","TTN")
# hcm_genes <- c("ACTC1","MYBPC3","MYH7","MYL2","MYL3","PLN","TNNI3","TNNT2","TPM1")
# fwrite(dt_vcf[ID %in% hcm_genes], 
#        "../data/hcm_cardiacG2P_canonical_CDS_point_mutations_GRCh38.vcf",
#        sep = "\t")
# fwrite(dt_vcf[ID %in% dcm_genes], 
#        "../data/dcm_cardiacG2P_canonical_CDS_point_mutations_GRCh38.vcf",
#        sep = "\t")

######

### Validate

# dt_am <- fread("../../../../../lms-ware-raw/live/resources/AlphaMissense/AlphaMissense_hg38_HCM_DCM_cardiacG2P_genes.tsv.gz")
# dt_am <- dt_am[,c("#CHROM", "POS", "REF", "ALT", "ensembl_transcript_id")]
# dt_am$`#CHROM` <- as.numeric(gsub("chr", "", dt_am$`#CHROM`))
# dt_am <- dt_am[ensembl_transcript_id %in% dt_out$ensembl_transcript_id]
# dt_test <- dt_vcf[dt_am, on = c("#CHROM", "POS", "REF", "ALT")]

# dt_ref21 <- fread("~/Dropbox/PhD/Data/Ensembl/Reference/Formatted/Human_REF_Ensembl_GRC38_v94_chr21.csv")
# colnames(dt_ref21) <- c("pos", "ref")
# dt_ref21$tmp <- "giggerdy"
# test <- dt_ref21[dt_out, on = c("pos", "ref")]



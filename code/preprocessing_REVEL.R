rm(list = ls())
graphics.off()

library(data.table)

# Import 
dt <- fread("/rds/general/project/lms-ware-raw/live/resources/REVEL/revel_with_transcript_ids.csv.gz")
dt_snp <- fread("../data/cardiacG2P_hcm_dcm_CDS_point_mutations_GRCh38.csv.gz")
# Format
dt <- dt[,c("chr", "grch38_pos", "ref", "alt", "REVEL", "Ensembl_transcriptid")]
colnames(dt) <- c("chromosome", "pos", "ref", "alt", "REVEL", "Ensembl_transcriptid")
dt_snp$chromosome <- gsub("chr", "", dt_snp$chromosome)
dt$pos <- as.integer(dt$pos)
# Merge
dt <- unique(merge(dt, dt_snp[,c("chromosome", "pos", "ref", "alt")]))
# Subset duplicates
dt_duplicates <- dt[duplicated(dt[, .(chromosome, pos, ref, alt)]) | duplicated(dt[, .(chromosome, pos, ref, alt)], fromLast = TRUE)]
# Split and expand the Ensembl_transcriptid column
dt[, Ensembl_transcriptid := strsplit(Ensembl_transcriptid, ";")]
dt_long <- dt[, .(chromosome, pos, ref, alt, REVEL, Ensembl_transcriptid = unlist(Ensembl_transcriptid)), by = .(chromosome, pos, ref, alt, REVEL)]
dt_long <- dt_long[,c("chromosome", "pos", "ref", "alt", "REVEL", "Ensembl_transcriptid")]
# Subset cardiacG2P transcripts
dt_long <- dt_long[Ensembl_transcriptid %in% dt_snp$ensembl_transcript_id]
colnames(dt_long)[colnames(dt_long) == "Ensembl_transcriptid"] <- "ensembl_transcript_id"
# Export
fwrite(dt_long, "../data/REVEL_cardiacG2P_GRCh38.csv.gz", compress = "gzip")


#####

# dt <- fread("/rds/general/project/lms-ware-raw/live/resources/REVEL/revel_with_transcript_ids.csv.gz")
# # dt_g2p <- fread("../data/biomart_cardiacG2P_transcript_ids_GRCh38.tsv")
# dt_snp <- fread("../data/cardiacG2P_hcm_dcm_CDS_point_mutations_GRCh38.csv.gz")
# # dt <- dt[,c("chr", "grch38_pos", "ref", "alt", "REVEL", "Ensembl_transcriptid")]
# dt <- dt[,c("chr", "grch38_pos", "ref", "alt", "REVEL")]
# colnames(dt) <- c("chromosome", "pos", "ref", "alt", "REVEL")
# dt <- merge(dt, dt_g2p[,c("ensembl_transcript_id","external_gene_name")], by = "ensembl_transcript_id", all.x = TRUE)
# 
# # dt <- dt[Ensembl_transcriptid %in% dt_g2p$ensembl_transcript_id]
# # colnames(dt) <- c("chromosome", "pos", "ref", "alt", "REVEL", "ensembl_transcript_id")
# dt <- merge(dt, dt_g2p[,c("ensembl_transcript_id","external_gene_name")], by = "ensembl_transcript_id", all.x = TRUE)
# fwrite(dt, "../data/REVEL_cardiacG2P_GRCh38.csv.gz", compress = "gzip")

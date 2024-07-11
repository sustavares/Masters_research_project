rm(list = ls())
graphics.off()

library(data.table)

dt <- fread("/rds/general/project/lms-ware-raw/live/resources/MVP/MVP_score_hg19.txt.gz")
dt_g2p <- fread("../data/biomart_cardiacG2P_transcript_ids_GRCh38.tsv")

# liftover canonical pos


# dt <- dt[,c("chr", "grch38_pos", "ref", "alt", "REVEL", "Ensembl_transcriptid")]
# dt <- dt[Ensembl_transcriptid %in% dt_g2p$ensembl_transcript_id]
# colnames(dt) <- c("chromosome", "pos", "ref", "alt", "REVEL", "ensembl_transcript_id")
# dt <- merge(dt, dt_g2p[,c("ensembl_transcript_id","external_gene_name")], by = "ensembl_transcript_id", all.x = TRUE)
# fwrite(dt, "../data/M-CAPL_cardiacG2P_GRCh38.csv.gz", compress = "gzip")

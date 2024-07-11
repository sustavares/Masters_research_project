# http://bejerano.stanford.edu/mcap/
# M-CAP only scores rare missense variants: hg19, ENSEMBL 75 missense, ExAC v0.3 where no super population has minor allele frequency above 1%. 
# If a missense variant has no M-CAP score, the M-CAP prediction should be assumed to be likely benign.

rm(list = ls())
graphics.off()

library(data.table)

dt <- fread("/Volumes/RDS/project/lms-ware-raw/live/resources/M-CAP/mcap_v1_4.txt.gz")
dt_g2p <- fread("../data/biomart_cardiacG2P_transcript_ids_GRCh38.tsv")

# liftover canonical pos


# dt <- dt[,c("chr", "grch38_pos", "ref", "alt", "REVEL", "Ensembl_transcriptid")]
# dt <- dt[Ensembl_transcriptid %in% dt_g2p$ensembl_transcript_id]
# colnames(dt) <- c("chromosome", "pos", "ref", "alt", "REVEL", "ensembl_transcript_id")
# dt <- merge(dt, dt_g2p[,c("ensembl_transcript_id","external_gene_name")], by = "ensembl_transcript_id", all.x = TRUE)
# fwrite(dt, "../data/M-CAPL_cardiacG2P_GRCh38.csv.gz", compress = "gzip")

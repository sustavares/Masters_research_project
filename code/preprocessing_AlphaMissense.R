# Check isoforms for missing canonical transcripts. 

rm(list = ls())
graphics.off()

library(data.table)

### Format by ensembl_transcript_id

### Import ==========

# AlphaMissense
# dt_am <- fread("../../../../../lms-ware-raw/live/resources/AlphaMissense/AlphaMissense_isoforms_hg38.tsv.gz")
dt_am <- fread("../../../../../lms-ware-raw/live/resources/AlphaMissense/AlphaMissense_hg38.tsv.gz")

# G2P transcript IDs
dt_g2p <- fread("../data/biomart_cardiacG2P_transcript_ids_GRCh38.tsv")

# # G2P all point mutations
# dt <- fread("../data/cardiacG2P_hcm_dcm_CDS_point_mutations_GRCh38.csv.gz")

### Format and merge ==========

# Subset and rename desired columns 
dt_am <- dt_am[, .(`#CHROM`, POS, REF, ALT, transcript_id, uniprot_id, am_pathogenicity, am_class)]
colnames(dt_am) <- c("chromosome", "pos", "ref", "alt", "transcript_id", "uniprot_id", "am_pathogenicity", "am_class")

colnames(dt_g2p)[colnames(dt_g2p) == "uniprotswissprot"] <- "uniprot_id"
dt_g2p_ensembl <- dt_g2p[,c("external_gene_name", "ensembl_transcript_id")]
dt_g2p_uniprot <- dt_g2p[,c("external_gene_name", "uniprot_id")]

# Split transcript_id into two columns by the . delimiter 
dt_am[, c("ensembl_transcript_id", "ensembl_transcript_version") := tstrsplit(transcript_id, ".", fixed = TRUE)]
dt_am$transcript_id <- NULL
dt_am$ensembl_transcript_version <- NULL

# Merege and reorder
dt_ensembl <- merge(dt_am, dt_g2p_ensembl, by = c("ensembl_transcript_id"))
dt_ensembl <- dt_ensembl[,c("chromosome","pos","ref","alt","am_pathogenicity","am_class","external_gene_name","ensembl_transcript_id")]

dt_uniprot <- merge(dt_am, dt_g2p_uniprot, by = c("uniprot_id"))
dt_uniprot <- dt_uniprot[,c("chromosome","pos","ref","alt","am_pathogenicity","am_class","external_gene_name","uniprot_id")]

# Merge
dt <- merge(dt_ensembl, dt_uniprot, 
            by = c("chromosome","pos","ref","alt","am_pathogenicity","am_class","external_gene_name"),
            all = TRUE)

# Label canonical transcripts
dt$transcript_is_canonical <- NA
dt$transcript_is_canonical[dt$ensembl_transcript_id %in% dt_g2p$ensembl_transcript_id[dt_g2p$transcript_is_canonical == 1]] <- 1
dt$transcript_is_canonical[dt$ensembl_transcript_id %in% dt_g2p$uniprot_id[dt_g2p$transcript_is_canonical == 1]] <- 1

### Export ==========

fwrite(dt, "../data/AlphaMissense_cardiacG2P_GRCh38.csv.gz", compress = "gzip")

#####

# ttn <- dt_am[external_gene_name == "TTN"]
# table(ttn$ensembl_transcript_id)

# abc <- unique(dt[,c("external_gene_name","ensembl_transcript_id","uniprot_id","transcript_is_canonical")])

# x <- as.data.table(table(dt$external_gene_name, dt$ensembl_transcript_id))
# x <- x[N != 0]
# x <- x[order(V1, N)]
# colnames(x) <- c("external_gene_name", "ensembl_transcript_id", "N")
# 
# dt_g2p <- fread("../data/biomart_cardiacG2P_transcript_ids_GRCh38.tsv")
# 
# test <- merge(dt_g2p, x, on = c("external_gene_name", "ensembl_transcript_id"), all.y = TRUE)
# 
# y <- fread("../data/AlphaMissense_all_isoforms_cardiacG2P_GRCh38.csv.gz")
# 
# x1 <- as.data.table(table(y$external_gene_name, y$ensembl_transcript_id))
# x1 <- x1[N != 0]
# x1 <- x1[order(V1, N)]
# colnames(x1) <- c("external_gene_name", "ensembl_transcript_id", "N")
# 
# table(x1$ensembl_transcript_id %in% x$ensembl_transcript_id)

# test <- merge(dt, dt_am, on = c("chromosome", "pos", "ref", "alt"), all.x = TRUE)
# 
# test <- test[gene_name %in% unique(dt$gene_name)]
# 
# # Assuming dt is your data.table
# duplicate_rows <- test[duplicated(.SD, by = c("chromosome", "pos", "ref", "alt")) | duplicated(.SD, by = c("chromosome", "pos", "ref", "alt"), fromLast = TRUE)]
# # Assuming dt is your data.table
# duplicate_rows <- test[duplicated(dt, by = c("chromosome", "pos", "ref", "alt")) | duplicated(test, by = c("chromosome", "pos", "ref", "alt"), fromLast = TRUE)]
# 
# 
# head(test)
# 
# test <- unique(test)
# 
# #####
# 
# 
# # Split transcript_id into two columns by the . delimiter 
# dt_am[, c("ensembl_transcript_id", "ensembl_transcript_version") := tstrsplit(transcript_id, ".", fixed = TRUE)]
# 
# colnames(dt_g2p)[colnames(dt_g2p) == "uniprotswissprot"] <- "uniprot_id"
# 
# dt_am <- dt_am[(ensembl_transcript_id %in% dt_g2p$ensembl_transcript_id) | (uniprot_id %in% dt_g2p$uniprot_id)]
# 
# # Merge by ensembl ID first
# test <- merge(dt_am, dt_g2p[,c("ensembl_transcript_id", "external_gene_name")], by = c("ensembl_transcript_id"), all.x = TRUE)
# 
# dt_na <- test[is.na(external_gene_name),]
# dt_na$external_gene_name <- NULL
# 
# # Merge by uniprot ID
# dt_na <- merge(dt_na, dt_g2p[,c("uniprot_id", "external_gene_name")], by = c("uniprot_id"), all.x = TRUE)
# 
# x <- unique(dt_g2p[,c("uniprot_id", "external_gene_name")])
# 
# 
# table(dt_am$external_gene_name, dt_am$ensembl_transcript_id)
# 
# table(dt_am$external_gene_name)

# dt_cb <- fread("../data/cardiacboost_cm_predict_hcm_dcm_cardiacG2P.csv")

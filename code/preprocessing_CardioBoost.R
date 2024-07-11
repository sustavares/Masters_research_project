rm(list = ls())
graphics.off()

library(data.table)

dt_cb <- fread("/Volumes/RDS/project/lms-ware-raw/live/resources/cardioboost_data/cm_prediction.txt")
dt_cb[,cds_pos := as.numeric(gsub("c\\.(\\d+).+", "\\1", cDot))] # extract cds_pos from cDot (ie c.1234A>G -> 1234)
dt_cb[, ensembl_transcript_id := gsub("^([A-Z]+[0-9]+)\\..*", "\\1", HGVSc)] # extract Ensembl transcript ID
dt_cb <- dt_cb[,c("gene", "ensembl_transcript_id", "cds_pos", "REF", "ALT", "pathogenicity", "isTrain")]
fwrite(dt_cb, "../data/CardioBoost_all_predictions.csv.gz")




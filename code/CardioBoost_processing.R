rm(list= ls ())

library(data.table)

# Import CardioBoost data

dt_cb <- fread("./data/pathogenicity_scores/CardioBoost_all_predictions.csv.gz")



# Subset desired columns
colnames(dt_cb) <- c("Gene","Ensembl_transcript_id","Pos_cds","Ref","Alt","CardioBoost_score","CardioBoost_isTrain")


dt_cb <- dt_cb[, c("Gene", "Ensembl_transcript_id", "Pos_cds", "Ref", "Alt", "CardioBoost_score")]



dcm_genes <- c("BAG3","DES","DSP","FLNC" ,"LMNA","MYH7","PLN","RBM20","SCN5A","TNNC1","TNNT2","TTN")
hcm_genes <- c("ACTC1","MYBPC3","MYH7","MYL2","MYL3","PLN","TNNI3","TNNT2","TPM1")

dt_cb <- dt_cb[classification %in% c("Benign", "Pathogenic", "Pathogenic/Likely Pathogenic", "Likely Pathogenic", "Benign/Likely Benign", "Likely Benign")]


# Rename and reorder columns 
colnames(dt_cb) <- c("Gene", "Pos", "Ref", "Alt", "CardioBoost_score")

dt_cb <- dt_cb[,c("Pos_cds", "Ref", "Alt", "Gene", "Ensembl_transcript_id", "CardioBoost_score")]

# Remove duplicates
dt_cb <- unique(dt_cb)

### Output
fwrite(dt_cb, "./data/pathogenicity_scores/CardioBoost_all_predictionsprocessed.csv.gz")




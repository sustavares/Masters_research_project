rm(list= ls ())

library(data.table)

# Import point mutations file

dt_mu <- fread("./data/mutations/cardiacG2P_hcm_dcm_canonical_CDS_point_mutations_GRCh38.csv.gz")


# Rename columns 

colnames(dt_mut) <- c("Chromosome", "Pos", "Ref", "Alt", "Gene", "Ensembl_transcript_id", "Strand", "Pos_cds", "Ref_coding", "Alt_coding", "Ref_codon", "Alt_codon", "Ref_aa", "Alt_aa", "Mutation_consequence")

dt_mu$Chromosome <- as.integer(sub("chr", "", dt_mut$Chromosome))

# Subset columns

dt_mu <- dt_mut[,c("Chromosome", "Pos", "Ref", "Alt", "Gene", "Ensembl_transcript_id", "Mutation_consequence")]

#  Output 

fwrite(dt_mu, "./data/mutations/cardiacG2P_hcm_dcm_canonical_CDS_point_mutations_GRCh38processed.csv.gz")



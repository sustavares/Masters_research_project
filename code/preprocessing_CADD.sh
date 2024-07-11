#!/bin/bash

# Script that uses tabix to subset the CADD tsv. Ensure bed file does not have "chr" suffix.

### load anaconda module and source environment
module load anaconda3/personal
source activate tabix

cd /rds/general/project/lms-ware-analysis/live/gjp/cardio_path_scores/code/

# Define file paths
FILE_RAW_BED="../data/cardiacG2P_hcm_dcm_all_transcript_CDS_splice.bed"
FILE_RAW_CADD="/rds/general/project/lms-ware-raw/live/resources/CADD/v1.7/whole_genome_SNVs.tsv.gz"
FILE_OUT="../data/CADD_v1.7_hcm_dcm_cardiacG2P_all_transcript_CDS_splice.tsv"

# Check format
gunzip -c ${FILE_RAW_CADD} | head

tabix -R ${FILE_RAW_BED} ${FILE_RAW_CADD} > ${FILE_OUT}
gzip ${FILE_OUT}




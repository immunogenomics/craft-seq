#!/bin/bash
# Date: 02/06/2025
# Author: Michelle Curtis, Yuriy Baglaenko
# This script performs alignment of R2 of DNA FASTQs, already demultiplexed into FASTQs for each cell, to an amplicon of interest using CRISPResso. 

OUT_DIR="../DNAOutput_R2/"
BATCH_BASE="../misc_data/CRISPRessoBatch_R2.txt" 
PLATE="CD45_1-1"
SEQ="TAACAGGTCCAGGAGAGCCTCAGATTATTTTTTGTAGAAGTGAAGCTGCACATCAAGGAGTAATTACCTGGAATCCCCCTCAAAGATCATTTCATAATTTTACCCTCTGTTATATAAAAGAGACAGGTAATTTGTGTAGAATTTAATTTCATCAGAAAAGAGAAATCAAGAATTTGAAAGTACTA" # Amplicon
FASTQ_DIR="../raw/CD45/Demultiplex/${PLATE}/"
OUTPUT="${OUT_DIR}${PLATE}"
mkdir -p $OUTPUT
cd $OUTPUT

PLATE_BATCH="${OUTPUT}/batch.batch"

head -1 $BATCH_BASE > $PLATE_BATCH

### Subset batch files to include only the files present in directory
ls $FASTQ_DIR*_R2.fastq.gz >> dna_filepaths.txt
sed '2,$s/.fastq/.fastq.gz/g' $BATCH_BASE | awk -v FASTQ_DIR="$FASTQ_DIR" -v OFS="\t" '$2=FASTQ_DIR$2' | awk 'NR!=1 {print}' | grep -f dna_filepaths.txt >> $PLATE_BATCH
rm dna_filepaths.txt

CRISPRessoBatch -o ${OUTPUT} --batch_settings ${PLATE_BATCH}  --skip_failed -p 4 -a ${SEQ} --bo ${OUTPUT}  

# Unzip all Alleles_frequency_tables for easy reading into R
find . -type f -name "Alleles_frequency_table.zip" -exec sh -c 'unzip -qq -d "${1%/*}" "$1"' _ {} \;
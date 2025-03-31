#!/bin/bash
# Date: 02/06/2025
# Author: Michelle Curtis, Yuriy Baglaenko
# This script demultiplexes paired DNA FASTQs into per-cell FASTQs for each plate, starting from plates with name format "merged_PLATE_DNA_R1_001.fastq.gz". This requires the cell barcode to be the first 10 nucleotides.

samples="../misc_data/CRISPR_DNABarcodes.txt"
mkdir -p Demultiplex

n="CD45_1-1"

# Remove trailing characters
n=$(echo "$n" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
r=`ls *$n*_R1_001.fastq.gz`
s=`ls *$n*_R2_001.fastq.gz`

echo "Demultiplexing $r and $s into Demultiplex/$n" 

paste <(zcat $r|paste - - - - ) <(zcat $s|paste - - - - ) | awk -v FS="\t" -v OFS="\n" 'FNR==NR {samples[$2]=$2; next} {barcode = substr($2,1,10); if(samples[barcode]) {print $1,$2,$3,$4>>""samples[barcode]"_R1.fastq"; print $5,$6,$7,$8>>""samples[barcode]"_R2.fastq"}}' $samples -

#Make a directory for the files per samples. Label with their code. And anything else you want. 
mkdir -p "Demultiplex/$n"

#Move all the files to the folder
mv -i *.fastq "./Demultiplex/$n"

#gip them all back up to save space
gzip ./Demultiplex/$n/*.fastq

#!/bin/bash
# Date: 02/06/2025
# Author: Michelle Curtis, Yuriy Baglaenko
# This script performs alignment of ADT FASTQs to an ADT reference using kallisto bustools. 

OUT_DIR="/data/srlab1/mcurtis/GSK/Experiments/20240924_FBXO11/ADTOutput/"
WHITELIST="/data/srlab2/yb966/GSKProject/Barcodes/RNAWhiteList.txt"
T2G="/data/srlab2/yb966/GSKProject/kite/featuremap/ADTFull/FeaturesMismatch.t2g"
INDEX="/data/srlab2/yb966/GSKProject/kite/featuremap/ADTFull/FeaturesMismatch.idx"

for ADT_barcode in $(cat ADTBarcodes.txt)
do

    r="/data/srlab1/GSK/Data/FBXO11/1727129250/fastq/merged_${ADT_barcode}_ADT_R1_001.fastq.gz"
    n=$ADT_barcode    
    s=`sed 's/_R1_001.fastq.gz/_R2_001.fastq.gz/g' <<< $r`
    
    echo ${OUT_DIR}${n}
    mkdir -p ${OUT_DIR}${n}
    
    echo "kb count -i $INDEX -x 0,0,14:0,14,26:1,0,0 -t 4 -g $T2G -w $WHITELIST -o ${OUT_DIR}${n} $r $s" >> ${OUT_DIR}${n}/Kallisto_cmd.txt
    kb count -i $INDEX -x 0,0,14:0,14,26:1,0,0 -t 4 -g $T2G -w $WHITELIST -o ${OUT_DIR}${n} $r $s

done



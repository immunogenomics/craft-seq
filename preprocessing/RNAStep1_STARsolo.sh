#!/bin/bash
# Date: 02/06/2025
# Author: Michelle Curtis, Yuriy Baglaenko
# This script performs alignment of RNA FASTQs to hg38 using STARSolo. 

INDEX="/data/srlab1/mcurtis/spatial_projects/data/cos-pitzalis/STAR_alignment/reference/mimic-GRCh38-2020-A-cellranger/"
WHITELIST="/data/srlab2/yb966/scRNA-DNA-ADT/BECD45-PE-BULK/AnalysisFiles/384RNABarcodes.txt"
OUT_DIR="/data/srlab1/mcurtis/GSK/Experiments/20240924_FBXO11/RNAOutput/"

mkdir -p ${OUT_DIR}/Log

function do_count()
{
    local n=$1
	local R1=$2
    local R2=$3
    
cat << EOF | bsub 

#!/bin/bash
#BSUB -o ${OUT_DIR}Log/${n}.o
#BSUB -e ${OUT_DIR}Log/${n}.e
#BSUB -J STAR_$n
#BSUB -q bigmem
#BSUB -M 32000
#BSUB -n 7

/PHShome/jn467/STAR-2.7.6a/source/STAR --genomeDir $INDEX --runThreadN 7 --readFilesIn ${R2} ${R1} --readFilesCommand zcat \
--outFileNamePrefix ${OUT_DIR}${n}/ --outSAMtype None --outSAMmode None \
--soloType CB_UMI_Simple --soloCBwhitelist $WHITELIST --soloCBstart 1 --soloCBlen 14 --soloUMIstart 15 --soloUMIlen 12 \
--soloCBmatchWLtype 1MM_multi_pseudocounts --soloFeatures GeneFull

EOF
}

cat lsf_params_STARRNA | while read n R1 R2
do
do_count $n $R1 $R2
done

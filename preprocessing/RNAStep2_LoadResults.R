#!/usr/bin/env Rscript
# Date: 02/06/2025
# Author: Michelle Curtis, Yuriy Baglaenko
# This script loads results from STARSolo and concatenates them into a single counts matrix for a given experiment across multiple plates.

source("../misc_funcs.R")
source("../utils.R")

library(stringr)
library(Matrix)

rna_dir = "./RNAOutput"

## List subdirectories with *Solo.out/Gene/filtered*
rna_paths = list.dirs(path = rna_dir, full.names = TRUE, recursive = TRUE)[grepl('filtered',
                                  list.dirs(path = rna_dir, full.names = TRUE, recursive = TRUE))]
rna_paths %>% head(2)

exprs = do.call(cbind, lapply(rna_paths, function(run){
    name = str_split(run, '/RNAOutput/') %>% map(2) %>% str_replace('/Solo.out/GeneFull/filtered', '')
    exprs = read10x_mtx(run, name)
})
)

meta = data.frame(cell = colnames(exprs),
                  plate = exprs %>% colnames %>% str_split('_') %>% map(2) %>% unlist)
rownames(meta) = colnames(exprs)

saveRDS(exprs, 'exprsRNA_PTPRC_all.rds')
write_csv(meta, 'metadata_PTPRC_all.csv')
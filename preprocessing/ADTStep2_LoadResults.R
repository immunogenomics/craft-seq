#!/usr/bin/env Rscript
# Date: 02/06/2025
# Author: Michelle Curtis, Yuriy Baglaenko
# This script loads results from kallisto bustools and concatenates them into a single ADT counts matrix for a given experiment across multiple plates.

source("../misc_funcs.R")
source("../utils.R")

library(stringr)
library(Matrix)

adt_dir = "./ADTOutput"

## List ADT subdirectories
adt_paths = list.dirs(path = adt_dir, full.names = TRUE, recursive = TRUE)[
    endsWith(list.dirs(path = adt_dir, full.names = TRUE, recursive = TRUE), 'counts_unfiltered')]
adt_paths %>% head(2)

counts_ADT = do.call(cbind, lapply(adt_paths, function(run){
    name = str_split(run, '/ADTOutput/') %>% map(2) %>% str_replace('/counts_unfiltered', '')
    counts = read10x_kallisto(run, name)
})
)

saveRDS(counts_ADT, 'countsADT_PTPRC_all.rds')

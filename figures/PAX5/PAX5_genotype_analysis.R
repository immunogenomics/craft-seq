####################################################################################################
## Script purpose: Genotyping pipeline from Yuriy's Jupyter notebook
## Author: Zepeng Mu
## Date: Tue Apr  1 10:16:05 2025
####################################################################################################
library(matrixStats)
library(data.table)
library(viridis)
library(RColorBrewer)
library(pheatmap)
library(Matrix)
library(tidyverse)
library(Biostrings)
library(ggpubr)
library(ggrepel)
"%&%" <- function(a, b) paste0(a, b)

source("../DNA_filtering_Functions.R")

# Region A ----
AllelesPAX5_A <- fread("../data/allCells_allele_A.txt.gz")

#Fix naming of columns
AllelesPAX5_A <- AllelesPAX5_A %>% rename(bc = "DNA_cell") %>% separate(smp, into = c("Plate", "sample"))
AllelesPAX5_A <- dplyr::inner_join(fread("../data/DNABarcodes.csv", header = T, col.names = c("Well", "DNA_cell")),
                                   AllelesPAX5_A, by = "DNA_cell") # Merge with DNA barcodes

#Create a plate_well ID for easy cell identification by combining the barcode DNA with the well_ID
AllelesPAX5_A <- AllelesPAX5_A %>% mutate(plate_well = paste0(sample, DNA_cell))
#Group by this ID to create a grouped table. 
AllelesPAX5_A <- AllelesPAX5_A %>% dplyr::group_by(plate_well)

#Add the total # of reads per cell in a a new columbn for filtering. Future iterations of this analysis will simply multiply % of total * the number. 
AllelesPAX5_A <- AllelesPAX5_A %>% summarize("TotalReads" = sum(`#Reads`)) %>% dplyr::inner_join(AllelesPAX5_A)

#Fix name
AllelesPAX5_A <- dplyr::rename(AllelesPAX5_A, Barcode_DNA = sample, Well_ID = Well)
#Take top 10 alleles for visualiztion
AllelesPAX5_A <- AllelesPAX5_A %>% dplyr::group_by(plate_well) %>% dplyr::slice_head(n = 10)

Filter_Alleles_PAX5A <- Filtering_Cells_Read(AllelesPAX5_A)
Filter_Alleles_PAX5A  %>% ggplot(aes(`%Reads`)) + 
  geom_histogram(aes(y = ..density..),bins = 100) + theme_gy() + geom_vline(xintercept = 25, color = "red") + 
  facet_grid(~Barcode_DNA)

Filter_Alleles_PAX5A <- Filtering_Alleles(ungroup(Filter_Alleles_PAX5A), 25)

Filter_Alleles_PAX5A <- mutate(Filter_Alleles_PAX5A, Reference = "CAGCGGTGCTTCTCCTATGTGACTGGTTCTAACTACCCTTTCCCTTTCCTTTTGTTTCTGATCTGTTTCAGGACATGGAGGAGTGAATCAGCTTGGGGGGGTTTTTGTGAATGGACGGCCACTCCCGGATGTAGTCCGCCAGAGGATA")

#Trim length for visualization and plotting
Filter_Alleles_PAX5A <- 
  mutate(Filter_Alleles_PAX5A, 
         "Aligned_Sequence" = str_trunc(Aligned_Sequence, 75, "left", ellipsis = "")) %>% 
  mutate("Reference"= str_trunc(Reference, 75, "left", ellipsis = ""))%>% 
  mutate("Aligned_Sequence" = str_trunc(Aligned_Sequence, 40, "right", ellipsis = "")) %>% 
  mutate("Reference"= str_trunc(Reference, 40, "right", ellipsis = ""))

Filter_Alleles_PAX5A <- Filter_Alleles_PAX5A %>% 
  add_count(Aligned_Sequence) %>% 
  filter(n >= 10)

Filter_Alleles_PAX5A %>% Plotting_Alleles

Filter_Alleles_PAX5A_gen <- Filter_Alleles_PAX5A %>% 
  Genotyping_Cells() %>% 
  add_count(AllelicGenotype) %>% 
  filter(!str_detect(AllelicGenotype, fixed("L")))

Filter_Alleles_PAX5A_gen %>% 
  dplyr::inner_join(Filter_Alleles_PAX5A, by = "plate_well") %>% 
  dplyr::select(AllelicGenotype, Barcode_DNA) %>% 
  table()

summReadsA <- Filter_Alleles_PAX5A %>% 
  dplyr::group_by(plate_well) %>% 
  summarise(A.totalReads = sum(`#Reads`)) %>% 
  filter(A.totalReads >= 20)

allCellsAFltSumm <- Filter_Alleles_PAX5A_gen %>% 
  filter(plate_well %in% summReadsA$plate_well) %>% 
  separate(genotype, c("geno1", "geno2"), sep = "_", remove = F, fill = "right") %>% 
  mutate(geno2 = case_when(is.na(geno2) ~ geno1, T ~ geno2)) %>% 
  rowwise() %>% 
  mutate(A.hap1 = str_flatten(str_sub(geno1, c(11, 32, 33, 35), c(11, 32, 33, 35))),
         A.hap2 = str_flatten(str_sub(geno2, c(11, 32, 33, 35), c(11, 32, 33, 35))))

sum(table(c(allCellsAFltSumm$A.hap1, allCellsAFltSumm$A.hap2)))
sort(table(c(allCellsAFltSumm$A.hap1, allCellsAFltSumm$A.hap2)))

commonHap <- table(c(allCellsAFltSumm$A.hap1, allCellsAFltSumm$A.hap2))

allCellsAFlt <- Filter_Alleles_PAX5A_gen %>% 
  filter(plate_well %in% summReadsA$plate_well) %>% 
  separate(genotype, c("geno1", "geno2"), sep = "_", remove = F, fill = "right") %>% 
  mutate(geno2 = case_when(is.na(geno2) ~ geno1, T ~ geno2),
         smp = str_extract(plate_well, "^[ABM][1234]"),
         well = str_sub(plate_well, start = 3, end = nchar(plate_well))) %>% 
  rowwise() %>% 
  mutate(A.hap1 = str_flatten(str_sub(geno1, c(11, 32, 33, 35), c(11, 32, 33, 35))),
         A.hap2 = str_flatten(str_sub(geno2, c(11, 32, 33, 35), c(11, 32, 33, 35))),
         A.T_84 = (str_sub(A.hap1, 1, 1) == "C") + (str_sub(A.hap2, 1, 1) == "C"),
         A.T_105 = (str_sub(A.hap1, 2, 2) == "C") + (str_sub(A.hap2, 2, 2) == "C"),
         A.T_106 = (str_sub(A.hap1, 3, 3) == "C") + (str_sub(A.hap2, 3, 3) == "C"),
         A.T_108 = (str_sub(A.hap1, 4, 4) == "C") + (str_sub(A.hap2, 4, 4) == "C")) %>% 
  filter(A.hap1 %in% names(commonHap) & A.hap2 %in% names(commonHap))

fwrite(allCellsAFlt, col.names = T, row.names = F, sep = "\t", quote = F,
       file = "../data/allCells_haplo_geno_A_12plates.txt.gz")

# Region B ----
AllelesPAX5_B <- fread("../data/allCells_allele.txt.gz")

#Fix naming of columns
AllelesPAX5_B <- AllelesPAX5_B %>% rename(bc = "DNA_cell") %>% separate(smp, into = c("Plate", "sample"))
AllelesPAX5_B <- dplyr::inner_join(fread("../data/DNABarcodes.csv", header = T, col.names = c("Well", "DNA_cell")),
                                   AllelesPAX5_B, by = "DNA_cell") # Merge with DNA barcodes


#Create a plate_well ID for easy cell identification by combining the barcode DNA with the well_ID
AllelesPAX5_B <- AllelesPAX5_B %>% mutate(plate_well = paste0(sample, DNA_cell))
#Group by this ID to create a grouped table. 
AllelesPAX5_B <- AllelesPAX5_B %>% dplyr::group_by(plate_well)

#Add the total # of reads per cell in a a new columbn for filtering. Future iterations of this analysis will simply multiply % of total * the number. 
AllelesPAX5_B <- AllelesPAX5_B %>% 
  summarize(TotalReads = sum(`#Reads`)) %>% 
  dplyr::inner_join(AllelesPAX5_B)

#Fix name
AllelesPAX5_B <- dplyr::rename(AllelesPAX5_B, Barcode_DNA = sample, Well_ID = Well)
#Take top 10 alleles for visualiztion
AllelesPAX5_B <- AllelesPAX5_B %>% dplyr::group_by(plate_well) %>% dplyr::slice_head(n = 10)

Filter_Alleles_PAX5B <- Filtering_Cells_Read(AllelesPAX5_B)

Filter_Alleles_PAX5B %>% ggplot(aes(`%Reads`)) + 
  geom_histogram(aes(y = ..density..), bins = 100) +
  theme_gy() +
  geom_vline(xintercept = 20, color = "red") + 
  facet_grid(~Barcode_DNA)

Filter_Alleles_PAX5B <- Filtering_Alleles(ungroup(Filter_Alleles_PAX5B), 20)

Filter_Alleles_PAX5B <- mutate(Filter_Alleles_PAX5B, Reference = "GTTTGCCTGGGAGATCAGGGACCGGCTGCTGGCAGAGCGGGTGTGTGACAATGACACCGTGCCTAGCGTCAGTTCCATCAACAGGTGAGGGGCTCGTGCCTGTGGGGGTTGGGGATTTGGAGGGATGGCAGGGCATCCTGGAGGCTCT")

#Trim length for visualization and plotting
Filter_Alleles_PAX5B <- 
  mutate(Filter_Alleles_PAX5B, 
         "Aligned_Sequence" = str_trunc(Aligned_Sequence, 84, "right", ellipsis = "")) %>% 
  mutate("Reference"= str_trunc(Reference, 84, "right", ellipsis = ""))%>% 
  mutate("Aligned_Sequence" = str_trunc(Aligned_Sequence, 30, "left", ellipsis = "")) %>% 
  mutate("Reference"= str_trunc(Reference, 30, "left", ellipsis = ""))

Filter_Alleles_PAX5B <- Filter_Alleles_PAX5B %>% 
  add_count(Aligned_Sequence) %>% 
  filter(n >= 10)

Filter_Alleles_PAX5B %>% Plotting_Alleles

Filter_Alleles_PAX5B_gen <- Filter_Alleles_PAX5B %>% 
  Genotyping_Cells() %>% 
  add_count(AllelicGenotype)

Filter_Alleles_PAX5B_gen %>% 
  dplyr::inner_join(Filter_Alleles_PAX5B, by = "plate_well") %>% 
  dplyr::select(AllelicGenotype, Barcode_DNA) %>% 
  table()

summReadsB <- Filter_Alleles_PAX5B %>% 
  dplyr::group_by(plate_well) %>% 
  summarise(B.totalReads = sum(`#Reads`)) %>% 
  filter(B.totalReads >= 20)

allCellsBFltSumm <- Filter_Alleles_PAX5B_gen %>% 
  filter(plate_well %in% summReadsB$plate_well) %>% 
  filter(str_count(genotype, fixed("_")) < 2) %>% 
  separate(genotype, c("geno1", "geno2"), sep = "_",remove = F, fill = "right") %>% 
  mutate(geno2 = case_when(is.na(geno2) ~ geno1, T ~ geno2)) %>% 
  rowwise() %>% 
  mutate(B.hap1 = str_flatten(str_sub(geno1, c(10, 19, 20, 24), c(10, 19, 20, 24))),
         B.hap2 = str_flatten(str_sub(geno2, c(10, 19, 20, 24), c(10, 19, 20, 24))))

sum(table(c(allCellsBFltSumm$B.hap1, allCellsBFltSumm$B.hap2)))
sort(table(c(allCellsBFltSumm$B.hap1, allCellsBFltSumm$B.hap2)))

commonHapB <- table(c(allCellsBFltSumm$B.hap1, allCellsBFltSumm$B.hap2))

allCellsBFlt <- Filter_Alleles_PAX5B_gen %>% 
  filter(str_count(genotype, fixed("_")) < 2) %>% 
  separate(genotype, c("geno1", "geno2"), sep = "_",remove = F, fill = "right") %>% 
  mutate(geno2 = case_when(is.na(geno2) ~ geno1, T ~ geno2),
         smp = str_extract(plate_well, "^[ABM][1234]"),
         well = str_sub(plate_well, start = 3, end = nchar(plate_well))) %>% 
  rowwise() %>% 
  mutate(B.hap1 = str_flatten(str_sub(geno1, c(10, 19, 20, 24), c(10, 19, 20, 24))),
         B.hap2 = str_flatten(str_sub(geno2, c(10, 19, 20, 24), c(10, 19, 20, 24))),
         B.T_64 = (str_sub(B.hap1, 1, 1) == "C") + (str_sub(B.hap2, 1, 1) == "C"),
         B.T_73 = (str_sub(B.hap1, 2, 2) == "C") + (str_sub(B.hap2, 2, 2) == "C"),
         B.T_74 = (str_sub(B.hap1, 3, 3) == "C") + (str_sub(B.hap2, 3, 3) == "C"),
         B.T_78 = (str_sub(B.hap1, 4, 4) == "C") + (str_sub(B.hap2, 4, 4) == "C")) %>% 
  filter(B.hap1 %in% names(commonHapB) & B.hap2 %in% names(commonHapB))

fwrite(allCellsBFlt, col.names = T, row.names = F, sep = "\t", quote = F,
       file = "../data/allCells_haplo_geno_B_12plates.txt.gz")

# Combine A and B ----
idxSort <- fread("../data/idxSort.txt")

allA <- fread("../data/allCells_haplo_geno_A_12plates.txt.gz")
allB <- fread("../data/allCells_haplo_geno_B_12plates.txt.gz")

cmbHaploGeno <- allA %>% 
  dplyr::select(-(genotype:n)) %>% 
  dplyr::full_join(allB %>% dplyr::select(-(plate_well:n)), by = c("well", "smp")) %>% 
  mutate(smp = "PAX5_"%&%smp) %>% 
  dplyr::left_join(idxSort, by = c("well", "smp"))

cmbHaploGeno <- cmbHaploGeno %>% 
  filter(str_detect(smp, "PAX5_A") & !is.na(A.T_84) | str_detect(smp, "PAX5_B") & !is.na(B.T_64) |
           str_detect(smp, "PAX5_M") & !is.na(A.T_84) & !is.na(B.T_64))

fwrite(cmbHaploGeno, col.names = T, row.names = F, sep = "\t", quote = F,
       file = "../data/allCells_haplo_geno_AB_12plates.txt.gz")

## Heatmap ====
allCellsABMGenoNew <- fread("../data/allCells_haplo_geno_AB_12plates.txt.gz")

genoMtrxAM <- allCellsABMGenoNew %>%
  filter(cellType != "NTC" & str_detect(smp, "PAX5_A|PAX5_M")) %>% 
  dplyr::select(A.T_84:A.T_108) %>%
  mutate(across(.cols = A.T_84:A.T_108, .fns = as.integer)) %>%
  as.matrix()

genoMtrxBM <- allCellsABMGenoNew %>%
  filter(cellType != "NTC" & str_detect(smp, "PAX5_B|PAX5_M")) %>% 
  dplyr::select(B.T_64:B.T_78) %>%
  mutate(across(.cols = B.T_64:B.T_78, .fns = as.integer)) %>%
  as.matrix()

genoMtrx <- allCellsABMGenoNew %>% 
  dplyr::select(A.T_84:A.T_108, B.T_64:B.T_78) %>% 
  mutate(across(A.T_84:B.T_78, as.integer)) %>% 
  as.matrix()

rownames(genoMtrx) <- allCellsABMGenoNew$smp%&%"-"%&%allCellsABMGenoNew$well
genoMtrx <- t(genoMtrx)

genoEditedMtrx <- allCellsABMGenoNew %>% 
  filter(cellType == "Edited") %>% 
  dplyr::select(A.T_84:A.T_108, B.T_64:B.T_78) %>% 
  mutate(across(A.T_84:B.T_78, as.factor)) %>% 
  as.matrix()

genoEditedMtrx <- t(genoEditedMtrx)

haTop <- HeatmapAnnotation(
  Condition = allCellsABMGenoNew$cellType,
  Editing = str_sub(allCellsABMGenoNew$smp, 6, 6),
  N = anno_barplot(apply(genoMtrx, 2, function(x) sum(x > 0, na.rm = T)),
                   gp = gpar(fill = "grey40"), bar_width = 1, axis_param = list(labels_rot = 0)),
  col = list(Condition = c("NTC" = "lightgreen", "Edited" = "purple"),
             Editing = c("A" = "yellow2", "B" = "green4", "M" = "blue2")),
  annotation_legend_param = list(direction = "horizontal", nrow = 1),
  simple_anno_size = unit(0.2, "in")
)

meanEdit <- t(apply(genoMtrx, 1, function(x) 100 * table(na.omit(x)) / length(na.omit(x))))

haRight <- rowAnnotation(
  "Edit %" = anno_barplot(meanEdit,bar_width = 1, gp = gpar(fill = c("wheat", "darkorange", "red3")),
                          axis_param = list(labels_rot = 0)),
  annotation_legend_param = list(direction = "horizontal", nrow = 1),
  simple_anno_size = unit(0.5, "in")
)

genoMtrx0 <- genoMtrx
genoMtrx0[is.na(genoMtrx0)] <- 0

tmpH <- draw(Heatmap(genoMtrx0, cluster_rows = F, cluster_columns = T, show_column_names = F,
                     clustering_method_columns = "ward.D2",
                     clustering_distance_columns = "euclidean"))

editHeatmap <- Heatmap(
  genoMtrx, na_col = "grey",
  col = c("0" = "wheat", "1" = "darkorange", "2" = "red3"),
  name = "Genotype",
  cluster_rows = F, cluster_columns = F,
  column_order = column_order(tmpH),
  row_names_side = "left", show_column_names = F,
  row_split = c(rep("A", 4), rep("B", 4)),
  column_split = str_sub(allCellsABMGenoNew$smp, 6, 6),
  top_annotation = haTop,
  right_annotation = haRight,
  column_title = NULL, row_title = NULL,
  raster_quality = 96, use_raster = F,
  heatmap_legend_param = list(direction = "horizontal", nrow = 1)
)

pdf("../figs/editingHeatmap.pdf", width = 8, height = 4)
draw(editHeatmap, heatmap_legend_side = "top", annotation_legend_side = "top")
dev.off()

## Only edited cells ----
myCol <- circlize::colorRamp2(breaks = c(0, 1), colors = c("white", "red2"))

genoMtrxA <- allCellsABMGenoNew %>%
  filter(cellType != "NTC" & str_detect(smp, fixed("PAX5_A"))) %>% 
  dplyr::select(A.T_84:A.T_108) %>%
  mutate(across(.cols = A.T_84:A.T_108, .fns = as.integer)) %>%
  as.matrix()

hA <- Heatmap(cor(genoMtrxA), col = myCol, cluster_rows = F, cluster_columns = F,
              show_column_names = F,
              cell_fun = function(j, i, x, y, width, height, fill) {
                grid.text(signif(cor(genoMtrxA)[i, j], 2), x, y, rot = -45, gp = gpar(fontsize = 8))
              },
              width = unit(1.2, "in"), height = unit(1.2, "in"))

pdf("../figs/corrA.pdf", width = 5, height = 5)
draw(hA)
dev.off()

genoMtrxB <- allCellsABMGenoNew %>%
  filter(cellType != "NTC" & str_detect(smp, fixed("PAX5_B"))) %>% 
  dplyr::select(B.T_64:B.T_78) %>%
  mutate(across(.cols = B.T_64:B.T_78, .fns = as.integer)) %>%
  as.matrix()

hB <- Heatmap(cor(genoMtrxB), col = myCol, cluster_rows = F, cluster_columns = F,
              show_column_names = F,
              cell_fun = function(j, i, x, y, width, height, fill) {
                grid.text(signif(cor(genoMtrxB)[i, j], 2), x, y, rot = -45, gp = gpar(fontsize = 8))
              },
              width = unit(1.2, "in"), height = unit(1.2, "in"))

pdf("../figs/corrB.pdf", width = 5, height = 5)
draw(hB)
dev.off()

genoMtrxM <- allCellsABMGenoNew %>%
  filter(cellType != "NTC" & str_detect(smp, fixed("PAX5_M"))) %>% 
  dplyr::select(A.T_84:A.T_108, B.T_64:B.T_78) %>%
  mutate(across(.cols = A.T_84:B.T_78, .fns = as.integer)) %>%
  as.matrix()

hM <- Heatmap(cor(genoMtrxM), col = myCol, cluster_rows = F, cluster_columns = F,
              show_column_names = F,
              cell_fun = function(j, i, x, y, width, height, fill) {
                grid.text(signif(cor(genoMtrxM)[i, j], 2), x, y, rot = -45, gp = gpar(fontsize = 8))
              },
              width = unit(2.4, "in"), height = unit(2.4, "in"))

pdf("../figs/corrM.pdf", width = 8, height = 8)
draw(hM)
dev.off()

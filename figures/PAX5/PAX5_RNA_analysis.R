####################################################################################################
## Script purpose: RNA analysis on combined PAX5 for all 12 plates
## Author: Zepeng Mu
## Date: Tue Feb 25 07:05:06 2025
####################################################################################################
library(tidyverse)
library(data.table)
library(Seurat)
library(SeuratObject)
library(matrixStats)
library(sparseMatrixStats)
library(SCpubr)
library(Phoenix)
library(tidyplots)
library(harmony)
library(ggpubr)
"%&%" <- function(a, b) paste0(a, b)

basicProcess <- function(object = NULL, nfeatures = 1000) {
  object <- NormalizeData(object, normalization.method = "LogNormalize", scale.factor = 10000)
  object <- FindVariableFeatures(object, selection.method = "vst", nfeatures = nfeatures)
  
  allGenes <- rownames(object)
  object <- ScaleData(object, features = allGenes)
  
  object <- RunPCA(object, features = VariableFeatures(object = object))
  
  DimPlot(object, reduction = "pca") + NoLegend()
  
  object <- FindNeighbors(object, dims = 1:15)
  
  return(object)
}

# Load data ----
smps <- "PAX5_"%&%c("A"%&%1:4, "B"%&%1:4, "M"%&%1:4)

pax5Cmb <- read_rds("../data/pax5Cmb_12plates.rds")

# RNA ----
VlnPlot(pax5Cmb, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 1)

pax5Cmb@meta.data %>% 
  mutate(nCount_RNA = log10(nCount_RNA)) %>% 
  tidyplot(x = nCount_RNA, y = nFeature_RNA) %>% 
  add_data_points(size = 0.2) %>% 
  add_title("RNA") %>% 
  adjust_font(fontsize = 12)

pax5Cmb <- basicProcess(pax5Cmb, nfeatures = 1500)
pax5Cmb <- FindClusters(pax5Cmb, resolution = 0.2, graph.name = "RNA_snn", cluster.name = "OldCluster")
pax5Cmb <- RunUMAP(pax5Cmb, dims = 1:15, n.neighbors = 20, spread = 0.6, min.dist = 0.9)

g1 <- do_DimPlot(pax5Cmb, group.by = "orig.ident", pt.size = 0.3, plot.title = "Plates") +
  labs(x = "RNA UMAP1", y = "RNA UMAP2") +
  theme_zm() +
  theme(aspect.ratio = 1, legend.position = "bottom")

g2 <- do_DimPlot(pax5Cmb, group.by = "seurat_clusters", pt.size = 0.3, plot.title = "Clusters",
                 colors.use = c("0" = "darkgreen", "1" = "darkred", "2" = "blue2")) +
  labs(x = "RNA UMAP1", y = "RNA UMAP2") +
  theme_zm() +
  theme(aspect.ratio = 1, legend.position = "bottom")

pax5Cmb$condition <- ifelse(pax5Cmb$cellType == "NTC", "NTC", str_sub(pax5Cmb$orig.ident, "6", "6"))

## Harmony ----
pax5Cmb[["RNA"]] <- split(pax5Cmb[["RNA"]], f = pax5Cmb$orig.ident)
pax5Cmb <- IntegrateLayers(pax5Cmb, method = HarmonyIntegration, orig.reduction = "pca",
                           new.reduction = "harmony")

pax5Cmb <- FindNeighbors(pax5Cmb, reduction = "harmony", dims = 1:15, graph.name = c("harmony_nn", "harmony_snn"))
pax5Cmb <- FindClusters(pax5Cmb, resolution = 0.2, cluster.name = "harmony_clusters", graph.name = "harmony_snn")
pax5Cmb <- RunUMAP(pax5Cmb, reduction = "harmony", dims = 1:15, n.neighbors = 15,
                   reduction.name = "har_umap", spread = 0.8, min.dist = 1)

do_DimPlot(pax5Cmb, pt.size = 0.3, plot.title = "Cluster", reduction = "har_umap",
           group.by = "harmony_clusters",
           colors.use = c("0" = "darkgreen", "1" = "darkred", "2" = "blue2")) +
  theme_zm(base_size = 12) +
  labs(x = "RNA UMAP 1", y = "RNA UMAP 2") +
  theme(aspect.ratio = 1,
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.direction = "horizontal",
        legend.position = "inside",
        legend.position.inside = c(0.7, 0.1))

ggsave(width = 3, height = 3, filename = "../figs/umap_har_clu.pdf")

pax5Cmb <- JoinLayers(pax5Cmb)

DefaultAssay(pax5Cmb) <- "RNA"

pax5Cmb$batch <- str_sub(pax5Cmb$orig.ident, 6, 6)
pax5Cmb$condition <- pax5Cmb$batch%&%"-"%&%pax5Cmb$cellType

## Cluster analysis ----
g1 <- do_DimPlot(pax5Cmb,
                 pt.size = 0.1, plot.title = "Condition", reduction = "har_umap",
                 group.by = "cellType", colors.use = c("NTC" = "grey", "Edited" = "orange2")) +
  theme_zm(base_size = 12) +
  labs(x = "RNA UMAP 1", y = "RNA UMAP 2") +
  theme(aspect.ratio = 1,
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.direction = "horizontal",
        legend.position = "inside",
        legend.position.inside = c(0.7, 0.1))

ggsave(plot = g1, width = 3, height = 3, filename = "../figs/har_umap_condition.pdf")

g2 <- do_DimPlot(pax5Cmb, pt.size = 0.3, plot.title = "Cluster", reduction = "har_umap",
                 group.by = "harmony_clusters",
                 colors.use = c("0" = "darkgreen", "1" = "darkred", "2" = "blue2")) +
  theme_zm(base_size = 12) +
  labs(x = "RNA UMAP 1", y = "RNA UMAP 2") +
  theme(aspect.ratio = 1,
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.direction = "horizontal",
        legend.position = "inside",
        legend.position.inside = c(0.7, 0.1))

ggsave(g2, width = 3, height = 3, filename = "../figs/umap_har_clu.pdf")

g4 <- pax5Cmb@meta.data %>%
  dplyr::inner_join(geno %>% dplyr::select(smp, well, A.edit.binary, B.edit.binary),
                    by = c("orig.ident" = "smp", "well")) %>% 
  mutate(edited = case_when(A.edit.binary == "FALSE" & B.edit.binary == "FALSE" ~ "Unedited",
                            A.edit.binary != "FALSE" & B.edit.binary == "FALSE" ~ "A",
                            A.edit.binary == "FALSE" & B.edit.binary != "FALSE" ~ "B",
                            A.edit.binary != "FALSE" & B.edit.binary != "FALSE" ~ "A and B")) %>% 
  dplyr::group_by(edited,  harmony_clusters) %>%
  tally() %>%
  ggplot(aes(x = harmony_clusters, y = n, fill = edited)) +
  theme_zm(base_size = 12) +
  geom_col(position = "fill", width = 0.7) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(limits = c("Unedited", "A", "B", "A and B"),
                    values = c("grey", "yellow2", "green4", "blue2"), name = "") +
  labs(x = "Clusters", y = "Proportion")

ggsave(plot = g4, width = 3, height = 2.5, filename = "../figs/barplot_har_clu_edited.pdf")

## Markers ----
Idents(pax5Cmb) <- pax5Cmb$harmony_clusters

diffGene <- FindMarkers(pax5Cmb,
                        assay = "RNA", logfc.threshold = 0,
                        min.pct = 0.25, ident.1 = "0", ident.2 = "1")

diffGeneSig <- diffGene %>%
  rownames_to_column("gene") %>% 
  filter(p_val_adj < 0.01 & abs(avg_log2FC) > 0.5)

diffGene1 <- FindMarkers(pax5Cmb,
                         assay = "RNA", logfc.threshold = 0,
                         min.pct = 0.25, ident.1 = "0", ident.2 = "2")

diffGeneSig1 <- diffGene1 %>%
  rownames_to_column("gene") %>% 
  filter(p_val_adj < 0.01 & abs(avg_log2FC) > 0.5)

## CD79B FCER2 IGHG3 VPREB3
tmpGene <- "CD79B"
tmpStat <- data.frame(group.by = c("1", "2"),
                      p = signif(c(diffGene[tmpGene, "p_val_adj"], diffGene1[tmpGene, "p_val_adj"]), 3),
                      y.position = max(pax5Cmb@assays[["RNA"]]["data"][tmpGene, ]) * c(0.8, 0.7))

do_ViolinPlot(pax5Cmb, group.by = "harmony_clusters",
              features = tmpGene, font.size = 10, pt.size = 0.1, line_width = 0, order = F,
              plot_boxplot = F, colors.use = c("0" = "darkgreen", "1" = "darkred", "2" = "blue2")) +
  geom_text(data = tmpStat, aes(x = c(1, 2), y = y.position, label = p), size = 8, size.unit = "pt") +
  theme_zm() +
  labs(x = "Cluster", y = str_glue("*{tmpGene}*")) +
  theme(legend.position = "none", axis.title.y = ggtext::element_markdown())

ggsave(str_glue("../figs/violin_{tmpGene}.pdf"), width = 2, height = 1.5)

####################################################################################################
## Script purpose: Mixed effect model with interaction for PAX5 multiplexing analysis
## Author: Zepeng Mu
## Date: Tue Apr  1 12:56:14 2025
####################################################################################################
library(tidyverse)
library(data.table)
library(Seurat)
library(pbmcapply)
library(Matrix)
library(matrixStats)
library(sparseMatrixStats)
library(lme4)
library(lmerTest)
library(MASS)
library(aod)
library(qvalue)
library(Phoenix)
library(ggrastr)
library(ggridges)
library(SCpubr)
library(clusterProfiler)
"%&%" <- function(a, b) paste0(a, b)

ggQQplot <- function(pList, colVec = c("black", "orange", "red", "purple"),
                     sizeVec = rep(0.5, length(colVec)),
                     legendLabel = NULL, title = "", openRange = F, sampling = 1) {
  
  if (is.null(legendLabel)) {
    legendLabel = colVec
  }
  
  if (length(pList) > length(colVec)) {
    stop("pList needs to be shorter than colVec!!!")
  }
  
  if (openRange) {
    pList <- lapply(pList, function(x) x[x > 0 | x < 1])
  }
  
  if (sampling < 1) {
    pList <- lapply(pList, function(x) sample(x, length(x) * sampling))
  }
  
  g <- ggplot() +
    theme_zm() +
    labs(title = title,
         x = "Empirical -log<sub>10</sub> (P)",
         y = "Observed -log<sub>10</sub> (P)") +
    theme(axis.title.x = ggtext::element_markdown(),
          axis.title.y = ggtext::element_markdown())
  
  for (i in 1:length(pList)) {
    xValue <- -log10(ppoints(length(pList[[i]])))
    yValue <- pList[[i]]
    yValue[yValue == 0] <- min(yValue[yValue > 0])
    yValue <- sort(-log10(yValue), decreasing = T)
    g <- g + ggrastr::geom_point_rast(mapping = aes(!!xValue, !!yValue, color = !!legendLabel[i]),
                                      size = sizeVec[i])
  }
  
  g <- g +
    geom_abline(slope = 1, intercept = 0, col = "red", lty = 2) +
    scale_color_manual(limits = legendLabel, values = colVec, name = "")
  
  return(g)
}

safe.glm.nb <- function(tmpDt, model, modelNull, gene, varName, test = c("LRT", "Wald", "Both"), testVar = NULL) {
  tryCatch({
    res <- glm.nb(formula = model, data = tmpDt)
    
    resSumm <- summary(res)$coefficients
    
    tmpDf <- data.frame(gene = gene)
    for (tmpVar in varName) {
      addDf <- data.frame(effect = resSumm[tmpVar, "Estimate"],
                          stdErr = resSumm[tmpVar, "Std. Error"],
                          z = resSumm[tmpVar, "z value"],
                          pval = resSumm[tmpVar, "Pr(>|z|)"])
      
      colnames(addDf) <- c("effect."%&%tmpVar, "stdErr."%&%tmpVar, "z."%&%tmpVar, "pval."%&%tmpVar)
      tmpDf <- cbind(tmpDf, addDf)
    }
    
    if (test == "LRT") {
      resNull <- glm.nb(formula = modelNull, data = tmpDt)
      cmp <- anova(res, resNull)
      tmpDf <- tmpDf %>% add_column(LRstat = cmp$`LR stat.`[2], pvalLRT = cmp$`Pr(Chi)`[2])
    } else if (test == "Wald") {
      cmp <- aod::wald.test(Sigma = vcov(res), b = coef(res), Terms = testVar)
      tmpDf <- tmpDf %>% add_column(chisq = cmp$result$chi2[1], pvalWald = cmp$result$chi2[3])
    } else if (test == "Both") {
      resNull <- glm.nb(formula = modelNull, data = tmpDt)
      cmp <- anova(res, resNull)
      cmp1 <- aod::wald.test(Sigma = vcov(res), b = coef(res), Terms = testVar)
      tmpDf <- tmpDf %>% add_column(LRstat = cmp$`LR stat.`[2], pvalLRT = cmp$`Pr(Chi)`[2],
                                    chisq = cmp1$result$chi2[1], pvalWald = cmp1$result$chi2[3])
    }
    
    tmpDf
  },
  error = function(e) {
    # Return a data frame with NA values if not converging
    tmpDfErr <- data.frame(gene = gene)
    for (tmpVar in varName) {
      addDf <- data.frame(effect = NA_real_, stdErr = NA_real_, z = NA_real_, pval = NA_real_)
      colnames(addDf) <- c("effect."%&%tmpVar, "stdErr."%&%tmpVar, "z."%&%tmpVar, "pval."%&%tmpVar)
      tmpDfErr <- cbind(tmpDfErr, addDf)
    }
    
    if (test == "LRT") {
      tmpDfErr <- tmpDfErr %>% add_column(LRstat = NA_real_, pvalLRT = NA_real_)
    } else if (test == "Wald") {
      tmpDfErr <- tmpDfErr %>% add_column(chisq = NA_real_, pvalWald = NA_real_)
    } else if (test == "Both") {
      resNull <- glm.nb(formula = modelNull, data = tmpDt)
      tmpDfErr <- tmpDfErr %>% add_column(LRstat = NA_real_, pvalLRT = NA_real_,
                                          chisq = NA_real_, pvalWald = NA_real_)
    }
    
    tmpDfErr
  })
}

na.replace <- function(x, value = 0) {
  x[is.na(x)] <- value
  return(x)
}

# Load data ----
genes <- fread("../data/genes.txt", header = F)$V1 # Coding genes
pax5Cmb <- read_rds("../data/pax5Cmb_processed_12plates.rds")
allCellsABMGeno <- fread("../data/amplicon/allCells_haplo_geno_AB_12plates.txt.gz")

allCellsABMGenoFlt <- allCellsABMGeno %>%
  mutate(across(c(A.T_84:A.T_108, B.T_64:B.T_78), na.replace)) %>% 
  rowwise() %>% 
  mutate(A.edit.binary = as.factor(sum(A.T_84, A.T_105, A.T_108) > 0),
         A.edit.cnt = sum(A.T_84, A.T_105, A.T_108),
         B.edit.binary = as.factor(sum(B.T_74,B.T_78) > 0),
         B.edit.cnt = sum(B.T_74,B.T_78),
         total.edit = sum(A.T_84, A.T_105, A.T_108, B.T_74, B.T_78))

# Running glmer.nb ----
metaDf <- pax5Cmb@meta.data
rnaCnt <- pax5Cmb[["RNA"]]$counts

## A.binary ----
metaDfFiltered <- metaDf %>%
  rownames_to_column("cellName") %>%
  dplyr::select(-(FSC:cellType)) %>%
  dplyr::inner_join(allCellsABMGenoFlt, by = c("orig.ident" = "smp", "well")) %>%
  filter(cellType != "NTC") %>%
  mutate(nCount_RNA = log(nCount_RNA), orig.ident = factor(orig.ident))

rnaCntNew <- rnaCnt[, metaDfFiltered$cellName]
rnaCntNewNorm <- t(t(rnaCntNew) / colSums2(rnaCntNew) * 1e4)

summCnt <- data.frame(gene = rownames(rnaCntNew),
                      propZero = rowSums(rnaCntNew == 0) / ncol(rnaCntNew),
                      meanCnt = rowMeans2(rnaCntNew),
                      meanNormCnt = rowMeans2(rnaCntNewNorm))

resLinear <- fread("../results/allCells_glmNB_wald_12plates.txt.gz")
resLinear <- resLinear %>% 
  mutate(qval.A.T_84 = qvalue(pval.A.T_84)$qvalues,
         qval.A.T_105 = qvalue(pval.A.T_105)$qvalues,
         qval.A.T_108 = qvalue(pval.A.T_108)$qvalues,
         qval.B.T_74 = qvalue(pval.B.T_74)$qvalues,
         qval.B.T_78 = qvalue(pval.B.T_78)$qvalues)

pvalCutoff <- 0.05 / (8 * nrow(resLinear))

resLinearSig <- resLinear %>% 
  dplyr::select(gene, starts_with(fixed("pval."))) %>% 
  pivot_longer(cols = -1, names_pattern = "([a-zE]+).(.*)$", names_to = c("term", "pos")) %>% 
  pivot_wider(id_cols = c(gene, pos), names_from = "term", values_from = "value") %>% 
  filter(pval <= pvalCutoff & pos == "B.T_78")

B29Score1 <- pbmclapply(resLinearSig$gene, function(x) {
  return(rnaCntNew[x, ] * resLinear$effect.B.T_78[resLinear$gene == x])
}, mc.cores = 4) %>% Reduce(rbind, .)

B29Score1 <- colMeans(B29Score1)
hist(B29Score1, breaks = 50)

metaDfFiltered <- metaDfFiltered %>% mutate(B29Score = B29Score1[metaDfFiltered$cellName])

metaDfFiltered %>% 
  mutate(rsdl = residuals(lm(B29Score ~ nCount_RNA + orig.ident))) %>% 
  ggplot(aes(x = factor(B.T_78), y = rsdl)) +
  theme_zm() +
  geom_boxplot() +
  labs(x = "B.29", y = "B29 gene score")

ggsave(filename = "../figs/B29_B29Score.pdf", width = 3, height = 3)

metaDfFiltered %>% 
  mutate(rsdl = residuals(lm(B29Score ~ nCount_RNA + orig.ident))) %>% 
  ggplot(aes(x = factor(A.edit.binary), y = rsdl)) +
  theme_zm() +
  geom_boxplot() +
  labs(x = "A edited", y = "B29 gene score")

ggsave(filename = "../figs/A.edit.binary_B29Score.pdf", width = 3, height = 3)

metaDfFiltered %>% 
  mutate(rsdl = residuals(lm(B29Score ~ nCount_RNA + orig.ident + B.T_78))) %>% 
  ggplot(aes(x = B.T_78, y = rsdl, color = A.edit.binary)) +
  theme_zm() +
  geom_point_rast(size = 0.3, shape = 21, alpha = 0.4,
                  position = position_jitter(width = 0.3, height = 0)) +
  geom_smooth(method = "lm", linewidth = 0.5) +
  scale_color_manual(limits = c("FALSE", "TRUE"), values = c("skyblue2", "red2"), name = "A edited") +
  labs(x = "B.29", y = "B.29 gene score\nconditioning on B.29") +
  theme(aspect.ratio = 1)

ggsave(filename = "../figs/A.edit.binary_B29Score_interaction.pdf", width = 4, height = 3)

metaDfFiltered %>% 
  mutate(rsdl = residuals(lm(B29Score ~ nCount_RNA + orig.ident + B.T_78))) %>% 
  ggplot(aes(x = rsdl, color = A.edit.binary)) +
  theme_zm() +
  geom_hline(yintercept = c(0, 1), lty = 2, lwd = 0.4, color = "grey") +
  stat_ecdf(lwd = 0.7) +
  scale_color_manual(limits = c("FALSE", "TRUE"), values = c("skyblue2", "red2"), name = "A edited") +
  labs(y = "CDF", x = "B.29 gene score conditioning on\nB.29 genotype") +
  coord_cartesian(xlim = c(-3, 4.5)) +
  theme(aspect.ratio = 1,
        legend.position = "inside",
        legend.position.inside = c(0.7, 0.3))

ggsave(filename = "../figs/cdf_A.edit.binary_B29Score_interaction.pdf", width = 3, height = 3)

metaDfFiltered %>% 
  mutate(rsdl = residuals(lm(B29Score ~ nCount_RNA + orig.ident))) %>% 
  ggplot(aes(x = factor(B.T_78), y = rsdl, color = A.edit.binary)) +
  theme_zm() +
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.7, outliers = F) +
  geom_point_rast(size = 0.1, position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.8),
                  raster.dpi = 720, alpha = 0.3) +
  scale_color_manual(limits = c("FALSE", "TRUE"), values = c("skyblue2", "red2"), name = "A edited") +
  labs(x = "B.29", y = "B29 gene score") +
  theme(legend.position = "inside", legend.position.inside = c(0.5, 0.9),
        legend.direction = "horizontal")

ggsave(filename = "../figs/boxplot_A.edit.binary_B29Score_interaction.pdf", width = 2.5, height = 3)

fitFull <- lm(B29Score ~ B.T_78 + A.edit.binary + B.T_78:A.edit.binary + nCount_RNA + orig.ident, data = metaDfFiltered)
fitNull <- lm(B29Score ~ B.T_78 + A.edit.binary + nCount_RNA + orig.ident, data = metaDfFiltered)

anova(fitNull, fitFull)

summary(fitFull)$coefficients

fit0 <- lm(B29Score ~ B.T_78 + A.edit.binary + nCount_RNA + orig.ident,
           data = metaDfFiltered %>% filter(B.T_78 == 0))

fit1 <- lm(B29Score ~ B.T_78 + A.edit.binary + nCount_RNA + orig.ident,
           data = metaDfFiltered %>% filter(B.T_78 == 1))

fit2 <- lm(B29Score ~ B.T_78 + A.edit.binary + nCount_RNA + orig.ident,
           data = metaDfFiltered %>% filter(B.T_78 == 2))

list(summary(fit0)$coefficients[2, ], summary(fit1)$coefficients[2, ], summary(fit2)$coefficients[2, ]) %>% 
  Reduce(rbind, .) %>% 
  as.data.frame() %>% 
  mutate(B.29 = factor(c(0, 1, 2)),
         CIlow = Estimate - `Std. Error` * 1.96, CIup = Estimate + `Std. Error` * 1.96,
         pText = str_glue("p = {signif(`Pr(>|t|)`, 3)}")) %>% 
  ggplot(aes(x = Estimate, y = B.29, color = -log10(`Pr(>|t|)`), fill = -log10(`Pr(>|t|)`))) +
  theme_zm() +
  geom_vline(xintercept = 0, lty = 2, color = "grey", lwd = 0.4) +
  geom_point(size = 2, show.legend = T) +
  geom_linerange(aes(xmin = CIlow, xmax = CIup), linewidth = 0.6, show.legend = F) +
  ggrepel::geom_text_repel(aes(label = pText), show.legend = F) +
  scale_color_gradient(low = "grey70", high = "red2", name = "-log10(P)") +
  scale_fill_gradient(low = "grey70", high = "red2") +
  guides(fill = "none") +
  labs(x = "Effect of A edit on B.29 gene score") +
  theme(legend.key.height = rel(0.6), legend.key.width = rel(0.8))

ggsave(filename = "../figs/forest_A.edit.binary_B29Score_interaction.pdf", width = 3, height = 2)

metaDfFiltered %>% 
  mutate(rsdl = residuals(lm(B29Score ~ nCount_RNA + orig.ident))) %>% 
  ggplot(aes(x = nCount_RNA, y = rsdl)) +
  theme_zm() +
  geom_point(size = 1.5, alpha = 0.3, shape = 21) +
  geom_smooth()

modelFull <- cnt ~ A.edit.binary + B.T_78 + A.edit.binary:B.T_78 + nCount_RNA + orig.ident
modelNull <- cnt ~ A.edit.binary + B.T_78 + nCount_RNA + orig.ident

resAll <- pbmclapply(summCnt %>% filter(propZero < 0.7 & meanCnt >= 2 & gene %in% genes) %>% pull(gene), function(x) {
  tmpDt <- metaDfFiltered %>% mutate(cnt = rnaCntNew[x, ])
  tmpFit <- safe.glm.nb(
    tmpDt, modelFull, modelNull, x,
    varName = c("A.edit.binaryTRUE", "B.T_78", "A.edit.binaryTRUE:B.T_78"),
    test = "LRT"
  )
}, mc.cores = 8) %>% Reduce(rbind, .)

resAll <- resAll %>% filter(!is.na(pval.A.edit.binaryTRUE))

fwrite(resAll, quote = F, sep = "\t", row.names = F, col.names = T,
       file = "../../results/PAX5_R/allCells_glmNB_score_interaction_12plates.txt.gz")

resAll <- fread("../../results/PAX5_R/allCells_glmNB_score_interaction_12plates.txt.gz")

resAll <- resAll %>% 
  mutate(onesideP = case_when(sign(effect.B.T_78) * sign(`effect.A.edit.binaryTRUE:B.T_78`) > 0 ~ pvalLRT / 2,
                              T ~ 1 - pvalLRT / 2))

pvalCutoff <- 0.05 / (8 * nrow(resLinear))

resLinearSig <- resLinear %>% 
  dplyr::select(gene, starts_with(fixed("pval."))) %>% 
  pivot_longer(cols = -1, names_pattern = "([a-zE]+).(.*)$", names_to = c("term", "pos")) %>% 
  pivot_wider(id_cols = c(gene, pos), names_from = "term", values_from = "value") %>% 
  filter(pval <= pvalCutoff & pos == "B.T_78")

resAllFiltered <- resAll %>% filter(gene %in% resLinearSig$gene)

table(sign(resAllFiltered$effect.B.T_78) * sign(resAllFiltered$`effect.A.edit.binaryTRUE:B.T_78`))

resAllFiltered$qvalLRT <- qvalue(resAllFiltered$onesideP)$qvalues
sum(resAllFiltered$qvalLRT <= 0.05)

fwrite(resAllFiltered, quote = F, row.names = F, col.names = T, sep = "\t",
       file = "../results/113_gene_interaction.txt")

resAllFiltered %>% 
  dplyr::select(gene, starts_with(fixed("pval."))) %>% 
  pivot_longer(cols = -1, names_pattern = "([a-zE]+).(.*)$", names_to = c("term", "pos")) %>% 
  pivot_wider(id_cols = c(gene, pos), names_from = "term", values_from = "value") %>% 
  filter(pval <= 0.05 / (2 * nrow(.))) %>% 
  dplyr::group_by(pos) %>% 
  tally() %>% 
  filter(pos != "A.edit.binaryTRUE:B.T_78") %>% 
  ggplot(aes(x = pos, y = n)) +
  theme_zm() +
  geom_col(width = 0.5) +
  geom_text(aes(label = n)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(labels = c("A edited", "B.29")) +
  labs(x = "", y = "Number of significant genes") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

ggsave(filename = "../figs/bar_interaction_mainSig.pdf", width = 3, height = 3)

resAllFiltered %>% 
  mutate(qval.B.T_78 = p.adjust(pval.B.T_78, method = "BH"),
         qval.A.edit.binaryTRUE = p.adjust(pval.A.edit.binaryTRUE, method = "BH"),
         `qval.A.edit.binaryTRUE:B.T_78` = p.adjust(`pval.A.edit.binaryTRUE:B.T_78`, method = "BH")) %>% 
  dplyr::select(gene, starts_with(fixed("qval."))) %>% 
  pivot_longer(cols = -1, names_pattern = "([a-zE]+).(.*)$", names_to = c("term", "pos")) %>% 
  pivot_wider(id_cols = c(gene, pos), names_from = "term", values_from = "value") %>% 
  filter(qval <= 0.1 / 3) %>% 
  dplyr::group_by(pos) %>% 
  tally() %>% 
  ggplot(aes(x = pos, y = n)) +
  theme_zm() +
  geom_col(width = 0.8) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(limits = c("A.edit.binaryTRUE", "B.T_78", "A.edit.binaryTRUE:B.T_78"),
                   labels = c("A edited", "B.29", "A edited * B.29")) +
  labs(x = "", y = "Number of significant genes") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

## Examples ----
for (tmpGene in resAllFiltered$gene[resAllFiltered$qvalLRT <= 0.05]) {
  tmpDt <- metaDfFiltered %>% mutate(cnt = rnaCntNew[tmpGene, ])
  tmpPval <- resAll$onesideP[resAll$gene == tmpGene]
  
  metaDfFiltered %>% 
    mutate(rsdl = residuals(glm.nb(formula = cnt ~ nCount_RNA + orig.ident, data = tmpDt))) %>% 
    ggplot(aes(x = factor(B.T_78), y = rsdl, color = A.edit.binary)) +
    theme_zm() +
    geom_boxplot(outliers = F, width = 0.7, size = 0.4, show.legend = F, position = position_dodge(width = 0.8), lwd = 0.4) +
    geom_point_rast(position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.3), raster.dpi = 720,
                    size = 0.2, scale = 0.5, show.legend = F) +
    scale_color_manual(limits = c("FALSE", "TRUE"), values = c("skyblue2", "red2"), name = "A edited") +
    scale_x_discrete(labels = c("A/A", "A/G", "G/G")) +
    labs(title = "", subtitle = str_glue("Interaction LRT p={signif(tmpPval, 3)}"),
         x = "B.29", y = str_glue("RNA - {tmpGene}"))
  
  ggsave(str_glue("../figs/score_interaction_{tmpGene}.pdf"), height = 3, width = 3)
}

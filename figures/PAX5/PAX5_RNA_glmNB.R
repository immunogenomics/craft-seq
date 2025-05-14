####################################################################################################
## Script purpose: NB model for PAX5 multiplexing analysis with Wald test
## Author: Zepeng Mu
## Date: Fri Mar 21 15:42:09 2025
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
library(patchwork)
library(ggrastr)
library(ggrepel)
library(Phoenix)
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
    # Return a data frame with NA_real_ values if model not converging
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
genes <- fread("../data/genes.txt", header = F)$V1 # Protein-coding genes
pax5Cmb <- read_rds("../data/pax5Cmb_processed_12plates.rds")
allCellsABMGeno <- fread("../data/amplicon/allCells_haplo_geno_AB_12plates.txt.gz")

allCellsABMGenoFlt <- allCellsABMGeno %>% 
  mutate(across(c(A.T_84:A.T_108, B.T_64:B.T_78), na.replace))

# Running glmer.nb ----
metaDf <- pax5Cmb@meta.data
rnaCnt <- pax5Cmb[["RNA"]]$counts

## All cells ----
metaDfFiltered <- metaDf %>% 
  rownames_to_column("cellName") %>% 
  dplyr::select(-(FSC:cellType)) %>% 
  dplyr::inner_join(allCellsABMGenoFlt, by = c("orig.ident" = "smp", "well")) %>% 
  filter(cellType != "NTC") %>% 
  mutate(nCount_RNA = log(nCount_RNA), orig.ident = factor(orig.ident))

rnaCntNew <- rnaCnt[, metaDfFiltered$cellName]
rnaCntNewNorm <- t(t(rnaCntNew) / metaDfFiltered$nCount_RNA * 1e4)

summCnt <- data.frame(gene = rownames(rnaCntNew),
                      propZero = rowSums(rnaCntNew == 0) / ncol(rnaCntNew),
                      meanCnt = rowMeans2(rnaCntNew),
                      meanNormCnt = rowMeans2(rnaCntNewNorm))

modelFull <- cnt ~ A.T_84 + A.T_105 + A.T_106 + A.T_108 + B.T_64 + B.T_73 + B.T_74 + B.T_78 + nCount_RNA + orig.ident
modelNull <- cnt ~ nCount_RNA + orig.ident

resAll <- pbmclapply(summCnt %>% filter(propZero < 0.7 & meanCnt >= 2 & gene %in% genes) %>% pull(gene), function(x) {
  tmpDt <- metaDfFiltered %>% mutate(cnt = rnaCntNew[x, ])
  tmpFit <- safe.glm.nb(
    tmpDt, modelFull, modelNull, x,
    varName = c("A.T_84", "A.T_105", "A.T_106", "A.T_108", "B.T_64", "B.T_73", "B.T_74", "B.T_78"),
    test = "Both", testVar = 2:9
  )
  return(tmpFit)
}, mc.cores = 6) %>% Reduce(rbind, .)

resAll <- resAll %>% filter(!is.na(pval.A.T_84))

fwrite(resAll, quote = F, sep = "\t", row.names = F, col.names = T,
       file = "../results/allCells_glmNB_wald_12plates.txt.gz")

resAll$qvalLRT <- qvalue(resAll$pvalLRT)$qvalues

pvalCutoff <- 0.05 / (8 * nrow(resAll))

resAll %>% 
  dplyr::select(gene, starts_with(fixed("pval."))) %>% 
  pivot_longer(cols = -1, names_pattern = "([a-zE]+).(.*)$", names_to = c("term", "pos")) %>% 
  pivot_wider(id_cols = c(gene, pos), names_from = "term", values_from = "value") %>% 
  filter(pval <= pvalCutoff) %>% 
  dplyr::group_by(pos) %>% 
  tally() %>% 
  ggplot(aes(x = pos, y = n)) +
  theme_zm() +
  geom_col(width = 0.8) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(limits = c("A.T_84", "A.T_105", "A.T_106", "A.T_108", "B.T_64", "B.T_73", "B.T_74", "B.T_78"),
                   labels = c("A.14", "A.35", "A.36", "A.38", "B.15", "B.24", "B.25", "B.29")) +
  labs(x = "", y = "Number of\nsignificant genes") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

ggsave(filename = "../figs/bar_number_of_hits.pdf", width = 3, height = 3)

ggQQplot(list(resAll$pval.A.T_84, resAll$pval.A.T_105, resAll$pval.A.T_106, resAll$pval.A.T_108,
              resAll$pval.B.T_64, resAll$pval.B.T_73, resAll$pval.B.T_74, resAll$pval.B.T_78),
         colVec = c("red2", "orange2", "grey", "yellow3", "grey", "grey", "green3", "steelblue"),
         legendLabel = c("A.14", "A.35", "A.36", "A.38", "B.15", "B.24", "B.25", "B.29"))

### Examples ----
resAll %>% 
  dplyr::select(gene, starts_with("effect"), starts_with("stdErr"), starts_with("pval.")) %>% 
  pivot_longer(cols = -1, names_pattern = "([a-zE]+).(.*)$", names_to = c("term", "pos")) %>% 
  pivot_wider(id_cols = c(gene, pos), names_from = "term", values_from = "value") %>% 
  mutate(ciLow = effect - 1.96 * stdErr, ciUp = effect + 1.96 * stdErr) %>% 
  filter(gene %in% c("JCHAIN", "VPREB3", "CD79B")) %>% 
  ggplot(aes(y = pos, x = effect, fill = -log10(pval))) +
  theme_zm(base_family = "") +
  geom_point(size = 2.5, shape = 23, stroke = NA) +
  geom_linerange(aes(color = -log10(pval), xmin = ciLow, xmax = ciUp), lwd = 0.6) +
  scale_y_discrete(limits = rev(c("A.T_84", "A.T_105", "A.T_106", "A.T_108", "B.T_64", "B.T_73", "B.T_74", "B.T_78")),
                   labels = rev(c("A.14", "A.35", "A.36", "A.38", "B.15", "B.24", "B.25", "B.29"))) +
  scale_fill_gradient(low = "grey80", high = "red2") +
  scale_color_gradient(low = "grey80", high = "red2") +
  labs(x = "beta", y = "") +
  geom_vline(xintercept = 0, lty = 2, lwd = 0.4) +
  facet_wrap(~gene, scale = "free") +
  theme(legend.position = "bottom")

ggsave(filename = "../figs/forest.pdf", height = 3, width = 6)

for (tmpGene in c("JCHAIN", "VPREB3", "CD79B")) {
  tmpDt <- metaDfFiltered %>% mutate(cnt = rnaCntNew[tmpGene, ])
  tmpPval <- anova(
    glm.nb(formula = modelFull, data = tmpDt),
    glm.nb(formula = cnt ~ nCount_RNA + orig.ident, data = tmpDt)
  )
  
  metaDfFiltered %>% 
    mutate(rsdl = residuals(glm.nb(formula = cnt ~ A.T_84 + A.T_105 + A.T_106 + A.T_108 + B.T_64 + B.T_73 + B.T_74 + nCount_RNA + orig.ident, data = tmpDt))) %>% 
    ggplot(aes(x = factor(B.T_78), y = rsdl, color = factor(B.T_78))) +
    theme_zm() +
    geom_boxplot(outliers = F, width = 0.6, size = 0.4, show.legend = F) +
    geom_point_rast(position = position_jitter(width = 0.2), size = 0.2, scale = 0.5, show.legend = F) +
    scale_color_manual(limits = c("0", "1", "2"), values = c("wheat2", "darkorange", "red3")) +
    scale_x_discrete(labels = c("A/A", "A/G", "G/G")) +
    labs(title = "", subtitle = str_glue("LRT p={signif(tmpPval$`Pr(Chi)`[2], 3)}"),
         x = "B.29", y = str_glue("RNA - {tmpGene}"))
  
  ggsave(filename = str_glue("../figs/boxplot_{tmpGene}_B78_12plates.pdf"), height = 2.2, width = 1.6)
}

for (tmpGene in c("VPREB3")) {
  tmpDt <- metaDfFiltered %>% mutate(cnt = rnaCntNew[tmpGene, ])
  tmpPval <- anova(
    glm.nb(formula = modelFull, data = tmpDt),
    glm.nb(formula = cnt ~ A.T_105 + A.T_106 + A.T_108 + B.T_64 + B.T_73 + B.T_74 + B.T_78 + nCount_RNA + orig.ident, data = tmpDt)
  )
  
  metaDfFiltered %>% 
    mutate(rsdl = residuals(glm.nb(formula = cnt ~ A.T_84 + A.T_105 + A.T_106 + A.T_108 + B.T_64 + B.T_73 + B.T_74 + B.T_78 + nCount_RNA + orig.ident, data = tmpDt))) %>% 
    ggplot(aes(x = factor(A.T_84), y = rsdl, color = factor(A.T_84))) +
    theme_zm() +
    geom_boxplot(outliers = F, width = 0.6, size = 0.4, show.legend = F) +
    geom_point_rast(position = position_jitter(width = 0.2), size = 0.2, scale = 0.5, show.legend = F) +
    scale_color_manual(limits = c("0", "1", "2"), values = c("wheat2", "darkorange", "red3")) +
    scale_x_discrete(labels = c("A/A", "A/G", "G/G")) +
    labs(title = "", subtitle = str_glue("LRT p={signif(tmpPval$`Pr(Chi)`[2], 3)}"),
         x = "A.14", y = str_glue("RNA - {tmpGene}"))
  
  ggsave(filename = str_glue("../figs/boxplot_{tmpGene}_A84_12plates.pdf"), height = 2.2, width = 1.6)
}

# Box plots for all variants
mut <- c("A.T_84", "A.T_105", "A.T_106", "A.T_108", "B.T_64", "B.T_73", "B.T_74", "B.T_78")
mutName <- c("A.14", "A.35", "A.36", "A.38", "B.15", "B.24", "B.25", "B.29")
names(mutName) <- mut

for (tmpGene in c("JCHAIN", "VPREB3", "CD79B")) {
  tmpDt <- metaDfFiltered %>% mutate(cnt = rnaCntNew[tmpGene, ])
  tmpgList <- lapply(mut, function(tmpMut) {
    otherMut <- setdiff(mut, tmpMut)
    tmpModelNull <- as.formula(str_glue("cnt ~ {str_flatten(otherMut,collapse = '+')} + nCount_RNA + orig.ident"))
    
    tmpPval <- anova(glm.nb(formula = modelFull, data = tmpDt),
                     glm.nb(formula = tmpModelNull, data = tmpDt))
    
    tmpg <- tmpDt %>% 
      mutate(rsdl = residuals(glm.nb(formula = tmpModelNull, data = tmpDt)),
             myX := factor(!!sym(tmpMut))) %>% 
      ggplot(aes(x = myX, y = rsdl, color = myX)) +
      theme_zm(base_size = 8) +
      geom_boxplot(outliers = F, width = 0.6, size = 0.4, show.legend = F) +
      geom_point_rast(position = position_jitter(width = 0.2), size = 0.2, scale = 0.3, show.legend = F) +
      scale_color_manual(limits = c("0", "1", "2"), values = c("wheat2", "darkorange", "red3")) +
      scale_x_discrete(labels = c("A/A", "A/G", "G/G")) +
      labs(title = "", subtitle = str_glue("LRT p={signif(tmpPval$`Pr(Chi)`[2], 3)}"),
           x = mutName[tmpMut], y = str_glue("RNA - {tmpGene}")) +
      theme(plot.margin = unit(rep(0, 4), "in"))
    
    return(tmpg)
  })
  
  gout <- wrap_plots(tmpgList, nrow = 2) + plot_layout(axis_titles = "collect")
  ggsave(str_glue("../figs/boxplot_{tmpGene}_allVar_12plates.pdf"), plot = gout, height = 3, width = 4)
}

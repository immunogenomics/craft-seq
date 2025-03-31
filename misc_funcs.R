# Load helpful/necessary libraries
library(data.table)
library(magrittr)
library(patchwork)
#library(useful)
#library(pheatmap)
#library(RColorBrewer)
#library(singlecellmethods)
#library(Seurat)
#library(RANN)
#library(RcppAnnoy)
library(tidyverse)
#library(seaborn)

# set figure size
fig.size <- function(height = 4, width = 4){
    options(repr.plot.height = height, repr.plot.width = width)
}

table <- function(..., useNA = 'ifany'){
    base::table(..., useNA = useNA)
}

#Define theme_pres() function
theme_pres <- function(base_size = 23){ 
    theme_bw(base_size = base_size) + 
    theme(
        # If facetting, remove grey boxes
        strip.background = element_blank(),

        # Make axis components dark grey
        panel.border = element_rect(color = 'grey46'), 
        axis.text = element_text(color = 'grey46'),
        
        #text elements
        plot.title = element_text(hjust = 0.5, # center title
                                  size = (base_size + 2)), # increase title size
        legend.title = element_text(face = 'italic') # italics legend label
    )
}

theme_clean <- function(base_size = 23){ 
    theme_bw(base_size = base_size) + 
    theme(
        # If facetting, remove grey boxes
        strip.background = element_blank(),

        # Make axis components dark grey
        panel.border = element_rect(color = 'grey46'), 
        axis.text = element_text(color = 'grey46'),
        
        #text elements
        plot.title = element_text(hjust = 0.5, # center title
                                  size = (base_size + 2)), # increase title size
        legend.title = element_text(face = 'italic'), # italics legend label
        
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()
    )
}


#Define theme_g() function, a theme_grey with better fonts
theme_g <- function(base_size = 23){ 
    theme_grey(base_size = base_size) + 
    theme(
#         # If facetting, remove grey boxes
#         strip.background = element_blank(),

#         # Make axis components dark grey
#         panel.border = element_rect(color = 'grey46'), 
#         axis.text = element_text(color = 'grey46'),
        
        #text elements
        plot.title = element_text(hjust = 0.5, # center title
                                  size = (base_size + 2)), # increase title size
        legend.title = element_text(face = 'italic') # italics legend label
    )
}


bsub_generic <- function(script_use, save_use, email = FALSE, submit = FALSE, queue = "short", nice = NULL){
    str <- paste('-R select[hname!=cn001]', '-R select[hname!=cn002]', '-R select[hname!=cn003]' ,
                 '-R select[hname!=cn004]', '-R select[hname!=cn005]',
                 'Rscript', script_use, '-o', save_use)
    if(!is.null(nice)){
        str <- paste('-n', nice, str)
    }
    str <- paste('bsub -q', queue, str)
    if(email == FALSE){
        str <- paste('LSB_JOB_REPORT_MAIL=N', str) # Add flag to avoid spamming yourself
    }
    if(submit == FALSE){ # just print out the command
        str
    } else { # actually submit via lsf
        system(str)
    }
}

# return sparsity of a matrix
get_sparsity <- function(mat){
    if(class(mat) == 'dgCMatrix'){
        return(1 - (length(mat@x) / (mat@Dim[1] * mat@Dim[2])))
    } else if(class(mat) == 'matrix'){
        nnzero <- sum(mat == 0)
        return(nnzero / (dim(mat)[1] * dim(mat)[2]))
    }
}

# create color palette
#testPal <- colorRampPalette(brewer.pal(9, 'Set1'))

# general function for plotting points; default is umap. If other dims (such as PCA), enter them as a string in plot_other_dims
# e.g. c("PC_1", "PC_2")
plot_umap_general <- function(umap_dims, labels, pt_size = 0.1, pch = 19, clusters_keep = NULL,
                              do_points = TRUE, do_density = FALSE, title = NULL, legend_size = 3, 
                              legend_title = NULL, scale_man = FALSE, do_facet = FALSE,
                              do_labels = FALSE, text_label_size = 3.5, text_label_color = 'black',
                              override_legend = TRUE, plot_other_dims = NULL){
    if(!is.null(plot_other_dims)){
        axis_labels <- plot_other_dims
    } else{
        axis_labels <- c("UMAP1", "UMAP2")
    }
    plt_df <- cbind.data.frame(umap_dims[, axis_labels], labels = labels)
    plt <- plt_df %>% sample_frac %>% 
        ggplot(aes(x = get(axis_labels[1]), y = get(axis_labels[2]), col = labels, label = labels))
    if(do_points){
        plt <- plt + geom_point(pch = pch, size = pt_size)
    }
    if(do_density){
        plt <- plt + geom_density_2d()
    }
    if(scale_man){
        plt <- plt + scale_color_manual(values = testPal(labels %>% unique %>% length))
    }
    plt <- plt +
        labs(x = axis_labels[1], y = axis_labels[2], title = title, col = legend_title) +
        theme(plot.title = element_text(hjust = 0.5)) +
        theme_classic()
    if(override_legend){
        plt <- plt + guides(color = guide_legend(override.aes = list(stroke = 1, alpha = 1, 
                                                                     shape = 19, size = legend_size))) 
    }
    if(do_labels){
        label_coords <- plt_df %>%
            group_by(labels) %>%
            summarize(UMAP1 = mean(UMAP1), UMAP2 = mean(UMAP2))
        plt <- plt + geom_text(data = label_coords, col = text_label_color, size = text_label_size)
    }
    if(do_facet){
        plt <- plt + facet_wrap(~labels)
    }
    return(plt)
}

plot_violinFeature <- function(meta, group, data, feature, title = NULL){
    a <- cbind.data.frame(Cluster = meta[,group], Feature = data[feature,])
    plot <- a %>% ggplot(aes(x = Cluster, y = Feature, col = Cluster)) +
        geom_violin() +
        stat_summary(fun.y = mean, geom = 'point', col = 'black') +
        theme_classic() +
        labs(x = 'Cluster', y = feature, title = title)
    return(plot)
}

# threshold values in a vector (values) by a designated quantile
quanThresh <- function(values, quant){
    max <- quantile(values, probs = quant)
    new.values <- values
    new.values[new.values > max] <- max
    return(new.values)
}

# quick and dirty way to plot gene expression on a umap. Default input is gene expression data is gene X cell
plot_gene_hex <- function (data, umap_frame, gene, pt_size = 0.1, pch = 19, thresh = NULL, 
    do_points = TRUE, do_density = FALSE, title = NULL, transposed = FALSE, plot_other_dims = NULL) 
{
    if (!is.null(plot_other_dims)) {
        axis_labels <- plot_other_dims
    }
    else {
        axis_labels <- c("UMAP1", "UMAP2")
    }
    
    if (!is.null(thresh)) {
        if (transposed) {
            gene_data <- quanThresh(data[, gene], thresh)
        }
        else {
            gene_data <- quanThresh(data[gene, ], thresh)
        }
    }
    else {
        if (transposed) {
            gene_data <- data[, gene]
        }
        else {
            gene_data <- data[gene, ]
        }
    }
    gene_frame <- cbind.data.frame(gene = gene_data, umap_frame)
    plt <- gene_frame %>% sample_frac %>% ggplot(aes(x = get(axis_labels[1]), 
        y = get(axis_labels[2]), z = gene))
    if (do_points) {
        plt <- plt + stat_summary_hex()
    }
    if (do_density) {
        plt <- plt + geom_density_2d()
    }
    plt <- plt + theme_classic() + labs(x = axis_labels[1], y = axis_labels[2], 
        fill = "", title = gene) + scale_fill_gradient(low = "snow2", 
        high = "red")
    return(plt)
}

plot_gene_general <- function (data, umap_frame, gene, pt_size = 0.1, pch = 19, thresh = NULL, 
    do_points = TRUE, do_density = FALSE, title = NULL, transposed = FALSE, plot_other_dims = NULL) 
{
    if (!is.null(plot_other_dims)) {
        axis_labels <- plot_other_dims
    }
    else {
        axis_labels <- c("UMAP1", "UMAP2")
    }
    
    if (!is.null(thresh)) {
        if (transposed) {
            gene_data <- quanThresh(data[, gene], thresh)
        }
        else {
            gene_data <- quanThresh(data[gene, ], thresh)
        }
    }
    else {
        if (transposed) {
            gene_data <- data[, gene]
        }
        else {
            gene_data <- data[gene, ]
        }
    }
    gene_frame <- cbind.data.frame(gene = gene_data, umap_frame)
    plt <- gene_frame %>% sample_frac %>% ggplot(aes(x = get(axis_labels[1]), 
        y = get(axis_labels[2]), col = gene))
    if (do_points) {
        plt <- plt + geom_point(size = pt_size)
    }
    if (do_density) {
        plt <- plt + geom_density_2d()
    }
    plt <- plt + theme_classic() + labs(x = axis_labels[1], y = axis_labels[2], 
        col = "", title = gene) + scale_color_gradient(low = "snow2", 
        high = "red")
    return(plt)
}

add_gene_violinPlot <- function(plt, do_mean = TRUE, do_median = FALSE, stat_size = 2, stat_col = 'black',
                                x_angle = 45, x_hjust = 1, title = NULL, x_title = NULL, y_title = NULL){
    plt <- plt + geom_violin() + theme_classic()
    if(do_mean){
        plt <- plt + stat_summary(fun.y = mean, geom = 'point', size = stat_size, color = stat_col)
    }
    if(do_median){
        plt <- plt + stat_summary(fun.y = median, geom = 'point', size = stat_size, color = stat_col)
    }    
    plt <- plt + theme(axis.text.x = element_text(angle = x_angle, hjust = x_hjust))
    plt <- plt + labs(title = title, x = x_title, y = y_title)
    return(plt)
}


# return all elements (str) that are present (or not present) in a vector. 
grepElements <- function(str, vec, get = 'TRUE'){
    if(get == 'TRUE'){
        return(vec[which(grepl(str, vec) == 'TRUE')])
    } else{
        return(vec[which(grepl(str, vec) == 'FALSE')])
    }
}

# sort hierarchical clustering by distance (see Kam's pheatmap tutorial)
sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))

# function to rotate colnames in a pheatmap by 45 degrees
draw_colnames_45 <- function (coln, gaps, ...) {
    coord <- pheatmap:::find_coordinates(length(coln), gaps)
    x     <- coord$coord - 0.5 * coord$size
    res   <- grid::textGrob(
      coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"),
      vjust = 0.75, hjust = 1, rot = 45, gp = grid::gpar(...)
    )
    return(res)
}
assignInNamespace(
  x = "draw_colnames",
  value = "draw_colnames_45",
  ns = asNamespace("pheatmap")
)

# accessory functions to presto
top_auc_pos_markers <- function (markers_table, clus, minFC = 0, max_pval = 0.05) 
{
    a <- markers_table[which(markers_table$group == clus), ]
    b <- a[which(a$padj < max_pval),]
    c <- b[order(b$auc, decreasing = T),]
    d <- c[which(c$logFC > minFC), ]
    return(d)
}

top_clus_pos_markers <- function(markers_table, minFC = 0, group = NULL, max_pval = 0.05){
    a <- markers_table[order(markers_table$auc, decreasing = T),]
    b <- a[which(a$padj < max_pval),]
    c <- b[which(b$logFC > minFC),]
    if(!is.null(group)){
        return(c[which(c$group == group),])
    } else{
        return(c)
    }
}

wilcox.nest <- function(X, y, y_nest = NULL, groups_use = NULL, verbose = TRUE, ...){
#     if (ncol(X) != length(y)) 
#         stop("number of columns of X does not\n match length of y")

    if(!is.null(y_nest) & is.null(groups_use)){
        res_list <- lapply(unique(y), function(group){
            idx_use <- which(y %in% intersect(group, y))
            y <- y_nest[idx_use]
            X <- X[, idx_use]
            res <- wilcoxauc(X, y)
        })
        names(res_list) <- paste0('Group_', unique(y))
        return(res_list)
    } else if(!is.null(y_nest) & !is.null(groups_use)){
        res_list <- lapply(groups_use, function(group){
            idx_use <- which(y %in% intersect(group, y))
            y <- y_nest[idx_use]
            X <- X[, idx_use]
            res <- wilcoxauc(X, y)
        })
        names(res_list) <- paste0('Group_', unique(groups_use))
        return(res_list)
    }
}

wilcox.nestSeurat <- function (X, group_by = NULL, group_by_nest = NULL, assay = "data", 
    groups_use = NULL, seurat_assay = "RNA", ...) 
{
    X_matrix <- Seurat::GetAssayData(X, assay = seurat_assay, 
        slot = assay)
    if (is.null(group_by)) {
        y <- Seurat::Idents(X)
    } else {
        y <- Seurat::FetchData(X, group_by) %>% unlist %>% as.character()
    }
    if (is.null(group_by_nest)) {
        group_by_nest <- Seurat::FetchData(X, "orig.ident") %>% 
            unlist %>% as.character()
    } else {
        group_by_nest <- Seurat::FetchData(X, group_by_nest) %>% 
            unlist %>% as.character()
    }
    wilcox.nest(X_matrix, y, y_nest = group_by_nest, groups_use = NULL)
}

genes_heatmapMat <- function(data, meta, genes){
    a <- data.frame('sub_type' = as.factor(as.character(meta[, 'sub_type'])),
                    data[, unique(genes)])
    b <- a %>% group_by(sub_type) %>%
        summarize_all(mean) %>% column_to_rownames('sub_type')
}

genes_heatmapDisease <- function(data, meta, genes){
    a <- data.frame('binary_polyp' = as.factor(as.character(meta[, 'binary_polyp'])),
                   'sub_type' = as.factor(as.character(meta[, 'sub_type'])),
                    data[, unique(genes)])
    a[, 'groups'] <- paste0(a[, 'binary_polyp'], '_', a[, 'sub_type'])
    b <- a %>% select(-c('binary_polyp', 'sub_type')) %>% group_by(groups) %>%
        summarize_all(mean) %>% column_to_rownames('groups')
    return(b)
}

get_group_numbers <- function (data, groups, nest_groups, prefix = "clus") 
{
    prop_table <- sapply(sort(unique(data[, groups])), function(x) {
        idx <- which(data[, groups] == x)
        return(table(data[, nest_groups][idx]))
    })
    prop_table %<>% t
    return(prop_table)
#     final_table <- cbind(groups = paste0(prefix, sort(unique(data[, groups]))), prop_table)
#     rownames(final_table) <- paste0(prefix, sort(unique(data[, groups])))
#     return(final_table)
}

addGroupToMelt <- function(table, group, start_ind = 1){
    tbl <- table %>% melt
    if(start_ind == 0){
        tbl$group <- factor(paste0(start_ind:((group %>% unique %>% length) - 1)),
                            levels = mixedsort(paste0(start_ind:((group %>% unique %>% length) - 1))))
    } else{
        tbl$group <- factor(paste0(start_ind:(group %>% unique %>% length)),
                            levels = mixedsort(paste0(start_ind:(group %>% unique %>% length))))
    }
    return(tbl)
}

generate_trajectories <-function(input_cell_coordinates){
    # input_cell_coordinates should have dimensions Cells x Dimensions
    # Returns 2D coordinates for cells on trajectory graph
    n_cells = nrow(input_cell_coordinates)
    n_centers = round(2 * 100 * log(n_cells) / (log(n_cells) + log(100)))
    ddrtree_output = DDRTree(X = t(input_cell_coordinates), ncenter = n_centers,
                             # Default parameters
                             dimensions = 2, 
                             initial_method = NULL, 
                             maxIter = 20, 
                             sigma = 0.001, # bandwidth parameter
                             lambda = NULL, # regularization parameter for inverse graph embedding
                             param.gamma = 10, # regularization parameter for k-means
                             tol = 0.001, # relative objective difference
                             verbose = F)
    # output is a list with W, Z, stree, Y, history
    # W is the orthogonal set of d (dimensions), 
    # linear basis vector Z is the reduced dimension space
    # stree is the smooth tree graph embedded in the low dimension space 
    # Y represents latent points as the center of Z
    trajectory_points = t(ddrtree_output$Z)
    colnames(trajectory_points) <- c('DDRTree1', 'DDRTree2')
    
    principal_curve_output <- principal_curve(trajectory_points)
    trajectory_points = data.frame(trajectory_points)
    trajectory_points$pseudo <- principal_curve_output$lambda
    # Fits a principal curve which describes a smooth curve that passes through the middle 
    # of the data x in an orthogonal sense. This curve is a nonparametric generalization of 
    # a linear principal component.
    
    return(trajectory_points)
}

generate_trajectories_custom <-function(input_cell_coordinates, sigma_val, lambda_val, gamma_val){
    # input_cell_coordinates should have dimensions Cells x Dimensions
    # Returns 2D coordinates for cells on trajectory graph
    n_cells = nrow(input_cell_coordinates)
    n_centers = round(2 * 100 * log(n_cells) / (log(n_cells) + log(100)))
    ddrtree_output = DDRTree(X = t(input_cell_coordinates), ncenter = n_centers,
                             # Default parameters
                             dimensions = 2, 
                             initial_method = NULL, 
                             maxIter = 20, 
                             sigma = sigma_val, # bandwidth parameter
                             lambda = lambda_val, # regularization parameter for inverse graph embedding
                             param.gamma = gamma_val, # regularization parameter for k-means
                             tol = 0.001, # relative objective difference
                             verbose = F)
    # output is a list with W, Z, stree, Y, history
    # W is the orthogonal set of d (dimensions), 
    # linear basis vector Z is the reduced dimension space
    # stree is the smooth tree graph embedded in the low dimension space 
    # Y represents latent points as the center of Z
    trajectory_points = t(ddrtree_output$Z)
    colnames(trajectory_points) <- c('DDRTree1', 'DDRTree2')
    
    principal_curve_output <- principal_curve(trajectory_points)
    trajectory_points = data.frame(trajectory_points)
    trajectory_points$pseudo <- principal_curve_output$lambda
    # Fits a principal curve which describes a smooth curve that passes through the middle 
    # of the data x in an orthogonal sense. This curve is a nonparametric generalization of 
    # a linear principal component.
    
    return(trajectory_points)
}

plot_pseudotime <-function(trajectory_points, col = NULL, pt_size = 0.2, title = '', scale_virid = FALSE,
                           override_legend = TRUE, legend_size = 4, legend_title = 'pseudo', plot_gene = FALSE,
                           do_facet = FALSE, do_density = FALSE, do_points = TRUE, plot_other_dims = NULL){
    if(!is.null(plot_other_dims)){
        axis_labels <- plot_other_dims
    } else{
        axis_labels <- c("DDRTree1", "DDRTree2")
    }
    trajectory_points <- cbind(trajectory_points, col)
    plt <- trajectory_points %>% sample_frac %>% 
        ggplot(aes(x = get(axis_labels[1]), y = get(axis_labels[2]), col = col))
    if(scale_virid){
         plt <- plt + scale_color_viridis(option = 'plasma')
    }
    if(plot_gene){
         plt <- plt + scale_color_gradient(low = "snow2", high = "royalblue") + labs(col = col)
    }
    if(do_points){
        plt <- plt + geom_point(size = pt_size)
    }
    if(do_density){
        plt <- plt + geom_density_2d()
    }
    plt <- plt +
        labs(title = title, col = legend_title, x = axis_labels[1], y = axis_labels[2]) +
        theme_classic()
    if(override_legend){
        plt <- plt + guides(color = guide_legend(override.aes = list(stroke = 1, 
            alpha = 1, shape = 19, size = legend_size)))
    }
    if(do_facet){
        plt <- plt + facet_wrap(~col)
    }
    return(plt)
}

insertElems = function(vect, pos, elems) {
  l = length(vect)
  j = 0
  for (i in 1:length(pos)){
    if (pos[i]==1)
      vect = c(elems[j+1], vect)
    else if (pos[i] == length(vect)+1)
      vect = c(vect, elems[j+1])
    else
      vect = c(vect[1:(pos[i]-1+j)], elems[j+1], vect[(pos[i]+j):(l+j)])
    j = j+1
  }
  return(vect)
}

induceFC <- function(baseline_f, clus = 1, fc = 2){
    baseline_f <- baseline_f %>% prop.table
    clusxfc <- baseline_f[clus] * fc
    otherclusters <- (baseline_f[-clus]/sum(baseline_f[-clus])) * (1 - clusxfc)
    return(insertElems(otherclusters, pos = clus, elems = clusxfc))
}

# Add data to workbook for exporting to excel spreadsheet
addSheetToWorkbook <- function(wb, data, name){
    addWorksheet(wb, name)
    writeData(wb, name, data)
}

calcBinomCI <- function(p, n, z = 1.96) {
   return(z * sqrt((p * (1 - p))/n))
}

# Get the aesthetics of a specific ggplot. Helpful to recover colors and make them consistent across ggplots
getGGPlot_aes <- function(ggPlot, aes = "colour"){
    build <- ggplot_build(ggPlot)
    return(unique(build$data[[1]][[aes]]))
}

AnnoyBuildObj <- function(data.use, metric = "euclidean", n.trees = 50){
    l <- ncol(x = data.use)
    annoyObj <- switch(EXPR = metric, 
                       euclidean = new(Class = RcppAnnoy::AnnoyEuclidean, l), 
                       cosine = new(Class = RcppAnnoy::AnnoyAngular, l), 
                       manhattan = new(Class = RcppAnnoy::AnnoyManhattan, l), 
                       hamming = new(Class = RcppAnnoy::AnnoyHamming, l), 
                       stop("Enter valid distance metric"))
    for(i in seq(nrow(x = data.use))){
        annoyObj$addItem(i - 1, data.use[i, ])
    }
    annoyObj$build(n.trees)
    return(annoyObj)
}

AnnoyGetNN <- function(obj, query, k.param = 30, search.k = -1, get.distance = TRUE, mc.cores = 1){
    n <- nrow(x = query)
    nn.idx <- mclapply(X = 1:n, function(x){
        nn <- obj$getNNsByVectorList(query[x, ], k.param, search.k, get.distance)
        list(nn$item + 1, nn$distance)
    }, mc.cores = mc.cores)
    idx <- do.call(rbind, lapply(nn.idx, "[[", 1))
    if(get.distance){
        dist <- do.call(rbind, lapply(nn.idx, "[[", 2))
    } else{
        dist <- matrix(nrow = n, ncol = k.param)
    }
    return(list(nn.idx = idx, nn.dists = dist))
}

BuildSNN <- function(data.use, k.param = 30, prune.SNN = 1/15, nn.eps = 0, method = "annoy", n.trees = 50,
                     search.k = -1, metric = "euclidean", mc.cores = 2){
    switch(EXPR = method, rann = {
        knn <- RANN::nn2(data = data.use, k = k.param, searchtype = "standard", eps = nn.eps)
    }, annoy = {
        annoyObj <- AnnoyBuildObj(data.use = data.use, metric = metric, n.trees = n.trees)
        knn <- AnnoyGetNN(obj = annoyObj, query = data.use, k.param = k.param, 
                          search.k = -1, get.distance = TRUE, mc.cores = mc.cores)
    })
    nn.ranked <- knn$nn.idx
    snn_res <- Seurat:::ComputeSNN(nn_ranked = nn.ranked, prune = prune.SNN)
    rownames(snn_res) <- row.names(data.use)
    colnames(snn_res) <- row.names(data.use)
    return(snn_res)
}

theme_gy <- function (base_size = 23) 
{
    theme_bw(base_size = base_size) + 
    theme(plot.title = element_text(hjust = 0.5, size = (base_size + 2)), 
          legend.title = element_text(face = "italic"),
          axis.title = element_text(size = 26),
          axis.text = element_text(size = 26))
}

read10x_mtx <- function(run, suffix, min_counts=1) {
    #min_counts initially 100
    barcode.loc <- list.files(run, pattern = 'barcodes.tsv(.gz)?', full.names = TRUE)
    gene.loc <- list.files(run, pattern = 'features.tsv(.gz)?', full.names = TRUE)
    matrix.loc <- list.files(run, pattern = 'matrix.mtx(.gz)?', full.names = TRUE)
    
    data <- readMM(file = matrix.loc) %>% as("CsparseMatrix")# %>% Matrix::t()
    cell.names <- readLines(barcode.loc)
    cell.names <- gsub("-1$", "", cell.names)    
    
    if (!missing(suffix)) {
        cell.names <- paste(cell.names, suffix, sep = "_")
    }
    
    genes = fread(gene.loc, header = FALSE)
    if (ncol(genes)>1){
        gene.names <- genes$V2
    } else {
        gene.names <- genes$V1
    }    
    
    row.names(data) <- gene.names
    colnames(data) <- cell.names

    data <- as(data, "CsparseMatrix")
    data <- data[, Matrix::colSums(data) >= min_counts]
    data <- data[which(!is.na(row.names(data))), ]
    data <- as(sumOverRowNames(data), "CsparseMatrix")
    return(data)
}

read10x_kallisto = function (run, suffix, min_counts = 1) {
    barcode.loc <- list.files(run, pattern = "cells_x_genes.barcodes.txt(.gz)?", 
        full.names = TRUE)
    gene.loc <- list.files(run, pattern = "cells_x_genes.genes.txt(.gz)?", 
        full.names = TRUE)
    matrix.loc <- list.files(run, pattern = "cells_x_genes.mtx(.gz)?", 
        full.names = TRUE)
    data <- readMM(file = matrix.loc) %>% as("dgCMatrix") %>% t
    cell.names <- readLines(barcode.loc)
    cell.names <- gsub("-1$", "", cell.names)
    if (!missing(suffix)) {
        cell.names <- paste(cell.names, suffix, sep = "_")
    }
    gene.names <- fread(gene.loc, header = FALSE)$V1
    row.names(data) <- gene.names
    colnames(data) <- cell.names
    data <- as(data, "dgCMatrix")
    data <- data[, Matrix::colSums(data) >= min_counts]
    data <- data[which(!is.na(row.names(data))), ]
    data <- as(sumOverRowNames(data), "dgCMatrix")
    return(data)
}

process_rna = function(exprs, meta, fields, filters, topn = 3000, gene_exp_pct = 0.05){
    ### Filter metadata so fields are met
    if (!missing(fields)){
        meta = meta[colSums((meta[fields] %>% t) == filters) == length(filters), ]
    }
    
    # Filter out genes not expressed in at least some % of cells
    exprs = exprs[rowSums(exprs[, rownames(meta)[rownames(meta) %in% colnames(exprs)]] > 0) 
                  >= floor(nrow(meta)*gene_exp_pct), ]

    norm_exprs = exprs[, rownames(meta)[rownames(meta) %in% 
                                              colnames(exprs)]] %>% singlecellmethods::normalizeData(method = 'log')
    var_genes <- vargenes_vst(norm_exprs, #meta[colnames(norm_exprs), ]$plate, 
                              topn = topn)
    length(var_genes %>% unique)
    
    scale_exprs <- norm_exprs[var_genes, ] %>% singlecellmethods::scaleData()
    
    return(scale_exprs)
    
}

process_adt = function(exprs, meta, fields, filters, topn = 3000){
    ### Filter metadata so fields are met
    if (!missing(fields)){
        meta = meta[colSums((meta[fields] %>% t) == filters) == length(filters), ]
    }
    
    norm_exprs = exprs[, rownames(meta)[rownames(meta) %in% 
                        colnames(exprs)]] %>% singlecellmethods::normalizeData(method = 'cellCLR')
    var_genes = rownames(norm_exprs)
    scale_exprs <- norm_exprs[var_genes, ] %>% singlecellmethods::scaleData()
    
    return(scale_exprs)
    
}

normalize_rna = function(exprs, meta, fields, filters, scale = FALSE){
    ### Filter metadata so fields are met
    if (!missing(fields)){
        meta = meta[colSums((meta[fields] %>% t) == filters) == length(filters), ]
    }
    
    norm_exprs = exprs[, rownames(meta)[rownames(meta) %in% 
                    colnames(exprs)]] %>% singlecellmethods::normalizeData(method = 'log')
    
    if (scale){
        scale_exprs <- norm_exprs %>% singlecellmethods::scaleData() %>% as.data.frame
    }else{
        scale_exprs <- norm_exprs %>% as.data.frame #norm_exprs[var_genes, ] #%>% singlecellmethods::scaleData()
    }
    return(scale_exprs)
}

#Alleles Function for Import - Updated- Fixed Naming to Automate and Changed to Capture the top 10 alleles per cell. 
AllelesFunction <- 
function (filepath, number_of_alleles = 10) #Filepath is the folder with all the allele tables. 
{
    filepath <- sub("/$", "", filepath)
    files <- list.files(path = filepath, pattern = "Alleles_frequency_table.txt$", full.names = T, 
        recursive = T) 
        
    matrix <- rbindlist(lapply(files, function(f) {
    # Extract Plate_Well ID from the file path
    plate_id <- str_extract(f, "[^/]+(?=/CRISPRessoBatch_on_batch/CRISPResso_on_[^/]+)") %>%
                paste0("_", str_extract(f, "(?<=CRISPResso_on_)[^/]+"))
    
        # Read file and add Plate_Well column
        mat_file <- fread(f, nrows = number_of_alleles, select = c(1, 5, 6, 7, 8, 9))
        mat_file[, Plate_Well := plate_id]  # Assign the extracted Plate_Well ID

        return(mat_file)  # Return the modified table
    }), use.names = TRUE, fill = TRUE)  # Bind all tables together

    matrix = cbind(matrix[, (ncol(matrix)):(ncol(matrix))], matrix[, 1:(ncol(matrix)-1)])
    
    return(matrix)
}

##Filtering Function and spreading out the two alleles. 
FilterAlleles <- function (matrix, n=10, p=40) 
{
    test <- matrix %>% dplyr::filter(`%Reads` > n) %>% dplyr::filter(`#Reads` > 
        p) %>% mutate(Allele = ifelse(duplicated(Plate_Well) == 
        F, "Allele1", "Allele2")) %>% mutate(AlleleReads = ifelse(duplicated(Plate_Well) == 
        F, "Allele1read", "Allele2read")) %>% spread(AlleleReads, 
        `#Reads`) %>% spread(Allele, Aligned_Sequence) %>% select(Plate_Well, 
        Allele1, Allele1read, Allele2, Allele2read, Reference)
    Return <- left_join(drop_na(test[, c(1, 2, 3, 6)]), drop_na(test[, 
        c(1, 4, 5, 6)])) %>% mutate(OnlyOneAlleleRecovered = ifelse(is.na(Allele2) == 
        T, T, F)) %>% mutate(Allele2 = ifelse(is.na(Allele2) == 
        T, Allele1, Allele2)) %>% mutate(Allele2read = ifelse(is.na(Allele2read) == 
        T, Allele1read, Allele2read))
    return(Return)
}

##Wait, modify the genotypes function
GenotypesFunction <- 
function (Return) 
{
    ## define a unique set of a_b sorted on (rev freq)
    test <- Return %>% filter(OnlyOneAlleleRecovered == F) %>% 
    mutate(Freq = paste0(Allele1,"_",Allele2)) %>% add_count(Freq) %>% arrange(n) %>%
    dplyr::select(Allele1, Allele2, Freq, n) %>% unique
   
    # Fix matrix. This can be fixed to not scan for HETS pretty quickly. 
    #For each combination of Allele 1 and 2
    for(i in (test$Freq)) {
    a <- str_split_fixed(i,2, pattern = "_")[1] #pull out allele1
    b <- str_split_fixed(i,2, pattern = "_")[2] #pull out allele2
    
    #Scan matrix of alleles one row at a time
    for(j in 1:dim(Return)[1]) { 
        
    if(Return$Allele1[j] == b & Return$Allele2[j] == a) { #if the allele combo you are looking for is found in reverse. 
    Return$Allele1[j] <- a # fix it
    Return$Allele2[j] <- b # fix it
    }}
    
    }
    
    allele1_mat <- matrix(str_split(Return$Allele1, "", simplify = TRUE), 
        nrow = nrow(Return), byrow = FALSE)
    allele2_mat <- matrix(str_split(Return$Allele2, "", simplify = TRUE), 
        nrow = nrow(Return), byrow = FALSE)
    ref_mat <- matrix(str_split(Return$Reference, "", simplify = TRUE), 
        nrow = nrow(Return), byrow = FALSE)
    Ref <- matrix(paste0(ref_mat, ref_mat), nrow(allele1_mat), 
        ncol(allele1_mat))
    Full <- matrix(paste0(allele1_mat, allele2_mat), nrow(allele1_mat), 
        ncol(allele1_mat))
    #Full <- Full %>% apply(2, function(x) {
     #   if_else(x == "TA", "AT", if_else(x == "GA", "AG", if_else(x == 
      #      "CA", "AC", if_else(x == "TC", "CT", if_else(x == 
       #     "GC", "CG", if_else(x == "TG", "GT", if_else(x == 
        #    "A-", "-A", if_else(x == "T-", "-T", if_else(x == 
         #   "G-", "-G", if_else(x == "C-", "-C", x))))))))))
    #})
    Full <- apply(Full, 1, paste, collapse = "") %>% as.matrix()
    Ref <- apply(Ref, 1, paste, collapse = "") %>% as.matrix
    Output <- as_tibble(cbind(Return$Barcode_DNA, Return$Well_ID, 
        Full, Ref))
    colnames(Output) <- c("Barcode_DNA", "Well", "Alleles", "Reference")
    genotypes_df <- Output %>% dplyr::select(Alleles) %>% unique() %>% 
        tibble::rownames_to_column("genotype")
    Output <- left_join(Output, genotypes_df)
    return(Output)
}

#GenotypeHeatmap2
GenotypeHeatmap2 <- function (matrix, y) 
{
    CRISPRessoDNA_heatmap <- matrix
    CRISPRessoDNA_heatmap <- CRISPRessoDNA_heatmap %>% add_count(genotype) %>% 
        dplyr::rename(count = n)
    CRISPRessoDNA_heatmap <- dplyr::filter(CRISPRessoDNA_heatmap, 
        count > y) %>% arrange(desc(count))

 ALLELES <- apply(matrix(CRISPRessoDNA_heatmap$Alleles), 1, 
        function(x) {
            gsub("(.{2})", "\\1 ", x)
        }) %>% str_split_fixed(pattern = " ", n = str_length(CRISPRessoDNA_heatmap$Alleles[1])/2)
    REF <- apply(matrix(CRISPRessoDNA_heatmap$Reference), 1, 
        function(x) {
            gsub("(.{2})", "\\1 ", x)
        }) %>% str_split_fixed(pattern = " ", n = str_length(CRISPRessoDNA_heatmap$Alleles[1])/2)
x <- copy(ALLELES)
    x[ALLELES == REF] <- "R"
    y <- copy(ALLELES)
    x <- x %>% unique %>% as.data.frame(stringsAsFactors = FALSE) %>% tidyr::gather(position, 
        value) %>% dplyr::mutate(position = gsub("V", 
        "", position))
    rownames(y) <- CRISPRessoDNA_heatmap$genotype

 y <- y %>% unique %>% as.data.frame(stringsAsFactors = FALSE) %>% 
        rownames_to_column("genotypes") %>% tidyr::gather(position, 
        value, -genotypes) %>% dplyr::mutate(position = gsub("V", 
        "", position))
 x <- mutate(x, genotypes = y$genotypes)
   y <- y %>% dplyr::rename(label = value)
    x <- x %>% cbind(y$label)
    x <- mutate(x, Genotype = paste0("G", genotypes))
    z <- y$label %>% str_split_fixed(pattern = "", n = 2)
    colnames(z) <- c("A", "B")
    df <- cbind(z, x) %>% pivot_longer(cols = -c(position, value, 
        `y$label`, Genotype, genotypes), values_to = "allele_value", 
        names_to = "allele")
    return(df)
}

pca_umap = function(scale_exprs, meta, fields, filters, harmony_var, npc = 25){
    ### Filter metadata so fields are met
    if (!missing(fields)){
        meta = meta[colSums((meta[fields] %>% t) == filters) == length(filters), ]
    }
    
    pca_res = prcomp_irlba(scale_exprs %>% t, npc)
    rownames(pca_res$x) = colnames(scale_exprs)
    rownames(pca_res$rotation) = rownames(scale_exprs)
    pcs = pca_res$x

    fig.size(5, 8)
    p = fviz_screeplot(pca_res, #addlabels = TRUE,
                   ncp = npc) +
        theme_g()+
        ggtitle('Percentage explained variance by PC')+
        xlab('PC')+
        ylab('% explained variances')  
    print(p)
    
    
    if (!missing(harmony_var)){
        set.seed(200)
        umap_input <- RunHarmony(pcs[, 1:10], meta[rownames(pcs), ], 
                                     harmony_var, theta = rep(1, length(harmony_var)), 
                               plot_convergence = TRUE, max.iter.harmony = 10, epsilon.cluster = -Inf, 
                                 epsilon.harmony = -Inf, 
                                 max.iter.cluster = 10, do_pca = F, verbose = T)
    } else {umap_input = pcs[, 1:10]}
    
    
    set.seed(200)
    umap_res = umap(umap_input)
    umap_df = data.frame(umap_res)

    colnames(umap_df) = c('UMAP1', 'UMAP2')
    meta_umap = cbind(meta[rownames(umap_df), ], umap_df, umap_input[rownames(umap_df), ])

    return(meta_umap)
    
}

plot_exprs = function(scale_exprs, meta_umap, markers, fig_size = c(6, 8)){
    
    exprs_umap = cbind(meta_umap, scale_exprs[, rownames(meta_umap)] %>% t)
    
    
    fig.size(fig_size[1], fig_size[2])
    for (marker in markers){
        print(ggplot(exprs_umap)+
        geom_point(aes(x= UMAP1, y = UMAP2, col = eval(parse(text = marker))), size = 3)+
        labs(col = marker)+
        theme_g(30)+
        scale_color_viridis())
    }
}

do_cluster = function(meta_umap, resolution_list = seq(0.2, 1.2, 0.2)){
    pca_res = meta_umap %>% select(starts_with('PC'))
    
    set.seed(200)
    snn_ref <- FindNeighbors(as.matrix(pca_res))
    # resolution_list <- c(0.2, 0.4, 0.6, 0.8, 1, 1.2)
    
    clusters <- lapply(resolution_list, function(x) {
        FindClusters(snn_ref$snn, resolution = x, random.seed = 100)
    }) %>% bind_cols
    
    meta_umap_clust = cbind(meta_umap, clusters[rownames(meta_umap), ])
    return(meta_umap_clust)
    
    
}

### Hierarchially cluster rows of a matrix using pheatmap within groups, returning row order
hclust_groups = function(mat, groups, cluster_cols = FALSE){
    order = c()
    
    for (group in unique(groups)){
        mat_group = mat[groups==group, ]
        group_order = pheatmap(mat_group, cluster_cols = cluster_cols)$tree_row$order
        order = c(order, rownames(mat_group)[group_order])
    }
    
    return(order)
}

library(dendsort)
sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))

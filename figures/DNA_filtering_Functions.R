# set figure size
fig.size <- function(height = 4, width = 4){
    options(repr.plot.height = height, repr.plot.width = width)
}

theme_gy <- function (base_size = 23) 
{
    theme_bw(base_size = base_size) + 
    theme(plot.title = element_text(hjust = 0.5, size = (base_size + 2)), 
          legend.title = element_text(face = "italic"),
          axis.title = element_text(size = 26),
          axis.text = element_text(size = 26))
}

sumOverRowNames <- function(X) {
    name_factors <- factor(row.names(X))
    res <- presto::sumGroups(X, name_factors)
    row.names(res) <- levels(name_factors)#[1:nrow(res)]
    colnames(res) <- colnames(X)
    return(res)
    }

read10x_mtx <- function(run, suffix, min_counts=1) {
    #min_counts initially 100
    barcode.loc <- list.files(run, pattern = 'barcodes.tsv(.gz)?', full.names = TRUE)
    gene.loc <- list.files(run, pattern = 'features.tsv(.gz)?', full.names = TRUE)
    matrix.loc <- list.files(run, pattern = 'matrix.mtx(.gz)?', full.names = TRUE)
    
    data <- readMM(file = matrix.loc) %>% as("dgCMatrix")# %>% Matrix::t()
    cell.names <- readLines(barcode.loc)
    cell.names <- gsub("-1$", "", cell.names)    
    
    if (!missing(suffix)) {
        cell.names <- paste(cell.names, suffix, sep = "_")
    }
    
    gene.names <- fread(gene.loc, header = FALSE)$V2
    row.names(data) <- gene.names
    colnames(data) <- cell.names

    data <- as(data, "dgCMatrix")
    data <- data[, Matrix::colSums(data) >= min_counts]
    data <- data[which(!is.na(row.names(data))), ]
    data <- as(sumOverRowNames(data), "dgCMatrix")
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
        theme_gy()+
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
        theme_gy(30)+
        scale_color_viridis())
    }
}

do_cluster = function(meta_umap, resolution_list = seq(0.2, 1.2, 0.2)){
    pca_res = meta_umap %>% dplyr::select(starts_with('PC'))
    
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

## Function for calling genotypes. Input aligned and trimmed allele table. 
Genotyping_Cells <- function(alleles, cell_col) { 
    matrix <- alleles
    Genotypes <- matrix$Aligned_Sequence %>% unique
    
    #Collapse the values of all alleles detected into one string.
    matrix_genotypes <- matrix %>% group_by(!!sym(cell_col)) %>% arrange(Aligned_Sequence) %>% 
            summarize(genotype = str_c(Aligned_Sequence, collapse ="_"))
    #Add empty allelic genotype column
    matrix_genotypes<- mutate(matrix_genotypes, AllelicGenotype = "")
    
    #Add genotype ID for later deconvolution
    for (i in 1:length(Genotypes)){      
        matrix_genotypes <- 
            mutate(matrix_genotypes, 
            AllelicGenotype = str_c(AllelicGenotype, ifelse(grepl(genotype, pattern = fixed(Genotypes[i])), LETTERS[i],"")))}
    return(matrix_genotypes)
}

## Function for plotting Alleles
Plotting_Alleles <- function(Filter_Allele_Tibble, ref = NULL) { 
    if(!missing(ref)) {
        Filter_Allele_Tibble[['Reference']] = ref
    }

    #Identify unique alleles
    if ('Reference' %in% colnames(Filter_Allele_Tibble)){
        UniqueAlleles <- select(Filter_Allele_Tibble, Aligned_Sequence, Reference) %>% unique() %>% as_tibble()

    } else{
        UniqueAlleles <- select(Filter_Allele_Tibble, Aligned_Sequence) %>% unique() %>% as_tibble()
    }

    #Separate bases across the amplicon
    ALLELES <- apply(matrix(UniqueAlleles$Aligned_Sequence), 1, 
        function(x) {
            gsub("(.{1})", "\\1 ", x)
        }) %>% str_split_fixed(pattern = " ", n = str_length(UniqueAlleles$Aligned_Sequence[1]))

        if ('Reference' %in% colnames(Filter_Allele_Tibble)){

        REF <- apply(matrix(UniqueAlleles$Reference), 1, 
        function(x) {
            gsub("(.{1})", "\\1 ", x)
        }) %>% str_split_fixed(pattern = " ", n = str_length(UniqueAlleles$Reference[1]))

    }

    x <- copy(ALLELES)

    #Identify reference nucleotides
    if ('Reference' %in% colnames(Filter_Allele_Tibble)){
        x[ALLELES == REF] <- "R"
        }
    y <- copy(ALLELES)

    #Collapse into a dataframe. 
    x <- x %>% unique %>% as.data.frame(stringsAsFactors = FALSE) %>% 
        tidyr::gather(position, value) %>% dplyr::mutate(position = gsub("V", 
        "", position))
    y <- y %>% unique %>% as.data.frame(stringsAsFactors = FALSE) %>% 
        rownames_to_column("genotypes") %>% tidyr::gather(position, 
        value, -genotypes) %>% dplyr::mutate(position = gsub("V", 
        "", position))
    x <- mutate(x, genotypes = y$genotypes)
    y <- y %>% dplyr::rename(label = value)
    x <- x %>% cbind(y$label)
    x <- mutate(x, genotypes = LETTERS[as.numeric(genotypes)])

    #Plotting the data
    fig.size(10,10)
       g<- x %>% #the number is the filter on genotype 
    ggplot(aes(x = reorder(position, as.integer(position)), 
               y = 1, fill = `value`)) + 
        geom_tile(aes(color = value), width=1, height=1) +
            geom_text(data = dplyr::filter(x),
                      aes(label = `y$label`), angle = 0, size = (4)) + 
            geom_text(data = dplyr::filter(x, value != "R"),
                      aes(label = `y$label`), angle = 0, size = (4))+
            scale_fill_manual(values = 
                        c(R = 'white',
                          C = "#5194ed", 
                          T = "#fdb462", 
                          G = "#7fc97f", 
                          A = "#ef3b2c",
                         `-` ="grey80"))+ 
            scale_color_manual(values = 
                       c(R = 'grey'))+

            ggtitle("") + 
            theme_gy()+
            xlab("")+
            theme(aspect.ratio = 0.1)+ 
            ylab("") +
            theme(axis.text.y = element_text(size = 0, angle = 0, hjust = .5, vjust = .5),
              axis.text.x = element_text(size = 0),
            axis.title.x = element_text(size = 20, angle = 0, hjust = .5, vjust = .5),
            axis.title.y = element_text(size = 20, angle = 90, hjust = .5, vjust = .5),
              axis.ticks.y = element_blank(),
             legend.position = "none") + 
        scale_x_discrete(breaks=NULL) + 
        facet_grid(factor((genotypes),)~.)+ 
    theme(
      strip.background = element_blank(),
      strip.text.x = element_blank(), 
        panel.border=element_blank()

    ) + 
        annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1) + 
        annotate(geom = 'segment', y = -Inf, yend = -Inf, color = 'black', x = -Inf, xend = Inf, size = 1)

    return(g)

}


Filtering_Cells_Read <- function(alleles, plate_col, cell_col) { 
    mat <- alleles %>% mutate(TotalReads = `#Reads`/(`%Reads`/100))
    Plates <- unique(mat[[plate_col]])
    
    ReadFilters <- lapply(Plates, function(Plate) {
        data <- (mat[mat[[plate_col]]==Plate,] %>% group_by(!!sym(cell_col)) %>% 
            top_n(1, wt = `#Reads`))[c('TotalReads', cell_col)] %>% 
            unique
        totals <- data$TotalReads
        o <- order(totals, decreasing = T)
        stuff <- rle(totals[o])
        run.totals <- stuff$values
        run.rank <- cumsum(stuff$lengths) - (stuff$lengths - 
            1)/2
        y <- log10(run.totals)
        x <- (run.rank)
        fit <- smooth.spline(x, y, df = 5)
        d1 <- predict(fit, deriv = 1)$y
        d2 <- predict(fit, deriv = 2)$y
        curvature <- d2/(1 + d1^2)^1.5
       
        plot1 <- ggplot(data = tibble(x = 1:length(curvature), 
            y = curvature), aes(x, y)) + geom_point() + xlab("Rank") + 
            ylab("Calculated Curvature") + theme_gy() + theme(axis.text.x = element_text(hjust = 0.75)) + 
            ggtitle(Plate) + geom_vline(xintercept = which.min(curvature), 
            linetype = "dashed", color = "red", size = 1)
       
        plot2 <- ggplot(data = tibble(x, y), aes(x, y)) + geom_point() + 
            xlab("Rank") + ylab("log10(Reads) \n per cell") + 
            theme_gy() + theme(axis.text.x = element_text(hjust = 0.75)) + 
            ggtitle(Plate) + geom_hline(yintercept = y[which.min(curvature)], 
            linetype = "dashed", color = "red", size = 1)
        
        p = plot1 + plot2 + plot_layout(nrow = 1)
        print(p)
     # print(plot1) | print(plot2)
        return(cutoff = 10^y[which.min(curvature)])
    })
    names(ReadFilters) <- Plates
    
    return(bind_rows(lapply(Plates, function(x) {
        mat[(mat[[plate_col]]==x) & mat[['TotalReads']] >= ReadFilters[x],]
    })))
    }

# A function for filtering alleles per cell. Inspired by EmptryDroplets from Lun et al. Genome Bio 2019
Filtering_Alleles <- function (alleles, plate_col, cell_col, cutoff = 10) 
{
    matrix = alleles
    # matrix <- alleles %>%  mutate(plate_well = paste0(Barcode_DNA, 
    #     Well_ID)) # create a plate_well identifier. 
    cells <- unique(matrix[[plate_col]]) # pull out unique cells
    Plates <- unique(matrix[[plate_col]]) # Pull out unique plates
    
    #Define allele threshold cutoff. Ie number of alleles to keep per cell. 
    Allele_Threshold <- lapply(cells, function(x) {
        data <- matrix[matrix[[plate_col]]==x, ] %>% 
            dplyr::select(`%Reads`) # filter on cell
        totals <- data$`%Reads` # pull out reads
        o <- order(totals, decreasing = T) # order 
        stuff <- rle(totals[o]) # Compute the lengths and values of runs of equal values in a vector
        run.totals <- stuff$value # pull out reads
        run.rank <- cumsum(stuff$lengths) - (stuff$lengths - 
            1)/2 # define rank allowing for ties.
        y <- log10(run.totals) # transform #reads by log10
        x <- (run.rank) # show ranks. 
        return(which.min(diff(y)/diff(x))) # differentiate values and return minimum.
    })
    names(Allele_Threshold) <- cells # add back cell values
    Allele_Threshold <- lapply(Allele_Threshold, function(x) ifelse(length(x) == 0, 1, x)) # Fix empty which.min function outputs - happens when only 1 allele is recovered
    
    #Filter matrix based on the threshold defined. 
    Filter_matrix <- bind_rows(
        lapply(cells, function(x) {
        if ((Allele_Threshold[x] > 0) == T) {     # Account for only 1 allele recovered.
        matrix[matrix[[plate_col]]==x, ] %>% top_n(n = Allele_Threshold[x], 
            wt = `%Reads`)}
            else{
                matrix[matrix[[plate_col]]==x, ] %>% top_n(n = 1, 
            wt = `%Reads`)}
    }))

    #Filter alleles below a minimum threshold. For example, if Allele A represents less than 5% of total alleled recovered in that cell it's removed. This is background
    Filter_matrix <- dplyr::filter(Filter_matrix, `%Reads` >= cutoff)
    return(Filter_matrix)
}
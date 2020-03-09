## differential analysis by DESeq2 ##

# these are commonly used functions for plots and writing results to excel spreadsheet
# author: Zhiyi Liu
# version 01

library(DESeq2)
library(vsn)
library(ggplot2)
library(NMF)
library(ComplexHeatmap)
library(data.table)
library(stringr)

# generate the DESeqDataSet
de <- function(countMatrix, tmntMatrix, condition){
        des <- DESeqDataSetFromMatrix(countData = countMatrix,
                                      colData = tmntMatrix,
                                      design = ~ condition)        
        des <-des[rowSums(counts(des)) > 0,] # remove genes without any counts
        des <- estimateSizeFactors(des) #calculate the size factor and add it to the dataset
}

# get the normalized counts
logNorm <- function(des){
        normCounts <- counts(des, normalized = TRUE)
        log.normCounts <- log2(normCounts + 1) # transform normalized counts to log2 scale
}

# boxplot - log2 transformed counts statistics
box_plot <- function(countData, notch = T, 
                     title = "log2-transformed read counts",
                     ylab = "log2(read counts"){
        boxplot(countData, 
                notch = notch, 
                las = 3,
                main = title,
                ylab = ylab)
}

# mean-sd plot
msd_plot <- function(des){
        des.rlog <- rlog(des, blind = T) # obtain regularized log-transformed values
        rlog.normCounts <- assay (des.rlog)
        msd_plot <- meanSdPlot(rlog.normCounts,
                               ranks = F,
                               plot = F) # mean-sd plot for rlog-transformed data
        msd_plot <- msd_plot$gg +
                ggtitle("rlog-transformed read counts") +
                ylab("Standard deviation")
        print(msd_plot)
}

# dendrogram: to examine the biological replicates
den_plot <- function(des){
        des.rlog <- rlog(des, blind = T) # obtain regularized log-transformed values
        rlog.normCounts <- assay (des.rlog)
        dist.rlog <- as.dist(1 - cor(rlog.normCounts, method = "pearson"))
        plot(hclust(dist.rlog),
             labels = colnames(rlog.normCounts),
             main = "rlog transformed red counts\nDistance: Pearson correlation")
}

# PCA plot
pca_plot <- function(des){
        des.rlog <- rlog(des, blind = T)
        pca <- plotPCA(des.rlog)
        pca <- pca + theme_bw() + ggtitle("rlog transformed counts")
        print(pca)
}

########################## Differential Analysis########################
# DE analysis
de_analysis <- function(des){
        des <- DESeq(des)
        de.results <- results(des, independentFiltering = T, alpha = 0.05)
}


# heatmaps
ht_map <- function(de.results, log.normCounts){
        # a heat map needs a matrix of values, e.g. log2-transformed read counts
        de.results.sorted <- de.results[order(de.results$padj), ] # sort the genes according to adj p vals
        de.genes <- rownames(subset(de.results.sorted, padj < 0.05)) # genes with p < 0.05
        mat.de.genes <- log.normCounts[de.genes, ] # extract the normalized read counts for DE genes into a matrix
        # use: library(ComplexHeatmap)
        hm.mat.de.genes <- t(scale(t(mat.de.genes), center = T, scale = T)) # center the log counts by row, then scale
        ht <- Heatmap(hm.mat.de.genes, 
                      show_row_names = F, 
                      name = "Differentially expressed genes (p < 0.05)",
                      row_dend_reorder = T, 
                      cluster_columns = F,
                      clustering_distance_columns = "euclidean",
                      clustering_distance_rows = "pearson")
        draw(ht)
}

# volcano plots 
vol_plot <- function(de.results){
        de.results$sig <- -log10(de.results$padj) # compute significance by -log10
        # volcano plot of adjusted p-values
        cols <- densCols(de.results$log2FoldChange, de.results$sig)
        cols[de.results$pvalue ==0] <- "purple"
        de.results$pch <- 19 # plotting point shape to be dots
        de.results$pch[de.results$pvalue == 0] <- 6 # if p-val is 0 then plot shape is 6 (reversed triangle)
        plot(de.results$log2FoldChange,
             de.results$sig,
             col = cols, 
             panel.first = grid(),
             main = "Volcano plot",
             xlab = "Effect size:log2(fold-change)",
             ylab = "-log10(adjusted p-value",
             pch = de.results$pch, 
             cex = 0.4)
        abline(v = 0)
        abline(v = c(-1,1), col = "brown")
        alpha <- 0.05 # set the significant level for the line
        abline(h = -log10(alpha), col = "brown")
# plot the names of a reasonable number of genes with |foldchange| > 2 & p < alpha

# gn.selected <- abs(de.results$log2FoldChange) > 2 & de.results$padj < alpha
# text(de.results$log2FoldChange[gn.selected],
#      -log10(de.results$padj)[gn.selected],
#      lab = rownames(de.results)[gn.selected], cex = 0.6)
}

# prepare matrix for GO and KEGG analysis
gkmatrix <- function(de.results){
        de.genes <- rownames(subset(de.results, padj < 0.05))
        
        # get the normalized counts
        normCounts <- counts(des, normalized = TRUE)

        # extract the normalized read counts for DE genes into a matrix
        mat.de.genes <- normCounts[de.genes, ]
        
        diff.gene <- de.results[de.genes,] # get the differential genes with p < 0.05
        diff.gene <- as.data.frame(diff.gene) # convert to matrix
        diff.gene <- diff.gene[,c("log2FoldChange", "pvalue", "padj")]
        
        # convert the rowname into the fist column of the matrix
        setDT(diff.gene, keep.rownames = "id")[] # library(data.table)
        
        # split the MSTRG/ENSG ids from the gene symbols
        gene.ids <- str_split_fixed(diff.gene[,id], "\\|", n = 2) #library(stringr)
        gene.ids <- as.data.table(gene.ids)
        setnames(gene.ids, 1:2, c("id_number","gene_symbol"))
        diff.gene <- cbind(diff.gene, gene.ids)
        setcolorder(diff.gene, c("id", "id_number", "gene_symbol", "log2FoldChange", "pvalue", "padj"))
        
        diff.gene <- diff.gene[!(grepl("MSTRG", id_number) & gene_symbol == ""),] # remove rows only containing MSTRG ids but not gene symbols -- these are predicted genes and not be used for GO/KEGG analysis
}








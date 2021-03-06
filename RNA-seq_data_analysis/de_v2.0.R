setwd("/home/RNA-Seq/data_mnt/ftp/LCS7823/08DEanalysis")
#allow the use of string commands pipe %>%
library(magrittr)

# read in the counts table and the sample information table

readCounts <- read.csv("gene_count_matrix.csv")
sampleInfo <- read.csv("sampleInfo.csv")

# row.names correspond to the gene_id in the dataframe.
row.names(readCounts) <- readCounts$gene_id

# get rid of the first column to just retain the counts data
readCounts <- readCounts[,-1]

# clean up sampleInfo table to be supplied to the DEseq object
# sampleId <- sampleInfo[,1]
# sampleInfo <- data.frame(sampleInfo[,-1])
# row.names(sampleInfo) <- sampleId
# colnames(sampleInfo) <- c("condition")

# take subset of PLC to work on the DE analysis
plc.counts <- readCounts[,7:12]
# generate sampleInfo for PLC cells
plc.sampleInfo <- data.frame(condition = (c(rep("ctrl", 3), rep("drug", 3))), row.names = names(plc.counts))

hepG.counts <- readCounts[,1:6]
hepG.sampleInfo <- data.frame(condition = (c(rep("ctrl", 3), rep("drug", 3))), row.names = names(hepG.counts))

#define de function

de <- function(countMatrix, tmntMatrix){
        #generate the DESeqDataSet
        library(DESeq2)
        
        des <- DESeqDataSetFromMatrix(countData = countMatrix,
              colData = tmntMatrix,
              design = ~ condition)
        
        # remove genes without any counts
        des <-des[rowSums(counts(des)) > 0,]
        
        #calculate the size factor and add it to the dataset
        des <- estimateSizeFactors(des)
        sizeFactors(des)
        
        # get the normalized counts
        normCounts <- counts(des, normalized = TRUE)
        
        # transform normalized counts to log2 scale
        log.normCounts <- log2(normCounts + 1)
        
        # boxplot
        jpeg("Read counts statistics.jpg", width = 500, height = 500)
        
        par(mar = c(10, 4, 4, 2) + 0.1)
        boxplot(log.normCounts, notch = T, las = 3,
                main = "log2-transformed read counts",
                ylab = "log2(read counts")
        dev.off()
        
        # mean-sd plot
        library(vsn)
        library(ggplot2)
        
        jpeg("Mean-Sd_plot.jpg", width = 500, height = 500)
        # obtain regularized log-transformed values
        des.rlog <- rlog(des, blind = T)
        rlog.normCounts <- assay (des.rlog)
        # mean-sd plot for rlog-transformed data
        msd_plot <- meanSdPlot(rlog.normCounts,
                ranks = F,
                plot = F)
        
        msd_plot <- msd_plot$gg +
                        ggtitle("rlog-transformed read counts") +
                        ylab("Standard deviation")
        
        print(msd_plot)
        dev.off()
        
        # dendrogram to examine the biological replicates
        jpeg("Clustering_of_samples.jpg", width = 500, height = 500)
        dist.rlog <- as.dist(1 - cor(rlog.normCounts, method = "pearson"))
        plot(hclust(dist.rlog),
                labels = colnames(rlog.normCounts),
                main = "rlog transformed red counts\nDistance: Pearson correlation")
        dev.off()
        
        jpeg("PCA_plot.jpg", width = 500, height = 500)
        # PCA analysis
        pca <- plotPCA(des.rlog)
        pca <- pca + theme_bw() + ggtitle("rlog transformed counts")
        print(pca)
        dev.off()
        
        ########################## Differential Analysis########################3
        # check the levels of the condition for des object
        # str(colData(des)$condition) # 'ctrl' is the reference level so we are good
        
        # DE analysis
        des <- DESeq(des)
        de.results <- results(des, independentFiltering = T, alpha = 0.05)
        
        # frequences of p-values
        jpeg("Distribution_of_p-values.jpg", width = 500, height = 500)
        hist(de.results$pvalue,
                col = "blue", border = "white", xlab = "p values", ylab = "Number of genes",
                main = "Frequences of p-values")
        dev.off()
        
        # MA plot
        jpeg("MA_plot.jpg", width = 500, height = 500)
        plotMA(de.results, alpha = 0.05, main = "Drug-treated vs control cells", ylim = c(-5, 5))
        dev.off()
        
        
        # heatmaps
        library(NMF)
        
        # a heat map needs a matrix of values, e.g. log2-transformed read counts
        de.results.sorted <- de.results[order(de.results$padj), ] # sort the genes according to adj p vals
        
        # genes with p < 0.05
        de.genes <- rownames(subset(de.results.sorted, padj < 0.05))
        
        # extract the normalized read counts for DE genes into a matrix
        mat.de.genes <- log.normCounts[de.genes, ]
        
        library(ComplexHeatmap)
        
        hm.mat.de.genes <- t(scale(t(mat.de.genes), center = T, scale = T)) # center the log counts by row, then scale
        jpeg("Heatmap_of_differential_genes.jpg", width = 900, height = 800)
        ht <- Heatmap(hm.mat.de.genes, 
                show_row_names = F, 
                name = "Differentially expressed genes (p < 0.05)",
                row_dend_reorder = T, cluster_columns = F,
                clustering_distance_columns = "euclidean",
                clustering_distance_rows = "pearson")
        draw(ht)
        dev.off()
        
        # write the DE results into a spreadsheet
        library(xlsx)
        write.xlsx2(as.data.frame(de.results.sorted),
                    file = "Differential genes orfered by p values.xlsx",
                    sheetName = "Differentail genes ordered by p values",
                    col.names = T, row.names = T)
        
}






































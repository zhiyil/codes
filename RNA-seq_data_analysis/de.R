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
sampleId <- sampleInfo[,1]
sampleInfo <- data.frame(sampleInfo[,-1])
row.names(sampleInfo) <- sampleId
colnames(sampleInfo) <- c("condition")

# take subset of PLC to work on the DE analysis
plc.counts <- readCounts[,c(7,8,9,10,11,12)]
# generate sampleInfo for PLC cells
plc.sampleInfo <- data.frame(condition = (c(rep("ctrl", 3), rep("drug", 3))), row.names = names(plc.counts))

#generate the DESeqDataSet
library(DESeq2)

des <- DESeqDataSetFromMatrix(countData = plc.counts,
                              colData = plc.sampleInfo,
                              design = ~ condition)
#check the object if needed:
## colData(des) %>% head
## assay(des, "counts") %>% head
## rowData(des) %>% head
## counts(des) %>% str

# remove genes without any counts
des <-des[rowSums(counts(des)) > 0,]

# check colSums(counts(des)) and colSums(plc.counts), which should be the same

#calculate the size factor and add it to the dataset
des <- estimateSizeFactors(des)
sizeFactors(des)

# if needed, check colData: colData(des), which now contains the sizeFactors column

# get the normalized counts
normCounts <- counts(des, normalized = TRUE)

# transform normalized counts to log2 scale
log.normCounts <- log2(normCounts + 1)


boxplot(log.normCounts, notch = T, las = 3,
        main = "log2-transformed read counts - PLC cell data",
        ylab = "log2(read counts")

# install vsn
BiocManager::install("vsn")

# mean-sd plot
library(vsn)
library(ggplot2)
msdplot <- meanSdPlot(log.normCounts,
                      ranks = F,
                      plot = T)
msdplot$gg + 
        ggtitle("Sequencing depth normalized log2(read counts)") +
        ylab("Standard deviation")

# obtain regularized log-transformed values
des.rlog <- rlog(des, blind = T)
rlog.normCounts <- assay (des.rlog)

# mean-sd plot for rlog-transformed data
msd_plot <- meanSdPlot(rlog.normCounts,
                       ranks = F,
                       plot = F)
msd_plot$gg +
        ggtitle("rlog-transformed read counts for PLC cells") +
        ylab("Standard deviation")

# dendrogram to examine the biological replicates
dist.rlog <- as.dist(1 - cor(rlog.normCounts, method = "pearson"))
plot(hclust(dist.rlog),
     labels = colnames(rlog.normCounts),
     main = "rlog transformed red counts\nDistance: Pearson correlation")

# PCA analysis
pca <- plotPCA(des.rlog)
pca <- pca + theme_bw() + ggtitle("rlog transformed counts")
print(pca)

########################## Differential Analysis########################3
# check the levels of the condition for des object
str(colData(des)$condition) # 'ctrl' is the reference level so we are good

# DE analysis
des <- DESeq(des)
de.results <- results(des, independentFiltering = T, alpha = 0.05)
# examine de.results: 
# head(de.results)
# table(de.results$padj < 0.05)
# DE gene list can be obtained by:
# rownames(subset(de.results, padj < 0.05))

# frequences of p-values
hist(de.results$pvalue,
     col = "blue", border = "white", xlab = "p values", ylab = "Number of genes",
     main = "Frequences of p-values")

plotMA(de.results, alpha = 0.05, main = "Drug-treated vs control cells", ylim = c(-4, 4))

# heatmaps
library(NMF)

# a heat map needs a matrix of values, e.g. log2-transformed read counts
de.results.sorted <- de.results[order(de.results$padj), ] # sort the genes according to adj p vals

# genes with p < 0.05
de.genes <- rownames(subset(de.results.sorted, padj < 0.05))

# extract the normalized read counts for DE genes into a matrix
mat.de.genes <- log.normCounts[de.genes, ]

# heatmap - normalized read counts for DE gene with adj p < 0.05
# aheatmap(hm.mat.de.genes,
#          Rowv = T, Colv = T,
#          distfun = "euclidean", hclustfun = "average",
#          scale = "row") # values are "centered" before plotting

library(ComplexHeatmap)

hm.mat.de.genes <- t(scale(t(mat.de.genes), center = T, scale = T)) # center the log counts by row, then scale

Heatmap(hm.mat.de.genes, show_row_names = F, name = "Differentially expressed genes (p < 0.05)",
        row_dend_reorder = T, cluster_columns = F,
        clustering_distance_columns = "euclidean",
        clustering_distance_rows = "pearson")


# write the DE results into a spreadsheet
library(xlsx)

write.xlsx2(as.data.frame(de.results.sorted), 
          file = "Differential genes ordered by p values.xlsx",
          sheetName = "Differential genes ordered by p values", 
          col.names= T, row.names = T)


# map gene names
# library(org.Hs.eg.db)
# 
# # types of keywords that are available for query: keytypes(org.Hs.eg.db)
# 
# anno <- select(org.Hs.eg.db,
#                keys = de.genes, keytype = "ENSEMBL",
#                columns = c("SYMBOL","GENENAME", "GO", "ONTOLOGY"))
































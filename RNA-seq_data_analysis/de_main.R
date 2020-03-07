#!/usr/bin/Rscript

#setwd()

source("de_funcs.R")

# read in data files from spreadsheet
readCounts <- read.csv("gene_count_matrix.csv")
row.names(readCounts) <- readCounts$gene_id # row.names correspond to the gene_id in the dataframe.
readCounts <- readCounts[,-1] # get rid of the first column to just retain the counts data


### temp #####
sampleInfo <- read.csv("sampleInfo.csv")
compChar <- colnames(sampleInfo)[1] # set the comparison strings, will merge into file names that generated in the following spreadsheets and plots
sampleInfo <- as.data.frame(sampleInfo[,2])
row.names(sampleInfo) <- colnames(readCounts)
colnames(sampleInfo) <- "condition"
###### end of temp #######

# generate DESeq2DataSet
des <- de(countMatrix = readCounts, tmntMatrix = sampleInfo, condition = condition)

# log2 of normalized data
log.norm <- logNorm(des)

# set margin parameters for all plots here:
plot.margin <- c(10, 4, 4, 2) + 0.1

# function for making different types of plots

# jpgPlot <- function(fileName, compChar = compChar, plotType, width, height, parSetting){
#         if (plotType == "read_counts"){
#                 fileName <- compChar+" Read counts statistics.jpg"
#         }
#         
#         fileName <- compChar
#         jpeg(fileName)
#         
# }

# boxplot
jpeg(paste(compChar," Read counts statistics.jpg"), width = 500, height = 500)
par(mar = plot.margin)
box_plot(log.norm)
dev.off()

# mean-sd plot
jpeg(paste(compChar," Mean-Sd_plot.jpg"), width = 500, height = 500)
par(mar = plot.margin)
msd_plot(des)
dev.off()

# dendrogram
jpeg(paste(compChar," Clustering_of_samples.jpg"), width = 500, height = 500)
den_plot(des)
dev.off()

# PCA plot
jpeg(paste(compChar," PCA_plot.jpg"), width = 500, height = 500)
pca_plot(des)
dev.off()

# DE analysis
de.results <- de_analysis(des)

# p-value plot
jpeg(paste(compChar," Distribution_of_p-values.jpg"), width = 500, height = 500)
hist(de.results$pvalue,
     col = "blue", 
     border = "white", 
     xlab = "p values", 
     ylab = "Number of genes",
     main = "Frequences of p-values")
dev.off()

# MA plot
jpeg(paste(compChar," MA_plot.jpg"), width = 500, height = 500)
plotMA(de.results, alpha = 0.05, main = "MA plot", ylim = c(-5, 5))
dev.off()

# heat map
jpeg(paste(compChar," Heatmap_of_differential_genes.jpg"), width = 900, height = 800)
ht_map(de.results, log.norm)
dev.off()

# volcano plot
jpeg(paste(compChar," Volcano_plot.jpg"), width = 900, height = 900)
par(mar = plot.margin)
vol_plot(de.results)
dev.off()

# write the DE results into a spreadsheet
library(xlsx)
de.results.sorted <- de.results[order(de.results$padj), ]
write.xlsx2(as.data.frame(de.results.sorted),
            file = paste(compChar," Differential genes ordered by p values.xlsx"),
            sheetName = "Differentail genes",
            col.names = T, 
            row.names = T)














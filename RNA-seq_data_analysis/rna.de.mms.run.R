#!/usr/bin/Rscript

dir <- system("pwd", intern = TRUE)
setwd(dir)
source("de_funcs.R")
library(xlsx)

# read in counts data file from spreadsheet
totalCounts <- read.csv("gene_count_matrix.csv")
# readCounts <- read.csv("gene_count_matrix.csv")
row.names(totalCounts) <- totalCounts$gene_id # row.names correspond to the gene_id in the dataframe.
totalCounts <- totalCounts[,-1] # get rid of the first column to just retain the counts data

smplFile <- file.path(dir, "sampleInfo.xlsx")

# library(XLConnect)
# wb <- loadWorkbook(smplFile)
# sheetNames <- getSheets(wb)# get all sheet names
# detach("package:XLConnect", unload=T)
# sheetNames <- c("WT_TNFvsWT_Vehicle", "KO_VehiclevsWT_Vehicle", "KO_TNFvsKO_Vehicle", "KO_TNFvsWT_TNF", "All")

sheetNames <- c("KO_TNFvsKO_Vehicle", "KO_TNFvsWT_TNF", "All")

for (sheet in sheetNames) {
        sampleInfo <- read.xlsx(file = smplFile, sheetName = sheet, header = TRUE) # read in the corresponding sheet for sample information and comparison
        compChar <- colnames(sampleInfo)[1] # set the comparison strings, will merge into file names that generated in the following spreadsheets and plots, also will be used for making sub-directory by the name
        rnames <- as.character(sampleInfo[,1])
        sampleInfo <- as.data.frame(sampleInfo[,2])
        row.names(sampleInfo) <- rnames
        colnames(sampleInfo) <- "condition"
        
        readCounts <- totalCounts[,rnames] # subset relevant data by rnames from the total raw counts table
        workingDir <- file.path(dir, compChar)
        dir.create(workingDir)
        setwd(workingDir)
        
        ############################################################
        ######################## Main work #########################
        # generate DESeq2DataSet
        des <- de(countMatrix = readCounts, tmntMatrix = sampleInfo, condition = condition)
        
        # log2 of normalized data
        log.norm <- logNorm(des)
        
        # set margin parameters for all plots here:
        plot.margin <- c(10, 4, 4, 2) + 0.1
        
        
        # boxplot
        jpeg(paste(compChar,"Read counts statistics.jpg"), width = 500, height = 500)
        par(mar = plot.margin)
        box_plot(log.norm)
        dev.off()
        
        # mean-sd plot
        jpeg(paste(compChar,"Mean-Sd_plot.jpg"), width = 500, height = 500)
        par(mar = plot.margin)
        msd_plot(des)
        dev.off()
        
        # dendrogram
        jpeg(paste(compChar,"Clustering_of_samples.jpg"), width = 500, height = 500)
        den_plot(des)
        dev.off()
        
        # PCA plot
        jpeg(paste(compChar,"PCA_plot.jpg"), width = 500, height = 500)
        pca_plot(des)
        dev.off()
        
        # DE analysis
        de.results <- de_analysis(des)
        
        # p-value plot
        jpeg(paste(compChar,"Distribution_of_p-values.jpg"), width = 500, height = 500)
        hist(de.results$pvalue,
             col = "blue", 
             border = "white", 
             xlab = "p values", 
             ylab = "Number of genes",
             main = "Frequences of p-values")
        dev.off()
        
        # MA plot
        jpeg(paste(compChar,"MA_plot.jpg"), width = 500, height = 500)
        plotMA(de.results, alpha = 0.05, main = "MA plot", ylim = c(-5, 5))
        dev.off()
        
        # heat map
        jpeg(paste(compChar,"Heatmap_of_differential_genes.jpg"), width = 900, height = 800)
        ht_map(de.results, log.norm)
        dev.off()
        
        # volcano plot
        jpeg(paste(compChar,"Volcano_plot.jpg"), width = 900, height = 900)
        par(mar = plot.margin)
        vol_plot(de.results)
        dev.off()
        
        # write the DE results into a spreadsheet
        # library(xlsx)
        de.results.sorted <- de.results[order(de.results$padj), ]
        write.xlsx2(as.data.frame(de.results.sorted),
                    file = paste(compChar,"Differential genes ordered by p values.xlsx"),
                    sheetName = "Differential genes",
                    col.names = T, 
                    row.names = T)
        
        # go analysis - generate the matrix for go analysis
        diff.gene <- gkmatrix(de.results)
        
        ############### Go enrichment analysis ###################
        library(clusterProfiler)
        # library(org.Hs.eg.db)
        library(org.Mm.eg.db)
        
        gsymbols <- diff.gene$gene_symbol
        ego <- enrichGO(gsymbols,
                        OrgDb = org.Mm.eg.db,
                        keyType = "SYMBOL",
                        ont = "ALL",
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.1)
        
        ego_result <- as.data.frame(ego)
        
        ## write the GO results into Excel
        # library(xlsx)
        write.xlsx2(ego_result,
                    file = paste(compChar, "GO enrichment results on differential genes.xlsx"),
                    col.names = T, row.names = F)
        
        
        ### bubble plot ####
        library(enrichplot)
        jpeg(paste(compChar, "Bubble_plot_GO_terms.jpg"), width = 1000, height = 1300)
        par(mar = plot.margin)
        print(dotplot(ego, showCategory = 20))
        dev.off()
        
        ### network plotting ####
        plotGenes <- setNames(as.character(diff.gene$gene_symbol), 2^diff.gene$log2FoldChange) # get a vector containing fold chage values with names being the gene symbols
        
        jpeg(paste(compChar, "Network_plot_GO_terms.jpg"), width = 1000, height = 1000)
        par(mar = plot.margin)
        print(cnetplot(ego, categorySize = "pvalue", foldChange = plotGenes))
        dev.off()
        
        jpeg(paste(compChar, "Network_circular_plot_GO_terms.jpg"), width = 1400, height = 1200)
        par(mar = plot.margin)
        print(cnetplot(ego, foldChange = plotGenes, circular = T, colorEdge = T))
        dev.off()
        
        ########################################################################
        ## KEGG analysis ##
        
        # gene.kegg <- bitr(diff.gene$gene_symbol,
        #                   fromType = "SYMBOL",
        #                   toType = "ENTREZID",
        #                   OrgDb = org.Hs.eg.db) # use bitr to map gene symbols to entrezids
        gene.kegg <- mapIds(org.Mm.eg.db, diff.gene$gene_symbol, "ENTREZID", "SYMBOL") # here use mapIds from org.Hs.eg.db to do the mapping, !!warning!! 1:multiple mapping obtained
        
        ekg <- enrichKEGG(gene.kegg,
                          organism = "mmu",
                          keyType = "kegg",
                          pvalueCutoff = 0.05,
                          pAdjustMethod = "BH",
                          qvalueCutoff = 0.10)
        
        ekg_result <- as.data.frame(ekg)
        
        ## write the KEGG results into Excel
        # library(xlsx)
        write.xlsx2(ekg_result,
                    file = paste(compChar, "KEGG enrichment results on differential genes.xlsx"),
                    col.names = T, row.names = F)
        
        
        ### bubble plot ####
        library(enrichplot)
        jpeg(paste(compChar, "Bubble_plot_KEGG_terms.jpg"), width = 1000, height = 1300)
        par(mar = plot.margin)
        print(dotplot(ekg, showCategory = 20))
        dev.off()
        
        ### network plotting ####
        plotGenes <- setNames(as.character(diff.gene$gene_symbol), 2^diff.gene$log2FoldChange) # get a vector containing fold chage values with names being the gene symbols
        ekg.symbol <- setReadable(ekg, "org.Mm.eg.db", "ENTREZID") #convert entrezid to gene symbol, making the following cnetplot readable
        
        jpeg(paste(compChar, "Network_plot_KEGG_terms.jpg"), width = 1000, height = 1000)
        par(mar = plot.margin)
        print(cnetplot(ekg.symbol, categorySize = "pvalue", foldChange = plotGenes))
        dev.off()
        
        jpeg(paste(compChar, "Network_circular_plot_KEGG_terms.jpg"), width = 1400, height = 1200)
        par(mar = plot.margin)
        print(cnetplot(ekg.symbol, foldChange = plotGenes, circular = T, colorEdge = T))
        dev.off()
        
        writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
        }
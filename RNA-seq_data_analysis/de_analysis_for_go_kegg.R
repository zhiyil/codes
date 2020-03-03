setwd("/home/RNA-Seq/data_mnt/ftp/LCS7823_analysis/08DEanalysis")
#allow the use of string commands pipe %>%
library(magrittr)

# read in the counts table and the sample information table

readCounts <- read.csv("gene_count_matrix.csv")


# row.names correspond to the gene_id in the dataframe.
row.names(readCounts) <- readCounts$gene_id

# get rid of the first column to just retain the counts data
readCounts <- readCounts[,-1]

# clean up sampleInfo table to be supplied to the DEseq object


# take subset of PLC to work on the DE analysis
plc.counts <- readCounts[,7:12]
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

########################## Differential Analysis########################3
# check the levels of the condition for des object
# str(colData(des)$condition) # 'ctrl' is the reference level so we are good

# DE analysis
des <- DESeq(des)
de.results <- results(des, independentFiltering = T, alpha = 0.05)
# examine de.results: 
# head(de.results)
# table(de.results$padj < 0.05)
# DE gene list can be obtained by:
de.genes <- rownames(subset(de.results, padj < 0.05))

# extract the normalized read counts for DE genes into a matrix
mat.de.genes <- normCounts[de.genes, ]

diff.gene <- de.results[de.genes,] # get the differential genes with p < 0.05
diff.gene <- as.data.frame(diff.gene) # convert to matrix
diff.gene <- diff.gene[,c("log2FoldChange", "pvalue", "padj")]

# convert the rowname into the fist column of the matrix
library(data.table)
setDT(diff.gene, keep.rownames = "id")[]


# split the MSTRG/ENSG ids from the gene symbols
library(stringr)
gene.ids <- str_split_fixed(diff.gene[,id], "\\|", n = 2)
gene.ids <- as.data.table(gene.ids)
setnames(gene.ids, 1:2, c("id_number","gene_symbol"))
diff.gene <- cbind(diff.gene, gene.ids)
setcolorder(diff.gene, c("id", "id_number", "gene_symbol", "log2FoldChange", "pvalue", "padj"))

diff.gene <- diff.gene[!(grepl("MSTRG", id_number) & gene_symbol == ""),] # remove rows only containing MSTRG ids but not gene symbols -- these are predicted genes and not be used for GO/KEGG analysis


############### Go enrichment analysis ###################
library(clusterProfiler)
library(org.Hs.eg.db)

gsymbols <- diff.gene$gene_symbol
ego <- enrichGO(gsymbols,
                OrgDb = org.Hs.eg.db,
                keyType = "SYMBOL",
                ont = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.1)

ego_result <- as.data.frame(ego)

## write the GO results into Excel
library(xlsx)
write.xlsx2(ego_result,
            file = "GO enrichment results on differential genes.xlsx",
            col.names = T, row.names = F)


### bubble plot ####
library(enrichplot)
jpeg("Bubble_plot_GO_terms.jpg", width = 1000, height = 1300)
par(mar = c(10, 4, 4, 2) + 0.1)

dotplot(ego, showCategory = 20)

dev.off()

### network plotting ####
plotGenes <- setNames(as.character(diff.gene$gene_symbol), 2^diff.gene$log2FoldChange) # get a vector containing fold chage values with names being the gene symbols

jpeg("Network_plot_GO_terms.jpg", width = 1000, height = 1000)
par(mar = c(10, 4, 4, 2) + 0.1)
cnetplot(ego, categorySize = "pvalue", foldChange = plotGenes)
dev.off()

jpeg("Network_circular_plot_GO_terms.jpg", width = 1100, height = 800)
par(mar = c(10, 4, 4, 2) + 0.1)
cnetplot(ego, foldChange = plotGenes, circular = T, colorEdge = T)
dev.off()
































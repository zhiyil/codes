# universe gene list from the whole data table 

allGenes <- row.names(normCounts)

allGenes <- gsub("(MSTRG.*|ENSG.*)\\|", "", allGenes) #extract the gene names right after "|" and get rid of the MSTRG# or ENSG# parts.

# selected gene list for GO analysis from the differential analysis


# map gene names
# library(org.Hs.eg.db)
# 
# # types of keywords that are available for query: keytypes(org.Hs.eg.db)
# 
# anno <- select(org.Hs.eg.db,
#                keys = de.genes, keytype = "ENSEMBL",
#                columns = c("SYMBOL","GENENAME", "GO", "ONTOLOGY"))
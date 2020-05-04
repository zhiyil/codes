
library(ComplexHeatmap)
library(xlsx)

workingFiles <- list.files()
workingFiles <- workingFiles[-7]

for (f in workingFiles) {
        df <- read.xlsx(f, 1, header = T)
        df <- df[!rowSums(is.na(df)) == ncol(df),] # remove rows with all NAs
        df <- df[ ,!colSums(is.na(df)) == nrow(df)] # remove columns with all NAs
        df <- unique(df) # remove duplicated rows
        row.names(df) <- df [,1]
        df <- df[,-1]
        hm.mat <- as.matrix(t(scale(t(log2(df+1)), center = T, scale = T))) # center the log counts by row, then scale
        hm.mat <- hm.mat[complete.cases(hm.mat), , drop=F] # remove NAs after transformation, and use 'drop=F' to keep the data structure (i.e. the matrix), otherwise it will simplify to vector
        
        if (nrow(hm.mat) > 30) {genelabel = F} else {genelabel = T} # if gene # > 30 then will not show the gene ids in the heatmap
        basename_f <- sub(".xlsx","",f)
        jpeg(paste0(basename_f, ".heatmap.jpg"), width = 900, height = 800)
        
        ht <- Heatmap(hm.mat, 
              show_row_names = genelabel, 
              name = paste("Heat map of", basename_f),
              row_dend_reorder = T, 
              cluster_columns = F,
              clustering_distance_columns = "euclidean",
              clustering_distance_rows = "pearson")
        draw(ht)
        
        dev.off()
}
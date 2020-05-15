setwd("/Users/zhiyi/Documents/Magic\ Briefcase/LC\ Sciences/customer/7921/plots")

library(xlsx)
library(dplyr)
library(ggplot2)
library(stringr)
#library(viridis)

comparisons <- c("B_TNF_acut vs Ctrl", "B_TNF_chron vs Ctrl", "W_TNF_acut vs Ctrl", "W_TNF_chron vs Ctrl")

for (comparison in comparisons){
        df <- read.xlsx(file = "forPlot.xlsx", sheetName = comparison, header = T)
        colnames(df) <- c("Gene_Set_Name", "Num_genes_in_set", "Description", "GeneCount", "GeneRatio", "p_value", "FDR")
        df <- df %>% select(Description, GeneCount, GeneRatio, FDR)
        
        if (nrow(df)>15) {
                n = 15
                df <- df[1:n, ]
        } else {
                n=nrow(df)
        }
        
        p <- ggplot(df, aes(y=str_wrap(Description, 45), x=GeneRatio, size=GeneCount, color=FDR)) +
                geom_point(alpha=0.5) +
                scale_color_gradient(low="red", high="blue") +
                xlab("Gene Ratio") +
                ylab("Gene Set")
        ggsave(filename = paste0(comparison,".png"), plot = p, device = "png")
        ggsave(filename = paste0(comparison,".jpg"), plot = p, device = "jpg")
}




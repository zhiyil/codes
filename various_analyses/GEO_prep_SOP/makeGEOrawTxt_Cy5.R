#!/usr/bin/Rscript

wkdir <- system("pwd", intern = T)
setwd(wkdir)

filelist <- system("cat fileList.txt", intern = T)

library(xlsx)

for (f in filelist) {
        datafile <- read.xlsx(file = f, sheetName = "Raw Data", header = T)
        fbase <- colnames(datafile)[9] # extract sample ID information
        fbase <- sub("B.", "", fbase)
        fbase <- sub(".cy5", "", fbase) # get rid of A. and .cy3 to get real sample ID
        datafile <- datafile[,-c(7,8)] # get rid of empty columns 9 and 10
        write.table(datafile,file = paste0(fbase,".txt"), sep = "\t", row.names = F, col.names = F)
}


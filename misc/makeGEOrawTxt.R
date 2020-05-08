#!/usr/bin/Rscript

wkdir <- system("pwd", intern = T)
setwd(wkdir)

filelist <- read.table("fileList.txt")
filelist <- as.character(filelist[ ,1]) # convert to vector

library(xlsx)

for (f in filelist) {
        datafile <- read.xlsx(file = f, sheetName = "Raw Data", header = T)
        fbase <- colnames(datafile)[7] # extract sample ID information
        fbase <- sub("A.", "", fbase)
        fbase <- sub(".cy3", "", fbase) # get rid of A. and .cy3 to get real sample ID
        datafile <- datafile[,-c(9,10)] # get rid of empty columns 9 and 10
        write.table(datafile,file = paste0(fbase,".txt"), sep = "\t", row.names = F, col.names = F)
}


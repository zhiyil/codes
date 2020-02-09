library(xlsx)

setwd("/home/RNA-Seq/codes/misc")

df <- read.xlsx("raw.xlsx", sheetIndex = 1, header = T)

library(dplyr)
library(data.table)

dt <- data.table(df)
dt <- dt[Group != "control",] # take out rows containing "control" probes

setorder(dt, Reporter) # reorder the probes so all replicates are put together

dt1 <- dt[1:3232,]
dt2 <- dt[3233:nrow(dt),]

dt1.a <- dt[seq(4, nrow(dt1), 4),] # extract every 4th row, which contains the *_A probes

dt1 <- rbind(dt1, dt1.a, dt1.a) # duplicate the *_A rows and join back to original dt1
setorder(dt1, Reporter)

dt <- rbind(dt1, dt2) # rejoin the first part and second part of the table, dt1 and dt2

setorder(dt, Reporter) # reorder the table, now the table has a long form with each probe containing 3 replicates

dt.rep <- data.table(replicates = c(rep(1:3, nrow(dt)/3)))

dt <- cbind(dt, dt.rep) # add replicates column to dt

dt <- dt[,-4] # delete the last column "Group" - no use for following steps

setcolorder(dt, c("Reporter", "replicates", "BindingSig", "KinaseSig")) # check if the 3 replicates are in place for each probe

clean.dt <- dcast(dt, Reporter ~ replicates, value.var=c("BindingSig", "KinaseSig")) # convert to wide form of the data

#####################################################################################
# calculate the correlation and the p values

m <- as.data.frame(clean.dt) # use data.frame for calculating correlation, as data.table doesn't support row.names

row.names(m) <- m[,1]

m <- m[,-1] # get rid of the "Reporter" column and keep binding data only

library(Hmisc)

cor.res <-rcorr(m)

format.pval(cor.res$P, nsamll=3, digits = 2) # show digits of the P values

library(ComplexHeatmap)

mat.hp <- as.matrix(m)

### should scale the data based on Binding and Kinase subsets, then join the two matrix

# hm.mat.hp <- t(scale(t(mat.hp), center = T, scale = T))
# 
# 
# Heatmap(hm.mat.hp, show_row_names = F, row_dend_reorder = T, cluster_columns = F,
#         clustering_distance_columns = "euclidean",
#         clustering_distance_rows = "pearson")






library(xlsx)

df <- read.xlsx("raw.xlsx", sheetIndex = 1, header = T)

library(dplyr)
library(data.table)

dt <- data.table(df)
dt <- dt[Group != "control",] # take out rows containing "control" probes
dt <- dt[,-4] # delete the last column "Group" - no use for following steps
setorder(dt, Reporter) # reorder the probes so all replicates are put together

dt.a <- dt[seq(4, nrow(dt), 4),] # extract every 4th row, which contains the *_A probes

dt <- rbind(dt, dt.a, dt.a) # duplicate the *_A rows and join back to original dt

setorder(dt, Reporter) # reorder the table, now the table has a long form with each probe containing 3 replicates

dt.rep <- data.table(replicates = c(rep(1:3, nrow(dt)/3)))

dt <- cbind(dt, dt.rep) # add replicates column to dt





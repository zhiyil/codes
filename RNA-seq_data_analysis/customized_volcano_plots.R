de.results$sig <- -log10(de.results$padj) # compute significance by -log10
# volcano plot of adjusted p-values
cols <- densCols(de.results$log2FoldChange, de.results$sig)
cols[de.results$pvalue ==0] <- "purple"
de.results$pch <- 19 # plotting point shape to be dots
de.results$pch[de.results$pvalue == 0] <- 6 # if p-val is 0 then plot shape is 6 (reversed triangle)

# plot(de.results$log2FoldChange,
#      de.results$sig,
#      col = cols, 
#      panel.first = grid(),
#      main = "Volcano plot",
#      xlab = "Effect size:log2(fold-change)",
#      ylab = "-log10(adjusted p-value",
#      pch = de.results$pch,
#      xlim = c(-23, 23),
#      ylim = c(0, 45),
#      cex = 0.4)
# abline(v = 0)
# abline(v = c(-1,1), col = "brown")
# alpha <- 0.05 # set the significant level for the line
# abline(h = -log10(alpha), col = "brown")

voldata <- data.frame(miRNA=row.names(de.results),fold.change=de.results$log2FoldChange, log.adjpval=de.results$sig)
voldata <- voldata[complete.cases(voldata),]

# voldata[voldata$miRNA=="hsa-miR-127-3p",] # where is miR-127-3p?


ggplot(voldata, aes(fold.change, log.adjpval, color = ifelse(abs(fold.change) > 1 & log.adjpval > -log10(0.05), "red", "grey"))) +
        geom_point(alpha = 1/3) +
        xlab("Log2(fold change)") +
        ylab("-Log10(adjusted p value)") +
        xlim(-23, 23) +
        ylim(0, 45) +
        geom_hline(yintercept = -log10(0.05) , color = "brown", linetype = "dotted", size = 0.5) +
        geom_vline(xintercept = c(-1,1), color = "brown", linetype = "dotted", size = 0.5) +
        theme(text = element_text(size = 16), legend.position = "none") +
        scale_color_manual(values = c("grey", "red")) +
        annotate(geom = "curve", x = 15, y = 12, xend = 10.6, yend = 9.8, curvature = 0, arrow = arrow(length = unit(2, "mm"))) +
        annotate(geom = "text", x =15.2, y= 12, label = "miR-127-3p", hjust = "left")
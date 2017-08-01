# from: https://gist.github.com/stephenturner/806e31fce55a8b7175af
setwd("/home/mate/workspace/katamari/src/root/ed/bob/processed")
# res <- read.table("24h_bobdata_ed2_volcano.csv", header=TRUE)
res <- read.csv("D8-rescue2.csv", header=TRUE)
head(res)

# Make a basic volcano plot
# with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Protein expresion after 24h T4 stimulation in WT vs PTPN22 KO cells", xlab="Log2 fold change" ,ylab="log10 p-value (raw)", xlim=c(-3,2.5)))
with(res, plot(log2ratio, -log10(p.value), pch=20, main="D8 cells", xlab="Log2 fold change" ,ylab="log10 p-value (raw)"))


# Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
with(subset(res, p.value<.01 ), points(log2ratio, -log10(p.value), pch=20, col="red"))
with(subset(res, abs(log2ratio)>1), points(log2ratio, -log10(p.value), pch=20, col="orange"))
with(subset(res, p.value<.01 & abs(log2ratio)>1), points(log2ratio, -log10(p.value), pch=20, col="green"))
# install.packages("calibrate")
# Label points with the textxy function from the calibrate plot
library(calibrate)
# with(subset(res, pvalue<.05 & abs(log2FoldChange)>1), textxy(log2FoldChange, -log10(pvalue), labs=Gene, cex=.8))
with(subset(res, p.value<.01 & abs(log2ratio)>1), textxy(log2ratio, -log10(p.value), labs=Gene.names, cex=.8))
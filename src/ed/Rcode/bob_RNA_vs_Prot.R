# plot RNA seq fold changes against mass spec fold changes from Bob's 24H T4 stimulation datasets
setwd("/home/mate/workspace/katamari/src/root/ed/bob/processed")
fileFrame <- read.csv("24h_bobgenes_in_common.csv", header = TRUE)

# plot(fileFrame$l2FC, fileFrame$Average.fold.change.RNA.seq)
# with(fileFrame, plot(l2FC,Average.fold.change.RNA.seq, pch=20, xlim = c(-4.5,4.5), ylim = c(-9,7), main = "effect of Ptpn22 on RNA and protein levels in OT1 cells after 24H T4 stimulation", xlab = "Fold change in protein level (log2)", ylab = "Fold change in RNA level (log2)"))

with(fileFrame, plot(l2FC,Average.log.2.fold.change.RNA.seq, pch=20, main = "effect of Ptpn22 on RNA and protein levels in OT1 cells after 24H T4 stimulation", xlab = "Fold change in protein level (log2)", ylab = "Fold change in RNA level (log2)"))

# add pretty colours
with(subset(fileFrame, abs(l2FC)>2), points(l2FC, Average.log.2.fold.change.RNA.seq, pch=20, col="red"))
with(subset(fileFrame, abs(Average.log.2.fold.change.RNA.seq)>2), points(l2FC, Average.log.2.fold.change.RNA.seq, pch=20, col="orange"))
with(subset(fileFrame, abs(Average.log.2.fold.change.RNA.seq)>2 & abs(l2FC)>2), points(l2FC, Average.log.2.fold.change.RNA.seq, pch=20, col="purple"))

# Label points with the textxy function from the calibrate plot
library(calibrate)
with(subset(fileFrame, abs(l2FC)>8 | abs(Average.log.2.fold.change.RNA.seq)>3), textxy(l2FC, Average.log.2.fold.change.RNA.seq, labs=GeneName, cex=.8))
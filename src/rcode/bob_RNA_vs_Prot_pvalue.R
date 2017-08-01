# plot RNA seq raw p values against mass spec raw p values from Bob's 24H T4 stimulation datasets
setwd("/home/mate/workspace/katamari/src/root/ed/bob/processed")
fileFrame <- read.csv("24h_bobgenes_in_common.csv", header = TRUE)

# plot(fileFrame$l2FC, fileFrame$Average.fold.change.RNA.seq)
# with(fileFrame, plot(l2FC,Average.fold.change.RNA.seq, pch=20, xlim = c(-4.5,4.5), ylim = c(-9,7), main = "effect of Ptpn22 on RNA and protein levels in OT1 cells after 24H T4 stimulation", xlab = "Fold change in protein level (log2)", ylab = "Fold change in RNA level (log2)"))

with(fileFrame, plot(-log10(pValue),-log10(raw.P.value.in.RNA.seq), pch=20, xlim = c(0,4.2), main = "effect of Ptpn22 on RNA and protein levels in OT1 cells after 24H T4 stimulation", xlab = "raw p value of protein level changes (-log10)", ylab = "raw p value of RNA level changes (-log10)"))

# add pretty colours
with(subset(fileFrame, -log10(pValue)>1 | -log10(raw.P.value.in.RNA.seq)>1), points(-log10(pValue),-log10(raw.P.value.in.RNA.seq), pch=20, col="blue"))
with(subset(fileFrame, -log10(pValue)>2 | -log10(raw.P.value.in.RNA.seq)>2), points(-log10(pValue),-log10(raw.P.value.in.RNA.seq), pch=20, col="cadetblue"))
with(subset(fileFrame, -log10(pValue)>3 | -log10(raw.P.value.in.RNA.seq)>3), points(-log10(pValue),-log10(raw.P.value.in.RNA.seq), pch=20, col="chartreuse3"))
with(subset(fileFrame, -log10(pValue)>4 | -log10(raw.P.value.in.RNA.seq)>4), points(-log10(pValue),-log10(raw.P.value.in.RNA.seq), pch=20, col="orange"))


# Label points with the textxy function from the calibrate plot
library(calibrate)
with(subset(fileFrame, sqrt((-log10(pValue)*-log10(pValue)*2.6)+(-log10(raw.P.value.in.RNA.seq)*-log10(raw.P.value.in.RNA.seq))) > 4.3), textxy(-log10(pValue), -log10(raw.P.value.in.RNA.seq), labs=GeneName, cex=.8))
# with(subset(fileFrame, -log10(pValue)>3 | -log10(raw.P.value.in.RNA.seq)>4.5), textxy(-log10(pValue), -log10(raw.P.value.in.RNA.seq), labs=GeneName, cex=.8))

with(subset(fileFrame, GeneName == "Ptpn22"), textxy(-log10(pValue), -log10(raw.P.value.in.RNA.seq), labs=GeneName, cex=.8, col = "red"))
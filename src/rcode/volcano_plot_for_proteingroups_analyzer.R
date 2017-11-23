# plot LFQ foldchanges (wt/ko) vs P value from exosome data
setwd("/home/mate/code/ed/src/data/cav1ko/processed/")
fileFrame <- read.csv("combined.csv", header = TRUE)

# plot(fileFrame$l2FC, fileFrame$Average.fold.change.RNA.seq)
# with(fileFrame, plot(l2FC,Average.fold.change.RNA.seq, pch=20, xlim = c(-4.5,4.5), ylim = c(-9,7), main = "effect of Ptpn22 on RNA and protein levels in OT1 cells after 24H T4 stimulation", xlab = "Fold change in protein level (log2)", ylab = "Fold change in RNA level (log2)"))

with(fileFrame, plot(Log2.Fold.change,P.value, pch=20, xlim = c(-10,10), ylim = c(1.0,0.0), main = "exosome cav1KO vs WT", xlab = "wt/ko Fold change in LFQ (log2)", ylab = "P value"))

# add pretty colours
with(subset(fileFrame, Log2.Fold.change>5), points(Log2.Fold.change, P.value, pch=20, col="orange"))
with(subset(fileFrame, Log2.Fold.change<(-5)), points(Log2.Fold.change, P.value, pch=20, col="orange"))
with(subset(fileFrame, P.value<(0.05)), points(Log2.Fold.change, P.value, pch=20, col="blue"))
with(subset(fileFrame, P.value<(0.05) & Log2.Fold.change<(-5)), points(Log2.Fold.change, P.value, pch=20, col="red"))
with(subset(fileFrame, P.value<(0.05) & Log2.Fold.change>5), points(Log2.Fold.change, P.value, pch=20, col="red"))

# Label points with the textxy function from the calibrate plot

# labelling not done yet, this code needs to be modified before use
library(calibrate)
with(subset(fileFrame, FC.LFQ>(FC.absquant*2)+ 10 & FC.LFQ>10), textxy(FC.LFQ, FC.absquant, labs=Gene.names, cex=.8))
with(subset(fileFrame, FC.absquant> (-6) & FC.LFQ<(-10)), textxy(FC.LFQ, FC.absquant, labs=Gene.names, cex=.8))
with(subset(fileFrame, 3 < FC.absquant & 1 > FC.LFQ), textxy(FC.LFQ, FC.absquant, labs=Gene.names, cex=.8))

with(subset(fileFrame, LFQ.intensity.TCR_KO_1 > 6e10), textxy(Intensity.TCR_KO_1, LFQ.intensity.TCR_KO_1, labs=Gene.names, cex=.8))

with(subset(fileFrame, protein.count.KO_1 > 3e7 | Intensity.TCR_KO_1 > 2e10), textxy(Intensity.TCR_KO_1, protein.count.KO_1, labs=Gene.names, cex=.8))

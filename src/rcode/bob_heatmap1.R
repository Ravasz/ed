# generate heatmap from Bob's data (LFQ) to see how it looks
library(gplots)
library(made4)
setwd("/home/mate/workspace/katamari/src/root/ed/bob")
fileDf <- read.csv("bobdata2_ed2.csv", sep = " ", stringsAsFactors = FALSE, header = TRUE)
heatMT <- as.matrix(fileDf[, 4:9])
heatMT[heatMT == 0] <- sample(10, replace = TRUE)
# heatmap(heatMT, distfun = function(x) dist(x, method = "maximum"), labCol = c("KO1", "KO2", "KO3", "WT1", "WT2", "WT3"))
heatplot (heatMT, labRow = " ", labCol = c("KO1", "KO2", "KO3", "WT1", "WT2", "WT3"))

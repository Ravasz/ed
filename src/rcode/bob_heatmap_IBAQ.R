# generate heatmap from Bob's IBAQ data to see how it looks
library(gplots)
library(made4)
library(scales)
setwd("/home/mate/workspace/katamari/src/root/ed/bob")
fileDf <- read.csv("bob_ibaq.csv", sep = ",", stringsAsFactors = FALSE, header = TRUE)

logD <- as.matrix(fileDf[, 1:6])
logD <- log2(logD)
logD[logD == -Inf] <- sample(100, replace = TRUE) / 1000
# newScale <- rescale(fileDf, to = c(0, 2), from = range(fileDf, na.rm = TRUE))
# newScaleM <- as.matrix(newScale)
# heatmap(logD, distfun = function(x) dist(x, method = "euclidean"), labCol = c("KO1", "KO2", "KO3", "WT1", "WT2", "WT3"))
heatplot (logD, labRow = " ", labCol = c("KO1", "KO2", "KO3", "WT1", "WT2", "WT3"))

# create histogram of p value score distributions
funT <- function(x) tryCatch(t.test(x[1:3], x[4:6])$p.value, error = function(y) 1)
ttestV <- apply(fileDf, 1, funT)
ttestV[is.nan(ttestV)] <- 1
hist(ttestV, breaks = 200, main = "p value distribution using IBAQ intensities", xlab = "p value", col = "chartreuse3")

matrixM <- data.matrix(fileDf)
funW <- function(x) wilcox.test(x[1:3], x[4:6], exact = FALSE)$p-value
wilcox <- apply(fileDf, 1, wilcox.test(1:3, 4:6, exact = FALSE)$p-value)
ttestV[is.nan(ttestV)] <- 1
hist(ttestV, breaks = 200, main = "p value distribution using IBAQ intensities", xlab = "p value", col = "chartreuse3")



# create histogram of IBAQ score distributions
logCol <- log2(fileDf[ , 1])
logCol[logCol == -Inf] <- 0
hist(logCol, breaks = 200, main = "IBAQ scores in the KO1 sample", xlab = "log2(IBAQ)", col = "chartreuse3")
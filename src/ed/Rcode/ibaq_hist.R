# take the IBAQ values from Bob's dataset, put all measurements into a single dataframe, and plot a histogram with their log2 values

setwd("/home/mate/workspace/katamari/src/root/ed/bob")
fileDf <- read.csv("bob_lfq.csv", sep = ",", stringsAsFactors = FALSE, header = TRUE)

logD <- as.matrix(fileDf[, 1:6]) # the log calculation really works if in a matrix
logD <- log2(logD) # convert to log scale
logD[logD == -Inf] <- sample(100, replace = TRUE) / 1000 # replace negative infinite (0) values with very small numbers for statistical pruposes
dfLogD <- as.data.frame(logD) # convert back to data frame

df2 <- do.call(rbind, lapply(1:nrow(dfLogD), function(x) t(dfLogD[x,])))
rownames(df2) <- NULL
head(df2)

hist(df2, breaks = 200, main = "LFQ scores in all samples", xlab = "log2(LFQ)", col = "chartreuse3")

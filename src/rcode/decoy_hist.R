setwd("/home/mate/workspace/katamari/src/root/ed/bob")
fileDf <- read.csv("bob_decoy_lfq.csv", header = TRUE)
funT <- function(x) tryCatch(t.test(x[1:3], x[4:6])$p.value, error = function(y) 1)
fileDf$ttest <- apply(fileDf[,c(1:6)], 1, funT)
# converts table to strings: fileDf = format(fileDf,scientific=FALSE) 
fileDf$X0 <- NULL
write.table(fileDf, file = "bobdata2_decoy_lfq_2.csv")
fileDf$ttest[is.nan(fileDf$ttest)] <- 1
hist(fileDf$ttest, breaks = 200, main = "p value distribution using LFQ intensities", xlab = "p value", col = "chartreuse3")

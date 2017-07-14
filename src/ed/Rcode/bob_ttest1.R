setwd("/home/mate/workspace/katamari/src/root/ed/bob/processed")
fileDf <- read.csv("24h_bobdata_no0_ed.csv", header = TRUE)
funT <- function(x) tryCatch(t.test(x[1:3], x[4:6])$p.value, error = function(y) 1)
fileDf$ttest <- apply(fileDf[,c(5:10)], 1, funT)
fileDf = format(fileDf,scientific=FALSE) 
fileDf$X0 <- NULL
write.table(fileDf, file = "24h_bobdata_ed2.csv")
fileDf$ttest[is.nan(fileDf$ttest)] <- 1
# hist(fileDf$ttest, breaks = 200, main = "p value distribution using LFQ intensities", xlab = "p value", col = "chartreuse3")

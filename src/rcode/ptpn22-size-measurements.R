# create histogram of OST-PTPN22 protein size quantifications

# what an embarrassing bit of code this turned out to be. Should use some kind of looping here

setwd("/home/mate/code/ed/src/data/size_measurements")
fileDfone <- read.csv("ost-ptpn22-sizes1.csv", header = FALSE)
fileDftwo <- read.csv("ost-ptpn22-sizes2.csv", header = FALSE)
fileDftri <- read.csv("ost-ptpn22-sizes3.csv", header = FALSE)
fileDffour <- read.csv("ost-ptpn22-sizes4.csv", header = FALSE)
fileDffive <- read.csv("ost-ptpn22-sizes5.csv", header = FALSE)
fileDfsix <- read.csv("ost-ptpn22-sizes6.csv", header = FALSE)
fileDfsev <- read.csv("ost-ptpn22-sizes7.csv", header = FALSE)
fileDfeig <- read.csv("ost-ptpn22-sizes8.csv", header = FALSE)


hist(fileDfone$V1, breaks = 10, main = "OST-PTPN22 protein sizes", xlab = "protein size [KDa]", xlim = c(50,120), ylim = c(0,25), col = rgb(0,1,0,1))

hist(fileDftwo$V1, breaks = 10, col = rgb(0,0,1,1), add = T)
hist(fileDftri$V1, breaks = 20, col = rgb(1,0,0,1), add = T)
hist(fileDffour$V1, breaks = 20, col = rgb(0.5,0.5,0.5,1), add = T)
hist(fileDffive$V1, breaks = 10, col = rgb(0,0.5,0.5,1), add = T)
hist(fileDfsix$V1, breaks = 5, col = rgb(0.5,0,0.5,1), add = T)
hist(fileDfsev$V1, breaks = 2, col = rgb(0.5,0.6,0.1,1), add = T)
# hist(fileDfeig$V1, breaks = 1, col = rgb(0.7,0.3,0,1), add = T)
legend("topright", c("longest isoform", "second longest", "3rd" , "4th", "5th", "6th", "7th"), fill=c(rgb(0,1,0,1),rgb(0,0,1,1), rgb(1,0,0,1),rgb(0.5,0.5,0.5,1), rgb(0,0.5,0.5,1), rgb(0.5,0,0.5,1), rgb(0.5,0.6,0.1,1)))
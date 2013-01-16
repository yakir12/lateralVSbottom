args <- commandArgs(trailingOnly = TRUE)
library(MASS)
library(png)
data <-read.table("/home/yg32/Documents/PostDoc/Holodeck/Code/table.txt",header=TRUE,sep="\t")
path <- args[1]
# path <- "/home/yg32/Documents/PostDoc/Holodeck/Code/processedata/1_Checkerboard_2_2012_11_16_13_3"
files <- list.files(path,"*.png")
m <- length(files)
I <-readPNG(paste(path,files[1],sep="/"))
sz <- dim(I)
L <- array(0, dim=c(sz,m))
x <- rep(0,m)
for(i in 1:m){
  L[,,i] <- readPNG(paste(path,files[i],sep="/"))
  x[i] <-  sum(data[grep(files[i],data$file),7:9]*c(60*60,60,1))
}
x <- scale(x)
dim(L) <- c(prod(sz),m)
mad1 <- apply(L,1,median)
kill <- 1.4826*abs(sweep(L,1,mad1,"-"))
eps <- apply(kill,1,median)
keep4 <- eps > 1/3*median(eps[which(eps > 1e-4)])
keep <- which(keep4)
notkeep <- which(!keep4)
kill2 <-array(0,dim=dim(L))
for(i in keep){
  f <- rlm(x,L[i,])
  kill[i,] <- f$residuals
}
kill <- kill^2
maxkill <- apply(kill,2,max)
kill <- sweep(kill,2,maxkill,"/")
dim(kill) <- c(sz,m)
path2 <- sub('processedata','2processedata',path)
for(i in 1:m){
  writePNG(kill[,,i],target=paste(path2,files[i],sep="/"))
}
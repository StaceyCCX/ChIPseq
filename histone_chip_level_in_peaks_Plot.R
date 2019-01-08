rm(list = ls())
setwd("/cluster/home/cgduan/Project/histone_modification/tmp/")
k <- 0
#file_list <- commandArgs()
file_list <- list.files(pattern = 'Raw.*.')
#for (x in file_list) {txt <- read.table(x,skip = 1,encoding = 'UTF-8')
for (x in file_list) {txt <- read.table(x,skip = 1)
data <- data.frame(txt[,7:ncol(txt)])
#mean_colSums=data.frame(round(colSums(data)/nrow(data),2))
k=k+1
name <- paste('data',k,sep='_')
assign(name,data)}

#data <- cbind(get(ls(pattern='Raw.*.')))
data <- rbind(data_1,data_2,data_3,data_4)
pdf('histone_modification_level_in_ASI1_EDM2_peaks.pdf')
for (i in seq(1,10000,by=1000)) { plot(seq(1,1000),data[1,i:(i+999)],type = 'l',ann=F,xlab = 'position',ylab = 'chip level',ylim = c(0,5),xaxt='n',yaxt='n',lwd=5,col='grey')
xlabel=c(-5000,'MidPoint',+5000)
ylabel=c(1,2,3,4,5)
axis(2,at=c(1,2,3,4,5),labels=ylabel,font=2,cex.axis=1.5)
axis(1,at = c(1,500,1000),labels=xlabel,font=2,cex.axis=1.5)
title(xlab='genome position (bps)',ylab = 'ChIP level',cex.lab=1.5,font.lab=2,main='histone ChIP level in peak region')
lines(seq(1,1000),data[2,i:(i+999)],col='blue',lwd=5)
lines(seq(1,1000),data[3,i:(i+999)],col='red',lwd=5)
lines(seq(1,1000),data[4,i:(i+999)],col='green',lwd=5)}
dev.off()

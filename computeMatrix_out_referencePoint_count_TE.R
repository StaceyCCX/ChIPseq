###This R script used to plot line plot of computeMatrixs reference-point output.
###before run this R script,you sould run computeMatrix reference-point to get output gz files,commond like this:
###ls newBW_1_*|while read id;do computeMatrix reference-point --referencePoint TSS -p 8 -b 3000 -a 3000 -bs 50 -R $id -S ~/myProject/chipseq/ASInewBW_1_EDM2_peak_macs2_recall/bw/?????F*bw --missingDataAsZero --skipZeros -out ${id%%.bed}_ChIP_level.gz;done
###note:must be one bed file one gz file.

##################### start works ###################################

rm(list = ls())

txt1=read.table('newBW_1_random_by_total_TE_ChIP_level.gz',skip=1)
txt2=read.table('newBW_1_totalTE_not_overlap_gene_ChIP_level.gz',skip=1)
txt3=read.table('newBW_1_totalTE_which_overlap_gene_ChIP_level.gz',skip=1)

data1=data.frame(txt1[,7:ncol(txt1)])
data2=data.frame(txt2[,7:ncol(txt2)])
data3=data.frame(txt3[,7:ncol(txt3)])

mean_colSums1=data.frame(round(colSums(data1)/nrow(data1),2))
mean_colSums2=data.frame(round(colSums(data2)/nrow(data2),2))
mean_colSums3=data.frame(round(colSums(data3)/nrow(data3),2))

names(mean_colSums1)='y1'
names(mean_colSums2)='y2'
names(mean_colSums3)='y3'

data=cbind(mean_colSums1,mean_colSums2,mean_colSums3)
write.table(data,file='tmp.txt',row.names = FALSE,col.names = FALSE,sep = '\t',)

x=data.frame(seq(1,nrow(data)))
names(x)='x'
y=data
all=cbind(x,y)

pdf('newBW_ChIP_level_from_TE_start.pdf')
plot(all$x[1:120],all$y1[1:120],type = 'l',ann=F,ylim=c(1,5),xlab = 'position',ylab = 'chip level',xaxt='n',yaxt='n',lwd=5,col='grey')
xlabel=c(-3000,'TE start',+3000)
ylabel=c(1,2,3,4,5)
axis(2,at=c(1,2,3,4,5),labels=ylabel,font=2,cex.axis=1.5)
axis(1,at = c(1,60,120),labels=xlabel,font=2,cex.axis=1.5)
title(xlab='genome position (bps)',ylab = 'ChIP level',cex.lab=1.5,font.lab=2,main='EDM2 F1 ChIP level in TE region')
lines(all$y2[1:120],col='blue',lwd=5)
lines(all$y3[1:120],col='red',lwd=5)

plot(all$x[1:120],all$y1[121:240],type = 'l',ann=F,ylim=c(1,5),xlab = 'position',ylab = 'chip level',xaxt='n',yaxt='n',lwd=5,col='grey')
xlabel=c(-3000,'TE start',+3000)
ylabel=c(1,2,3,4,5)
axis(1,at = c(1,60,120),labels=xlabel,font=2,cex.axis=1.5)
axis(2,at=c(1,2,3,4,5),labels=ylabel,font=2,cex.axis=1.5)
title(xlab='genome position (bps)',ylab = 'ChIP level',cex.lab=1.5,font.lab=2,main='EDM2 F2 ChIP level in TE region')
lines(all$y2[121:240],col='blue',lwd=5)
lines(all$y3[121:240],col='red',lwd=5)

plot(all$x[1:120],all$y1[241:360],type = 'l',ann=F,ylim=c(1,5),xlab = 'position',ylab = 'chip level',xaxt='n',yaxt='n',lwd=5,col='grey')
xlabel=c(-3000,'TE start',+3000)
ylabel=c(1,2,3,4,5)
axis(1,at = c(1,60,120),labels=xlabel,font=2,cex.axis=1.5)
axis(2,at=c(1,2,3,4,5),labels=ylabel,font=2,cex.axis=1.5)
title(xlab='genome position (bps)',ylab = 'ChIP level',cex.lab=1.5,font.lab=2,main='ASI1 F1 ChIP level in TE region')
lines(all$y2[241:360],col='blue',lwd=5)
lines(all$y3[241:360],col='red',lwd=5)

plot(all$x[1:120],all$y1[361:480],type = 'l',ann=F,ylim=c(1,5),xlab = 'position',ylab = 'chip level',xaxt='n',yaxt='n',lwd=5,col='grey')
xlabel=c(-3000,'TE start',+3000)
ylabel=c(1,2,3,4,5)
title(xlab='genome position (bps)',ylab = 'ChIP level',cex.lab=1.5,font.lab=2,main='ASI1 F2 ChIP level in TE region')
axis(1,at = c(1,60,120),labels=xlabel,font=2,cex.axis=1.5)
axis(2,at=c(1,2,3,4,5),labels=ylabel,font=2,cex.axis=1.5)
lines(all$y2[361:480],col='blue',lwd=5)
lines(all$y3[361:480],col='red',lwd=5)

plot(all$x[1:120],all$y1[481:600],type = 'l',ann=F,ylim=c(1,5),xlab = 'position',ylab = 'chip level',xaxt='n',yaxt='n',lwd=5,col='grey')
xlabel=c(-3000,'TE start',+3000)
ylabel=c(1,2,3,4,5)
axis(1,at = c(1,60,120),labels=xlabel,font=2,cex.axis=1.5)
axis(2,at=c(1,2,3,4,5),labels=ylabel,font=2,cex.axis=1.5)
title(xlab='genome position (bps)',ylab = 'ChIP level',cex.lab=1.5,font.lab=2,main='edm2 F1 ChIP level in TE region')
lines(all$y2[481:600],col='blue',lwd=5)
lines(all$y3[481:600],col='red',lwd=5)

plot(all$x[1:120],all$y1[601:720],type = 'l',ann=F,ylim=c(1,5),xlab = 'position',ylab = 'chip level',xaxt='n',yaxt='n',lwd=5,col='grey')
xlabel=c(-3000,'TE start',+3000)
ylabel=c(1,2,3,4,5)
axis(1,at = c(1,60,120),labels=xlabel,font=2,cex.axis=1.5)
axis(2,at=c(1,2,3,4,5),labels=ylabel,font=2,cex.axis=1.5)
title(xlab='genome position (bps)',ylab = 'ChIP level',cex.lab=1.5,font.lab=2,main='edm2 F2 ChIP level in TE region')
lines(all$y2[601:720],col='blue',lwd=5)
lines(all$y3[601:720],col='red',lwd=5)
dev.off()

###This R script used to plot line plot of computeMatrixs scale-regions output.
###before run this R script,you sould run computeMatrix scale-regions to get output gz files,commond like this:
###ls 2_*bed|while read id;do computeMatrix scale-regions -p 8 -bs 50 -R $id -S ~/myProject/chipseq/ASI1_EDM2_peak_macs2_recall/bw/?????F*bw  -b 2000 -a 2000 --regionBodyLength 5000 --missingDataAsZero --skipZeros -o ${id%%.bed}_ChIP_for_TEcontainingGene.gz;done
###note: must be one bed file one gz file

##################### start works ###################################
txt1=read.table('2_random_select_gene_ChIP_for_TEcontainingGene.gz',skip=1)
txt2=read.table('2_nonTE_overlaped_gene_ChIP_for_TEcontainingGene.gz',skip=1)
txt3=read.table('2_TE_overlaped_gene_ChIP_for_TEcontainingGene.gz',skip=1)

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

pdf('newBW_ChIP_level_in_TEcontainingGene.pdf')
plot(all$x[1:180],all$y1[1:180],type = 'l',ann=F,ylim=c(1,5),xlab = 'position',ylab = 'chip level',xaxt='n',yaxt='n',lwd=5,col='grey')
xlabel=c('-2000','TSS','TTS','+2000')
ylabel=c(1,2,3,4,5)
axis(1,at = c(1,40,140,180),labels=xlabel,font=2,cex.axis=1.5)
axis(2,at=c(1,2,3,4,5),labels=ylabel,font=2,cex.axis=1.5)
title(xlab='genome position (bps)',ylab = 'ChIP level',cex.lab=1.5,font.lab=2,main='EDM2 F1 ChIP level in TE region')
lines(all$y2,col='blue',lwd=5)
lines(all$y3,col='red',lwd=5)

plot(all$x[1:180],all$y1[181:360],type = 'l',ann=F,ylim=c(1,5),xlab = 'position',ylab = 'chip level',xaxt='n',yaxt='n',lwd=5,col='grey')
title(xlab='genome position (bps)',ylab = 'ChIP level',cex.lab=1.5,font.lab=2,main='EDM2 F2 ChIP level in TE region')
xlabel=c('-2000','TSS','TTS','+2000')
ylabel=c(1,2,3,4,5)
axis(1,at = c(1,40,140,180),labels=xlabel,font=2,cex.axis=1.5)
axis(2,at=c(1,2,3,4,5),labels=ylabel,font=2,cex.axis=1.5)
lines(all$y2[181:360],col='blue',lwd=5)
lines(all$y3[181:360],col='red',lwd=5)

plot(all$x[1:180],all$y1[361:540],type = 'l',ann=F,ylim=c(1,5),xlab = 'position',ylab = 'chip level',xaxt='n',yaxt='n',lwd=5,col='grey')
xlabel=c('-2000','TSS','TTS','+2000')
ylabel=c(1,2,3,4,5)
axis(1,at = c(1,40,140,180),labels=xlabel,font=2,cex.axis=1.5)
axis(2,at=c(1,2,3,4,5),labels=ylabel,font=2,cex.axis=1.5)
title(xlab='genome position (bps)',ylab = 'ChIP level',cex.lab=1.5,font.lab=2,main='ASI1 F1 ChIP level in TE region')
lines(all$y2[361:540],col='blue',lwd=5)
lines(all$y3[361:540],col='red',lwd=5)

plot(all$x[1:180],all$y1[541:720],type = 'l',ann=F,ylim=c(1,5),xlab = 'position',ylab = 'chip level',xaxt='n',yaxt='n',lwd=5,col='grey')
xlabel=c('-2000','TSS','TTS','+2000')
ylabel=c(1,2,3,4,5)
axis(1,at = c(1,40,140,180),labels=xlabel,font=2,cex.axis=1.5)
axis(2,at=c(1,2,3,4,5),labels=ylabel,font=2,cex.axis=1.5)
title(xlab='genome position (bps)',ylab = 'ChIP level',cex.lab=1.5,font.lab=2,main='ASI1 F2 ChIP level in TE region')
lines(all$y2[541:720],col='blue',lwd=5)
lines(all$y3[541:720],col='red',lwd=5)

plot(all$x[1:180],all$y1[721:900],type = 'l',ann=F,ylim=c(1,5),xlab = 'position',ylab = 'chip level',xaxt='n',yaxt='n',lwd=5,col='grey')
xlabel=c('-2000','TSS','TTS','+2000')
ylabel=c(1,2,3,4,5)
axis(1,at = c(1,40,140,180),labels=xlabel,font=2,cex.axis=1.5)
axis(2,at=c(1,2,3,4,5),labels=ylabel,font=2,cex.axis=1.5)
title(xlab='genome position (bps)',ylab = 'ChIP level',cex.lab=1.5,font.lab=2,main='edm2 F1 ChIP level in TE region')
lines(all$y2[721:900],col='blue',lwd=5)
lines(all$y3[721:900],col='red',lwd=5)

plot(all$x[1:180],all$y1[901:1080],type = 'l',ann=F,ylim=c(1,5),xlab = 'position',ylab = 'chip level',xaxt='n',yaxt='n',lwd=5,col='grey')
xlabel=c('-2000','TSS','TTS','+2000')
ylabel=c(1,2,3,4,5)
axis(1,at = c(1,40,140,180),labels=xlabel,font=2,cex.axis=1.5)
axis(2,at=c(1,2,3,4,5),labels=ylabel,font=2,cex.axis=1.5)
title(xlab='genome position (bps)',ylab = 'ChIP level',cex.lab=1.5,font.lab=2,main='edm2 F2 ChIP level in TE region')
lines(all$y2[901:1080],col='blue',lwd=5)
lines(all$y3[901:1080],col='red',lwd=5)
dev.off()

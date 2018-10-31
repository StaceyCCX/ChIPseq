#computeMatrix scale-regions 
#-bs 100 -R tair10.chr1.sizes.bed -S /home/galaxy/myProject/chipseq/sicer_call_bw/*bw --regionBodyLength 163700 -b 0 -a 0 --missingDataAsZero --skipZeros 
#-o computeMatrix_out_TAIR10_chr1.gz

library(data.table)

txt<-fread('computeMatrix_out_TAIR10_chr1.gz',header=FALSE,skip = 1)
chr<-data.frame(txt[,7:ncol(txt)])
chr<-as.numeric(chr)
sample_name <- c('EDM2_WT_F1','EDM2_WT_F2','EDM2_WT_input','ASI1_WT_F1','ASI1_WT_F2','ASI1_WT_input','EDM2_F1','EDM2_F2','EDM2_input','ASI1_F1','ASI1_F2','ASI1_input','edm2_F1','edm2_F2','edm2_input')
x1=c(1,1638,3275,4912,6549,8186,9823,11460,13097,14734,16371,18008,19645,21282,22919)
x2=c(1637,3274,4911,6548,8185,9822,11459,13096,14733,16370,18007,19644,21281,22918,24555)
r <- data.frame(sample_name,x1,x2,stringsAsFactors = F)
for (i in 1:nrow(r)){
	cut_table_name <- r$sample_name[i]
	assign(cut_table_name,c(chr[r[i,2]:r[i,3]]))
	}

count_EDM2_F1<-data.frame(EDM2_F1,EDM2_input,EDM2_WT_F1,EDM2_WT_input)
names(count_EDM2_F1) <- c('chip','chip_input','wt','wt_input')
count_EDM2_F1$chip_signal <- log2(count_EDM2_F1$chip/count_EDM2_F1$chip_input)-log2(count_EDM2_F1$wt/count_EDM2_F1$wt_input)

count_EDM2_F2<-data.frame(EDM2_F2,EDM2_input,EDM2_WT_F2,EDM2_WT_input)
names(count_EDM2_F2) <- c('chip','chip_input','wt','wt_input')
count_EDM2_F2$chip_signal <- log2(count_EDM2_F2$chip/count_EDM2_F2$chip_input)-log2(count_EDM2_F2$wt/count_EDM2_F2$wt_input)

count_ASI1_F1<-data.frame(ASI1_F1,ASI1_input,ASI1_WT_F1,ASI1_WT_input)
names(count_ASI1_F1) <- c('chip','chip_input','wt','wt_input')
count_ASI1_F1$chip_signal <- log2(count_ASI1_F1$chip/count_ASI1_F1$chip_input)-log2(count_ASI1_F1$wt/count_ASI1_F1$wt_input)

count_ASI1_F2<-data.frame(ASI1_F2,ASI1_input,ASI1_WT_F2,ASI1_WT_input)
names(count_ASI1_F2) <- c('chip','chip_input','wt','wt_input')
count_ASI1_F2$chip_signal <- log2(count_ASI1_F2$chip/count_ASI1_F2$chip_input)-log2(count_ASI1_F2$wt/count_ASI1_F2$wt_input)

count_edm2_F1<-data.frame(edm2_F1,edm2_input,ASI1_WT_F1,ASI1_WT_input)
names(count_edm2_F1) <- c('chip','chip_input','wt','wt_input')
count_edm2_F1$chip_signal <- log2(count_edm2_F1$chip/count_edm2_F1$chip_input)-log2(count_edm2_F1$wt/count_edm2_F1$wt_input)


count_edm2_F2<-data.frame(edm2_F2,edm2_input,ASI1_WT_F2,ASI1_WT_input)
names(count_edm2_F2) <- c('chip','chip_input','wt','wt_input')
count_edm2_F2$chip_signal <- log2(count_edm2_F2$chip/count_edm2_F2$chip_input)-log2(count_edm2_F2$wt/count_edm2_F2$wt_input)

chrom_size<-30427671
normal_size<-1637 #len(longest_chr)/len(shortest_chr) ,close to ?00,devide bin size 100,eg. 30427671/18585056=1.637212,so normal_size=163700/100=1637
center_pos<-14811145
centromere<-round(center_pos/(chrom_size/normal_size))

pdf("computeMatrix_out_TAIR10_chr1.pdf")

plot(1:nrow(count_EDM2_F1),count_EDM2_F1$chip_signal,type = 'l',ylim=c(-10,10),ann	=F,lwd=5,col='black')
xlabel='centromere'
axis(1,at = centromere,labels=xlabel,font=2,cex.axis=1.5)
title(xlab='genome position (bps)',ylab = 'ChIP signal',ylim=c(-10,10),cex.lab=1.5,font.lab=2,main='EDM2 F1 ChIP signal in genomic region')

plot(1:nrow(count_EDM2_F2),count_EDM2_F2$chip_signal,type = 'l',ylim=c(-10,10),ann	=F,lwd=5,col='black')
xlabel='centromere'
axis(1,at = centromere,labels=xlabel,font=2,cex.axis=1.5)
title(xlab='genome position (bps)',ylab = 'ChIP signal',ylim=c(-10,10),cex.lab=1.5,font.lab=2,main='EDM2 F2 ChIP signal in genomic region')

plot(1:nrow(count_ASI1_F1),count_ASI1_F1$chip_signal,type = 'l',ylim=c(-10,10),ann	=F,lwd=5,col='black')
xlabel='centromere'
axis(1,at = centromere,labels=xlabel,font=2,cex.axis=1.5)
title(xlab='genome position (bps)',ylab = 'ChIP signal',ylim=c(-10,10),cex.lab=1.5,font.lab=2,main='ASI1 F1 ChIP signal in genomic region')

plot(1:nrow(count_ASI1_F2),count_ASI1_F2$chip_signal,type = 'l',ylim=c(-10,10),ann	=F,lwd=5,col='black')
xlabel='centromere'
axis(1,at = centromere,labels=xlabel,font=2,cex.axis=1.5)
title(xlab='genome position (bps)',ylab = 'ChIP signal',ylim=c(-10,10),cex.lab=1.5,font.lab=2,main='ASI1 F2 ChIP signal in genomic region')

plot(1:nrow(count_edm2_F1),count_edm2_F1$chip_signal,type = 'l',ylim=c(-10,10),ann	=F,lwd=5,col='black')
xlabel='centromere'
axis(1,at = centromere,labels=xlabel,font=2,cex.axis=1.5)
title(xlab='genome position (bps)',ylab = 'ChIP signal',ylim=c(-10,10),cex.lab=1.5,font.lab=2,main='ASI1/edm2 F1 ChIP signal in genomic region')

plot(1:nrow(count_edm2_F2),count_edm2_F2$chip_signal,type = 'l',ylim=c(-10,10),ann	=F,lwd=5,col='black')
xlabel='centromere'
axis(1,at = centromere,labels=xlabel,font=2,cex.axis=1.5)
title(xlab='genome position (bps)',ylab = 'ChIP signal',ylim=c(-10,10),cex.lab=1.5,font.lab=2,main='ASI1/edm2 F1 ChIP signal in genomic region')

dev.off()

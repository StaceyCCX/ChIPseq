### nl test.bed>test.saf # change bed file into saf annotation file
### featureCounts -T 20 -p -C -F SAF -a test.saf  -o test.txt *bam # count fragment in saf regions 
rm(list = ls())

filelist=list.files(pattern = 'featureCounts_right_.*.txt$')
datalist=lapply(filelist, function(x) read.table(x,header=F,stringsAsFactors = F))
WT1=function(x){
	x[7]
}
WT2=function(x){
	x[8]
}
EDM2=function(x){
	x[9]
}
ASI1=function(x){
	x[10]
}

WT1_list=sapply(datalist,WT1)
WT2_list=sapply(datalist,WT2)
EDM2_list=sapply(datalist,EDM2)
ASI1_list=sapply(datalist,ASI1)

Total_frag_WT1 <- 3637627
Total_frag_WT2 <- 3745665
Total_frag_EDM2 <- 4167406
Total_frag_ASI1 <- 10032809

nor_WT1_list <- sapply(WT1_list, function(x) as.numeric(x)/Total_frag_WT1*1000000)
nor_WT2_list <- sapply(WT2_list, function(x) as.numeric(x)/Total_frag_WT2*1000000)
nor_EDM2_list <- sapply(EDM2_list, function(x) as.numeric(x)/Total_frag_EDM2*1000000)
nor_ASI1_list <- sapply(ASI1_list, function(x) as.numeric(x)/Total_frag_ASI1*1000000)


pdf("boxplot_in_ASI1specific_EDM2specific_ASI1overlapEDM2.pdf")

ASI1_specific <- cbind(nor_WT1_list[[1]],nor_WT2_list[[1]],nor_EDM2_list[[1]],nor_ASI1_list[[1]])
n <- nrow(ASI1_specific)
set1 <- c(ASI1_specific[2:n,1],ASI1_specific[2:n,2],ASI1_specific[2:n,3],ASI1_specific[2:n,4])
fset1 <- factor(c(rep(1,length(ASI1_specific[2:n,1])),rep(2,length(ASI1_specific[2:n,2])),rep(3,length(ASI1_specific[2:n,3])),rep(4,length(ASI1_specific[2:n,4]))))
plot(fset1,set1,xlab='Samples',ylab='ChIP Fragment Counts',ylim=c(0,100),outline=F,xaxt='n',col=c('grey','grey','blue','green'))
xlabel=c('WT/EDM2','WT/ASI1','EDM2','ASI1')
axis(1,at = c(1,2,3,4),labels=xlabel,font=2,cex.axis=1)
title(cex.lab=1,font.lab=2,main='WT EDM2 ASI1 ChIP fragment count in ASI1 specific region')

EDM2_specific <- cbind(nor_WT1_list[[2]],nor_WT2_list[[2]],nor_EDM2_list[[2]],nor_ASI1_list[[2]])
n <- nrow(EDM2_specific)
set2 <- c(EDM2_specific[2:n,1],EDM2_specific[2:n,2],EDM2_specific[2:n,3],EDM2_specific[2:n,4])
fset2 <- factor(c(rep(1,length(EDM2_specific[2:n,1])),rep(2,length(EDM2_specific[2:n,2])),rep(3,length(EDM2_specific[2:n,3])),rep(4,length(EDM2_specific[2:n,4]))))
plot(fset2,set2,xlab='Samples',ylab='ChIP Fragment Counts',ylim=c(0,100),outline=F,xaxt='n',col=c('grey','grey','blue','green'))
xlabel=c('WT/EDM2','WT/ASI1','EDM2','ASI1')
axis(1,at = c(1,2,3,4),labels=xlabel,font=2,cex.axis=1)
title(cex.lab=1,font.lab=2,main='WT EDM2 ASI1 ChIP fragment count in EDM2 specific region')


ASI1_overlap_EDM2 <- cbind(nor_WT1_list[[3]],nor_WT2_list[[3]],nor_EDM2_list[[3]],nor_ASI1_list[[3]])
n <- nrow(ASI1_overlap_EDM2)
set3 <- c(ASI1_overlap_EDM2[2:n,1],ASI1_overlap_EDM2[2:n,2],ASI1_overlap_EDM2[2:n,3],ASI1_overlap_EDM2[2:n,4])
fset3 <- factor(c(rep(1,length(ASI1_overlap_EDM2[2:n,1])),rep(2,length(ASI1_overlap_EDM2[2:n,2])),rep(3,length(ASI1_overlap_EDM2[2:n,3])),rep(4,length(ASI1_overlap_EDM2[2:n,4]))))
plot(fset3,set3,xlab='Samples',ylab='ChIP Fragment Counts',ylim=c(0,100),outline=F,xaxt='n',col=c('grey','grey','blue','green'))
xlabel=c('WT/EDM2','WT/ASI1','EDM2','ASI1')
axis(1,at = c(1,2,3,4),labels=xlabel,font=2,cex.axis=1)
title(cex.lab=1,font.lab=2,main='WT EDM2 ASI1 ChIP fragment count in ASI1 overlap EDM2 region')

dev.off()

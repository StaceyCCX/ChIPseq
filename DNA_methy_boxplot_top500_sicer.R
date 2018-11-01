rm (list=ls())
setwd("/home/galaxy/myProject/DNA_methy_for_ASI1_ChIP/top500_sicer_ASI1_EDM2_methy_2reps")

filelist=list.files(pattern = '.*.rep2.mC.pv.txt') 
datalist=lapply(filelist, function(x) read.table(x,header=F)) 
subdata=function(x){
    x[4]
    } 
mC_list=sapply(datalist,subdata)
mC=c(mC_list[[1]],mC_list[[4]],mC_list[[7]],
	mC_list[[2]],mC_list[[5]],mC_list[[8]],
	mC_list[[3]],mC_list[[6]],mC_list[[9]])
fmC=factor(c(rep(1,length(mC_list[[1]])),rep(2,length(mC_list[[4]])),rep(3,length(mC_list[[7]])),
		rep(4,length(mC_list[[2]])),rep(5,length(mC_list[[5]])),rep(6,length(mC_list[[8]])),
		rep(7,length(mC_list[[3]])),rep(8,length(mC_list[[6]])),rep(9,length(mC_list[[9]]))
		))
pdf('/home/galaxy/myProject/DNA_methy_for_ASI1_ChIP/top500_sicer_ASI1_EDM2_methy_2reps_rep2_mC.pdf')
plot(fmC,mC,xlab='Samples',ylab='mC methylation level',outline=F,
                                col=c(rgb(0.54,0.59,0.48),rgb(0.54,0.59,0.48),rgb(0.54,0.59,0.48),
                                rgb(0.96,0.81,0),rgb(0.96,0.81,0),rgb(0.96,0.81,0),
                                rgb(0.90,0.51,0.03),rgb(0.90,0.51,0.03),rgb(0.90,0.51,0.03)))
dev.off()

# lgt_contig_length_boxplot.R

library(ggplot2)

bac_len<-c(119,76,101,78,74,126,136,146,97,273,109,183,125)
hg_len<-c(95,52,74,73,77,61,68,91,58,96,85,108,78)

df<-data.frame("bac",bac_len)
df2<-data.frame("hg",hg_len)

colnames(df)<-c("contig","length")
colnames(df2)<-c("contig","length")

total<-rbind(df,df2)

pdf(file="/Users/ksieber/Documents/lgt_contig_length.pdf")
qplot( factor(contig), length, data=total, geom="boxplot", xlab="contig", ylab="length", fill=factor(contig) ) + theme(legend.position="none")

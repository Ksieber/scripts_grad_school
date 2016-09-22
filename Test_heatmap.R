# Load libraries and source code
library(pvclust)
library(snow)
source("/Users/ksieber/scripts/dendroCol.R")
source("/Users/ksieber/scripts/heatmap.3.R")
source("/Users/ksieber/scripts/Heatmap_colors.R")
source("/Users/ksieber/scripts/colorTbl.R")

# Loading data from grid
MICRO<-read.table(file="/local/projects-t3/HLGT/TCGA/ksieber_dir/COAD/all_coad_merge2/test_LGTSeek2heatmaps/test_Rtable_micro.txt", header = TRUE, na.strings = "NA", sep="\t")
# Load data from MBP
MICRO<-read.table(file="/Users/ksieber/HLGT/COAD/all_coad_merge2/test_LGTSeek2heatmaps/test_Rtable_micro.txt", header=TRUE, sep="\t", na.strings="NA")

# Manipulate data into proper format
drop_rownames<-c("Reads")
MICRO2<-as.matrix(MICRO[,!(names(MICRO) %in% drop_rownames)])
rownames(MICRO2)<-MICRO[,1]
MICRO2[is.na(MICRO2)]<-0
MICRO_proportion<-t(apply(MICRO2, 1, function(x) 100 * x/sum(x)))

# Cluster with hclust and pvclust
hc<-hclust(dist(as.matrix(MICRO_proportion)))
# Single core:
pv<-pvclust(t(MICRO_proportion), nboot=1000)
## Alternative: multi-thread pvclust
cl<-makeCluster(3,type="MPI")
pv<-parPvclust(cl, t(MICRO_proportion), nboot=1000)
stopCluster(cl)

# Cut the clustering trees in multiple ways for comparison
cut_hc_K<-cutree(hc, k=20)
cut_hc_H<-cutree(hc, h=(max(hc$height)/1.5))
cut_pv_K<-cutree(pv$hclust, k=20)
cut_pv_H<-cutree(pv$hclust, h=(max(pv$hclust$height)/1.5))

# Create sidebars for heatmaps
hc_K_SB<-colorTbl[as.vector(cut_hc_K)]
hc_H_SB<-colorTbl[as.vector(cut_hc_H)]
pv_K_SB<-colorTbl[as.vector(cut_pv_K)]
pv_H_SB<-colorTbl[as.vector(cut_pv_H)]
all_SB<-cbind(hc_K_SB,hc_H_SB,pv_K_SB,pv_H_SB)
colnames(all_SB)<-c("hc_K", "hc_H", "pv_K", "pv_H")
all_SB<-as.matrix(t(all_SB))

# Create colored dendrogram for pvclust significant clusters
pv_sig_keys<- unlist(pvpick(pv, alpha=0.95, pv="au", type="geq", max.only=TRUE)$clusters)
dendrogram_colored <- dendrapply(as.dendrogram(pv$hclust), dendroCol, keys=pv_sig_keys, xPar="edgePar", bgr="black", fgr="red", pch=20)

# Heatmap for cluster with pvclust
pdf(file="/Users/ksieber/HLGT/COAD/all_coad_merge2/test_LGTSeek2heatmaps/test_pv-v-hc_pv-clustered_k20.pdf")
heatmap.3(MICRO_proportion, col=COLORS, Rowv=dendrogram_colored, RowSideColors=all_SB, RowSideColorsSize=4, na.color="#999999", trace="none", keysize=1, key=1	)
dev.off()

# Heatmap for cluster with hclust
pdf(file="/Users/ksieber/HLGT/COAD/all_coad_merge2/test_LGTSeek2heatmaps/test_pv-v-hc_hc-clustered_k20.pdf")
heatmap.3(MICRO_proportion, col=COLORS, Rowv=as.dendrogram(hc),	RowSideColors=all_SB, RowSideColorsSize=4, na.color="#999999", trace="none", keysize=1, key=1	)
dev.off()

# Heatmap of log transformed
# Run this for ln(proportion) of each OTU:
MICRO_log<-t(apply(MICRO2, 1, function(x) log(100 * x/sum(x))))
MICRO_log[!is.finite(MICRO_log)]<-NA
otu_cluster<-hclust(dist(as.matrix(t(MICRO_proportion))))
heatmap.3(MICRO_log, Rowv=as.dendrogram(hc), col=COLORS, na.color="#999999", Colv=as.dendrogram(otu_cluster))


# # Same thing, over multi-lines
# # pvclust
# heatmap.3(MICRO_proportion, 
# 	col=COLORS, 
# 	Rowv=dendrogram_colored,
# 	RowSideColors=all_SB,
# 	RowSideColorsSize=4,
# 	na.color="#999999",
# 	trace="none",
# 	keysize=1,
# 	key=1
# 	)
# # hclust:
# heatmap.3(MICRO_proportion, )
# 	col=COLORS, 
# 	Rowv=as.dendrogram(hc),
# 	RowSideColors=all_SB,
# 	RowSideColorsSize=4,
# 	na.color="#999999",
# 	trace="none",
# 	keysize=1,
# 	key=1
# 	)
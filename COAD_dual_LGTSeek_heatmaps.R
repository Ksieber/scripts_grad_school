# COAD_dual_LGTSeek_heatmaps !!! NOT FINISHED !!!
# Load libraries and source code
library(pvclust)
# library(snow)
source("/Users/ksieber/scripts/dendroCol.R")
source("/Users/ksieber/scripts/heatmap.3.R")
source("/Users/ksieber/scripts/Heatmap_colors.R")
source("/Users/ksieber/scripts/colorTbl.R")

## This fxn compares matrix X to Y, and creates COLUMNS that are 
## filled with NA in X that are missing compared to Y.
## Returns X. 
fill_empty_OTU_fxn <- function (x,y){ # X = filled in, relative to Y.
	for (i in 1:length(colnames(y))){ # foreach OTU in Y
		if(colnames(y)[i] %in% colnames(x)==FALSE){	# If the OTU isn't in X but is in Y
			empty<-matrix(data=NA, nrow=length(rownames(x)), ncol=1)			
			colnames(empty)<-c(colnames(y)[i])
			x<-cbind(x,empty)			
		}
	}
	return(x)
}

## This fxn compares matrix X to Y, and creates ROWS that are 
## filled with NA in X that are missing compared to Y.
## Returns X. 
fill_empty_ROW_fxn <- function (x,y){ 
	for (i in 1:length(rownames(y))){ 
		if(rownames(y)[i] %in% rownames(x)==FALSE){	
			empty<-matrix(data=NA, ncol=length(colnames(x)), nrow=1)			
			rownames(empty)<-c(rownames(y)[i])
			x<-rbind(x,empty)			
		}
	}
	return(x)
}

# Load data from MBP
MICRO<-read.table(file="/Users/ksieber/HLGT/COAD/all_coad_merge2/test_LGTSeek2heatmaps/test_Rtable_micro.txt", header=TRUE, sep="\t", na.strings="NA")
LGT<-read.table(file="/Users/ksieber/HLGT/COAD/all_coad_merge2/test_LGTSeek2heatmaps/test_Rtable_lgt.txt", header=TRUE, sep="\t", na.strings="NA")

MICRO<-read.table(file="/Users/ksieber/HLGT/COAD/all_coad_merge2/dual_heatmaps/input_rm-no-lca/test_micro_Rtable.txt", header=TRUE, sep="\t", na.strings="NA", stringsAsFactors=FALSE)
LGT<-read.table(file="/Users/ksieber/HLGT/COAD/all_coad_merge2/dual_heatmaps/input_rm-no-lca/test_lgt_Rtable.txt", header=TRUE, sep="\t", na.strings="NA", stringsAsFactors=FALSE)

# Manipulate data into proper format
# Microbiome first
drop_rownames<-c("Ids")
MICRO2<-as.matrix(MICRO[,!(names(MICRO) %in% drop_rownames)])
MICRO3<-matrix(as.numeric(unlist(MICRO2)), nrow=nrow(MICRO2))
rownames(MICRO3)<-MICRO[,1]
colnames(MICRO3)<-colnames(MICRO2)
# LGT
LGT2<-as.matrix(LGT[,!(names(LGT) %in% drop_rownames)])
LGT3<-matrix(as.numeric(unlist(LGT2)), nrow=nrow(LGT2))
rownames(LGT3)<-LGT[,1]
colnames(LGT3)<-colnames(LGT2)

# Ensure all matrices have all the data.
LGT3<-fill_empty_OTU_fxn(LGT3, MICRO3)
MICRO3<-fill_empty_OTU_fxn(MICRO3,LGT3)
LGT3<-fill_empty_ROW_fxn(LGT3,MICRO3)
MICRO3<-fill_empty_ROW_fxn(MICRO3,LGT3)

MICRO_counts <- MICRO3
MICRO_counts[is.na(MICRO_counts)] <- 0
LGT_counts <- LGT3
LGT_counts[is.na(LGT_counts)] <- 0

# Cut MICRO and LGT data to subsets of the top # of OTU
## 1. Cut based on MICRO abundance
top_MICRO_col_list <- order(colSums(MICRO_counts), decreasing="T")[1:20]
top_MICRO <- MICRO_counts[,top_MICRO_col_list]
top_LGT <- LGT_counts[,top_MICRO_col_list]

## 2. Cut based on LGT abundance instead
# top_LGT_col_list <- order(colSums(LGT_counts), decreasing="T")[1:20]
# top_MICRO <- MICRO_counts[,top_LGT_col_list]
# top_LGT <- LGT_counts[,top_LGT_col_list]

# Calculate proportion of each OTU / participant
MICRO_proportion <- (top_MICRO/rowSums(top_MICRO))*100
LGT_proportion <- (top_LGT/rowSums(top_LGT))*100

# Some participants may have no LGT/microbiome reads
# This masks them from proportion calc where they = NA
MICRO_proportion_for_cluster <- MICRO_proportion
MICRO_proportion_for_cluster[is.na(MICRO_proportion_for_cluster)]<-0
LGT_proportion_for_cluster <- LGT_proportion
LGT_proportion_for_cluster[is.na(LGT_proportion_for_cluster)]<-0

# Cluster with hclust by MICRO
hc<-hclust(dist(as.matrix(MICRO_proportion_for_cluster)))

# Cluster by LGT
hc<-hclust(dist(as.matrix(LGT_proportion_for_cluster)))

# Cut with set number of clusters
cut_hc<-cutree(hc, k=4)
cluster_colors<-colorTbl[as.vector(cut_hc)]
hc_SB<-as.matrix(t(cluster_colors))
rownames(hc_SB)<-c("Clustering")

# Cut based on relative size of clustering
# cut_hc<-cutree(hc, h=(max(hc$height)/1.5))
# cluster_colors<-colorTbl[as.vector(cut_hc)]
# hc_SB<-t(as.matrix(cluster_colors))
# rownames(hc_SB)<-c("Clustering")

# # Ln transform the data
MICRO_for_Log <- MICRO_proportion
MICRO_for_Log[ MICRO_for_Log == 0] <- NA
MICRO_Log <- log(MICRO_for_Log)
MICRO_Log[is.infinite(MICRO_Log)]<-NA
LGT_for_Log <- LGT_proportion
LGT_for_Log[ LGT_for_Log == 0] <- NA
LGT_Log <- log(LGT_proportion_for_cluster)
LGT_Log[is.infinite(LGT_Log)]<-NA

# # Log10 transform:
# MICRO_for_Log <- MICRO_proportion
# MICRO_for_Log[ MICRO_for_Log == 0] <- NA
# MICRO_Log <- log10(MICRO_for_Log)
# MICRO_Log[is.infinite(MICRO_Log)]<-NA
# LGT_for_Log <- LGT_proportion
# LGT_for_Log[ LGT_for_Log == 0] <- NA
# LGT_Log <- log10(LGT_proportion_for_cluster)
# LGT_Log[is.infinite(LGT_Log)]<-NA

# Heatmap:
heatmap.3(top_MICRO, col=COLORS, Rowv=as.dendrogram(hc), RowSideColors=hc_SB, RowSideColorsSize=4, na.color="#999999", trace="none", keysize=1, key=1 )
heatmap.3(MICRO_proportion, col=COLORS, Rowv=as.dendrogram(hc), RowSideColors=hc_SB, RowSideColorsSize=4, na.color="#999999", trace="none", keysize=1, key=1 )
heatmap.3(MICRO_Log, col=COLORS, Rowv=as.dendrogram(hc), RowSideColors=hc_SB, RowSideColorsSize=4, na.color="#999999", trace="none", keysize=1, key=1 )

# Heatmap to PDF
pdf(file="/Users/ksieber/HLGT/COAD/all_coad_merge2/test_LGTSeek2heatmaps/test_pv-v-hc_hc-clustered_k20.pdf")
heatmap.3(MICRO_proportion, col=COLORS, Rowv=as.dendrogram(hc), RowSideColors=hc_SB, RowSideColorsSize=4, na.color="#999999", trace="none", keysize=1, key=1 )
dev.off()

# LGT Heatmap
heatmap.3(LGT_proportion, col=COLORS, Rowv=as.dendrogram(hc), RowSideColors=hc_SB, RowSideColorsSize=4, na.color="#999999", trace="none", keysize=1, key=1 )
heatmap.3(LGT_Log, col=COLORS, Rowv=as.dendrogram(hc), RowSideColors=hc_SB, RowSideColorsSize=4, na.color="#999999", trace="none", keysize=1, key=1 )

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


## OLD pvclust code:
# Single core:
# pv<-pvclust(t(MICRO_proportion_for_cluster), nboot=100, method.hclust="complete", use.cor="all.obs")
# ## Alternative: multi-thread pvclust
# cl<-makeCluster(3,type="MPI")
# pv<-parPvclust(cl, t(MICRO_proportion), nboot=1000)
# stopCluster(cl)

# cut_pv_K<-cutree(pv$hclust, k=5)
# cut_pv_H<-cutree(pv$hclust, h=(max(pv$hclust$height)/1.5))

# pv_K_SB<-colorTbl[as.vector(cut_pv_K)]
# pv_H_SB<-colorTbl[as.vector(cut_pv_H)]
# all_SB<-cbind(hc_K_SB,hc_H_SB,pv_K_SB,pv_H_SB)
# colnames(all_SB)<-c("hc_K", "hc_H", "pv_K", "pv_H")
# all_SB<-as.matrix(t(all_SB))

# Create colored dendrogram for pvclust significant clusters
# pv_sig_keys<- unlist(pvpick(pv, alpha=0.95, pv="au", type="geq", max.only=TRUE)$clusters)
# dendrogram_colored <- dendrapply(as.dendrogram(pv$hclust), dendroCol, keys=pv_sig_keys, xPar="edgePar", bgr="black", fgr="red", pch=20)
# Heatmap for cluster with pvclust
# pdf(file="/Users/ksieber/HLGT/COAD/all_coad_merge2/test_LGTSeek2heatmaps/test_pv-v-hc_pv-clustered_k20.pdf")
# heatmap.3(MICRO_proportion, col=COLORS, Rowv=dendrogram_colored, RowSideColors=all_SB, RowSideColorsSize=4, na.color="#999999", trace="none", keysize=1, key=1	)
# dev.off()
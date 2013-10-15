
library("ALL", lib.loc="/home/ksieber/R/x86_64-unknown-linux-gnu-library/2.14");
library(gplots);

LCAMatrix=read.table(file="/local/projects-t3/HLGT/TCGA/ksieber_dir/Heatmap_LGT-vs-BacMM/LGT/STAD/4-new_lca_bacMM/lca_matrix_files/STAD_combined.sam2lca.liberal_family_lca.matrix", header = TRUE, na.strings = "NA");
tLCA=t(LCAMatrix);   	    ##transpose my matrix from metastats_matrix_creator.pl
######################################
###Heatmap colorTbl###
colorTbl <- rep(NA,0)
colorTbl[1] <- colors()[90]  #"darkorange"
colorTbl[2] <- colors()[553] #"red1"
colorTbl[3] <- colors()[47]  #"chartreuse"
colorTbl[4] <- colors()[12]  #"aquamarine4"
colorTbl[5] <- colors()[461] #"mediumblue"
colorTbl[6] <- colors()[526] #"lightsalmon2"
colorTbl[7] <- colors()[429] #"lightseagreen"
colorTbl[8] <- colors()[234] #"gray81"
colorTbl[9] <- colors()[652] #"yellow"
colorTbl[10] <- colors()[624] #"tan4"
colorTbl[11] <- colors()[550] #"purple1"
colorTbl[12] <- colors()[173] #"gray20"
colorTbl[13] <- colors()[103] #"darkseagreen1"
colorTbl[14] <- "gray"
colorTbl[15] <- "pink"
########################################
##Calculate the proportion for each Sample
RS=rowSums(tLCA);	    		##Calculate the total # of reads in each Sample (by row)
prop=(tLCA/RS)*100;			##Each row / total *100 = proportion of each Bacteria/Sample
################################
##For BacMM, selecting only the most abundant bacteria (sum of all %proportion for each bac > 10%)
bacSum=apply(prop, 2, sum)
propGrtrthan10 = prop[,bacSum >20]      ##Take only bacteria LCA w/ >20"%" combined across all patients
##Manually ordering the col by bacterial abundance
propOrd = propGrtrthan10[,order(colSums(propGrtrthan10),decreasing="T")]
propHier = hclust(dist(as.matrix(propOrd)), method="ward")  ## manually create a clustering of the data
##plot(hier)  ## show ~dendrogram of clustering
propCutHier = cutree(propHier,k=3)   ## grab only 5 clusterings
######################################## 
##Working with just counts##############
CS=colSums(tLCA);
counts=tLCA[,CS>5000];
countsOrd = counts[,order(colSums(counts),decreasing="T")];
countsHier = hclust(dist(as.matrix(countsOrd)), method="ward")  ## manually create a clustering of the data
##plot(hier)  ## show ~dendrogram of clustering
countsCutHier = cutree(countsHier,k=2)   ## grab only 5 clusterings
################################

####LOG transformation test#####
nrows <- nrow(propOrd);
ncols <- ncol(propOrd);
C <- array(0, dim=c(nrows,ncols));
for (i in 1:nrows){
  for (j in 1:ncols){
    C[i,j] = propOrd[i,j];
  }
}
Cnormed <- array(0,dim=c(nrows,ncols));
for (i in 1:nrows){
  for (j in 1:ncols){
    if (C[i,j] > 0){
      Cnormed[i,j] = log10(C[i,j]);
    } else {
      Cnormed[i,j] = 3;
    }
  }
}

for (i in 1:nrows){
  for (j in 1:ncols){
    if (Cnormed[i,j] == 3 ){
      Cnormed[i,j] = min(Cnormed) - 0.00001;
    }
  }
}

colnames(Cnormed) = colnames(propOrd);
rownames(Cnormed) = rownames(propOrd);
logOrd = Cnormed[,order(colSums(Cnormed),decreasing="T")];
logHier = hclust(dist(as.matrix(Cnormed)), method="ward");
logCutHier = cutree(logHier, k=3);

######LOG Transformation END###########
#######################################
##Setup a few variables for heatmap configuration
nrows = dim(prop)[1];                  ##nrows = # of rows in matrix 
ncols = dim(prop)[2];                  ##ncols = # of col. in matrix
sw1=.1/nrows;
sw2=.1/ncols;
reds = rev(rainbow(200))[1:10];        ##Setup colors to use in heatmap (from Skiff)
COLORS = c(rainbow(200)[32:200],reds); ##Setup colors to use in heatmap (from Skiff)
########################################
##Color Clustering for rows
propSideBars = cbind(colorTbl[propCutHier]);
rownames(propSideBars) = propCutHier;
countsSideBars = cbind(colorTbl[countsCutHier]);
rownames(countsSideBars) = countsCutHier;
logSideBars = cbind(colorTbl[logCutHier]);
rownames(logSideBars) = logCutHier;
#########################################  
#Set up color sidebar for patients by LGT/!LGT
srrlgt=read.table(file="/local/projects-t3/HLGT/TCGA/ksieber_dir/STAD/SRR_seqSourceSite.LGT.v2.txt", header = TRUE, na.strings = "NA");   ### This is the original used for the figures with LGT
SSSbar=cbind(colorTbl[srrlgt$SequencingSourceSite]);
rownames(SSSbar)=srrlgt$SRR;
LGTbar=cbind(colorTbl[srrlgt$SignificantLGT]);
rownames(LGTbar)=srrlgt$SRR;
propSIDEBAR2 = cbind(propSideBars,LGTbar,SSSbar);
countsSIDEBAR2 = cbind(countsSideBars,LGTbar,SSSbar);
logSIDEBAR2 = cbind(logSideBars,LGTbar,SSSbar);
##########################################
source("/home/ksieber/scripts/R/heatmap3.R")
##########################################
##Make COUNTS Heatmap

#pdf(file="/local/projects/HLGT/ksieber_dir/TCGA/Heatmap_LGT-vs-BacMM/LGT/STAD/4-new_lca_bacMM/STAD_combined.sam2lca.liberal_family_lca.counts.pdf");

heatmap.3(countsOrd,
         col=COLORS,                      # Use Chime coloring
         Colv=NA,                         # comment this out if you want to have columns clustered as well
         dendrogram="row",                # Draw the dendrogram only for the rows
         Rowv=as.dendrogram(countsHier),        # use our h. clustering
         RowSideColors=countsSIDEBAR2,          # this sets three color columns
         margins=c(10,5),                 # =c(bottom margin, right margin)
         ##sepwidth=c(sw1,sw2),             # Widthe between points
         NumRowSideColors=3,              # Number of side bars
         labRow=NA,                       # suppress row labels
         scale="none",                    # Normalize data
         trace="none",                    # Lines tracing each row and column data
         density.info="none",             # Density info as histogram of the distribution of values in the Key
         keysize = 0.8,                   # Size of the Key
         main="",                         # Title
         legend("bottomleft",
                legend=c("Clustering","LGT","SeqSite"),
                fill=c("darkorange","pink","lightseagreen"),
                border=FALSE,
                bty="y",
                y.intersp = 0.7,
                cex=0.7))

#dev.off();

#########################################
##Make PROPORTIONS Heatmap
#pdf(file="/local/projects/HLGT/ksieber_dir/TCGA/Heatmap_LGT-vs-BacMM/LGT/STAD/4-new_lca_bacMM/STAD_combined.sam2lca.liberal_family_lca.proportions.pdf");

heatmap.3(propOrd,
         col=COLORS,                      # Use Chime coloring
         Colv=NA,                         # comment this out if you want to have columns clustered as well
         dendrogram="row",                # Draw the dendrogram only for the rows
         Rowv=as.dendrogram(propHier),        # use our h. clustering
         RowSideColors=propSIDEBAR2,          # this sets three color columns
         margins=c(10,5),                 # =c(bottom margin, right margin)
         ##sepwidth=c(sw1,sw2),             # Widthe between points
         NumRowSideColors=3,              # Number of side bars
         labRow=NA,                       # suppress row labels
         scale="none",                    # Normalize data
         trace="none",                    # Lines tracing each row and column data
         density.info="none",             # Density info as histogram of the distribution of values in the Key
         keysize = 0.8,                   # Size of the Key
         main="",                         # Title
         legend("bottomleft",
                legend=c("Clustering","LGT","SeqSite"),
                fill=c("darkorange","pink","lightseagreen"),
                border=FALSE,
                bty="y",
                y.intersp = 0.7,
                cex=0.7))

#dev.off();
#########################################
#### LOg transformed
#pdf(file="/local/projects-t3/HLGT/TCGA/ksieber_dir/Heatmap_LGT-vs-BacMM/LGT/STAD/FOO.pdf");

heatmap.3(logOrd,
          col=COLORS,                      # Use Chime coloring
          Colv=NA,                         # comment this out if you want to have columns clustered as well
          dendrogram="row",                # Draw the dendrogram only for the rows
          Rowv=as.dendrogram(logHier),        # use our h. clustering
          RowSideColors=logSIDEBAR2,          # this sets three color columns
          margins=c(10,5),                 # =c(bottom margin, right margin)
          ##sepwidth=c(sw1,sw2),             # Widthe between points
          NumRowSideColors=3,              # Number of side bars
          labRow=NA,                       # suppress row labels
          scale="none",                    # Normalize data
          trace="none",                    # Lines tracing each row and column data
          density.info="none",             # Density info as histogram of the distribution of values in the Key
          keysize = 0.8,                   # Size of the Key
          main="",                         # Title
          legend("bottomleft",
                 legend=c("Clustering","LGT","SeqSite"),
                 fill=c("darkorange","pink","lightseagreen"),
                 border=FALSE,
                 bty="y",
                 y.intersp = 0.7,
                 cex=0.7))

#dev.off();

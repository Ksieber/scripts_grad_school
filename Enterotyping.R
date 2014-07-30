# Enterotyping Example: http://enterotype.embl.de/enterotypes.html

data=read.table("/local/projects-t3/HLGT/TCGA/ksieber_dir/MetaHIT_example/MetaHIT_SangerSamples.genus.txt", header=T, row.names=1, dec=".", sep="\t")
# data=read.table("/Users/Karsten/HLGT/MetaHIT_example/MetaHIT_SangerSamples.genus.txt", header=T, row.names=1, dec=".", sep="\t")
data=data[-1,]

dist_JSD <- function(inMatrix, pseudocount=0.000001, ...) {
	KLD <- function(x,y) { sum(x *log(x/y)) }
	JSD<- function(x,y) { sqrt(0.5 * KLD(x, (x+y)/2) + 0.5 * KLD(y, (x+y)/2)) }
	matrixColSize <- length( colnames( inMatrix ) )
	matrixRowSize <- length( rownames( inMatrix ) )
	colnames <- colnames( inMatrix )
	resultsMatrix <- matrix( 0, matrixColSize, matrixColSize )

	inMatrix = apply( inMatrix, 1:2, function(x) ifelse ( x==0, pseudocount, x ) )
	for ( i in 1:matrixColSize ) {
		for ( j in 1:matrixColSize ) {
			resultsMatrix[ i, j ] = JSD( as.vector( inMatrix[ , i ] ), as.vector( inMatrix[ , j ] ) )
		}
	}
	colnames -> colnames( resultsMatrix ) -> rownames( resultsMatrix )
	as.dist( resultsMatrix ) -> resultsMatrix
	attr( resultsMatrix, "method" ) <- "dist"
	return( resultsMatrix )
}

data_dist = dist_JSD( data )

pam.clustering=function(x,k) {  # x is a distance matrix and k the number of clusters
	require(cluster)
	cluster = as.vector( pam( as.dist(x), k, diss=TRUE ) $clustering)
	return( cluster) 
}
data_cluster = pam.clustering(data_dist, k=3 )

require(clusterSim)
nclusters = index.G1( t(data), data_cluster, d = data_dist, centrotypes = "medoids" )

for (k in 1:20 ) {
	if (k==1){
		nclusters[k]=NA
	} else {
		data_cluster_temp = pam.clustering( data_dist, k )
		nclusters[k] = index.G1(t(data), data_cluster_temp,  d = data_dist, centrotypes = "medoids")
	}
}

plot(nclusters, type="h", xlab="k clusters", ylab="CH index")
data_cluster = pam.clustering(data_dist, k=3 )

library(cluster)
obs_silhouette = mean( silhouette( data_cluster, data_dist )[,3] )

noise_removal <- function (dataframe, percent=0.01, top=NULL ){
	dataframe ->Matrix
	bigones <- rowSums( Matrix )*100/(sum( rowSums( Matrix ) ) ) > percent
	Matrix_1 <- Matrix[bigones,]
	print( percent )
	return( Matrix_1) 
}

library('ade4')
## Between Class Analysis: Original paper but not recommended w/ JSD
obs_pca = dudi.pca( data.frame(t(data)), scannf=F, nf=10)
obs_bet = bca(obs_pca, fac=as.factor(data_cluster), scannf=F, nf=k-1)
s.class(obs_bet$ls, fac=as.factor(data_cluster), grid=F) 

## PCoA: recommended w/ JSD
obs_pcoa=dudi.pco( data_dist, scannf=F, nf=3 )
s.class( obs_pcoa$li, fac=as.factor(data_cluster), grid=F ) ## With Grid & centroid connectors
s.class( obs_pcoa$li, fac=as.factor(data_cluster), grid=F, cell=0, cstar=0, col=c(3,2,4) ) ## NO Grid & centroid connectors
# Trying enterotyping strat on Isize
# data=read.table("/local/projects-t3/HLGT/TCGA/ksieber_dir/tmp/Isize_test.txt", header=T, row.names=1, sep="\t")

##################
# Load Libraries #
##################
require(boot)

##################
# Load functions #
##################
calc_JSD <- function(inMatrix, pseudocount=0.000001, ...) {
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
calc_JSD_boot_fxn = function (x_df, index) {
	tmp_df <- data.frame( x_df[,0], x_df[index,1] )
	return ( calc_JSD(tmp_df) ) 
}

##################
#### Run Data ####
##################
# Counts=read.table("/Users/ksieber/HLGT/tmp/Isize_test.txt", header=T, row.names=1, sep="\t")
# Counts=read.table("/Users/ksieber/HLGT/STAD/calc_integration_paper_data/Fig_LGT_reads/patient-a/opti_n/jsd_test/Counts_Isize_all_N.txt", header=T, row.names=1, sep="\t")
Table=read.table("/Users/ksieber/HLGT/STAD/calc_integration_paper_data/Fig_LGT_reads/patient-a/opti_n/jsd_working_but_not_right/Patient-a_minus-SRR203140.13068890_psort/Counts_Isize_all_N.txt", header=T, row.names=1, sep="\t")
Sample = data.frame(Table[,1],Table[,2])
data=prop.table(as.matrix(Table), margin=2)
data_dist = dist_JSD( data )
JSD = as.matrix(data_dist)[,1]

for (i in 2:(length( colnames( data ) ) ) ) {
	JSD_boot <- boot( as.matrix(data[,1],data[,i]), calc_JSD_boot_fxn, R=1000, parallel="multicore", ncpus=4 )
	JSD_ci <- boot.ci(JSD_boot, type="norm")	
}



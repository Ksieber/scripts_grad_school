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
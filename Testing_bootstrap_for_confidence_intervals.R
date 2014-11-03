## Testing bootstrap for confidence intervals
library(boot)
source("/Users/ksieber/Documents/vegdist/vegdist.R")
dyn.load("/Users/ksieber/Documents/vegdist/vegdist.so")

pop1 <- c(10,12,13,14)
pop2 <- c(91,92,93,94)

df <- data.frame(pop1,pop2)
## JSdist <- vegdist( df, method = "jensen-shannon", useShrinkage = TRUE )

calc.JSdist.fxn <- function(x_df,index){
	foo <- data.frame( x_df[,1], x_df[index,2] )
	return (vegdist(foo, method = "jensen-shannon", useShrinkage = TRUE ) )
}

JSdist_boot <- boot(df, calc.JSdist.fxn, R=100, parallel="multicore", ncpus=4)
JSdist_ci <- boot.ci(JSdist_boot,type="stud")
print(JSdist_ci)



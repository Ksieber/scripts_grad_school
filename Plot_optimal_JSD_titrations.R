# Plot_optimal_JSD_titrations.R
# Plot optimal JSD titrations with discontinuous and continuous x-axis
library(ggplot2);
##############################################################################################################
# Create JSD calibration figure 
JSD=read.table(file="/Users/ksieber/HLGT/STAD/calc_integration_paper_data/Fig_LGT_reads/optimal_jsd_possible/c5_downsampled/merged_reads-JSD.txt");
colnames(JSD)<-c("reads","JSD")
# Create whisker plot with discountinious x-axis
# pdf(file="/Users/ksieber/HLGT/STAD/calc_integration_paper_data/Fig_LGT_reads/optimal_jsd_possible/JSD_calibration_v2_dis.pdf")
JSD_Discontinuous_plot <- ggplot(JSD, aes(factor(reads),JSD)) + ylim(0,1)
JSD_Discontinuous_plot + geom_boxplot()
# dev.off();

# Create continous x-axis scatter plot
# pdf(file="/Users/ksieber/HLGT/STAD/calc_integration_paper_data/Fig_LGT_reads/optimal_jsd_possible/JSD_calibration_v2_cont.pdf")
JSD_Continuous_plot <- ggplot(JSD, aes(x=reads, y=JSD, group=reads)) + ylim(0,1) + scale_x_continuous(limits=c(0,100000))
JSD_Continuous_plot + geom_boxplot()
# dev.off();




##############################################################################################################

__END__

## This is code I tried to use, but it was ugly


# I got this from: http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
# Merge plots into 1 Figure (doesn't look good, wouldn't recommend)
# pdf(file="/Users/ksieber/HLGT/STAD/calc_integration_paper_data/Fig_LGT_reads/optimal_jsd_possible/JSD_calibration_v2_dis-cont.pdf")
# multiplot(JSD_Discontinuous_plot, JSD_Continuous_plot, cols=2)
# dev.off();

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
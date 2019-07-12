#'
#'@title Plot a growth transition matrix as an image.
#'
#'@description Function to plot a growth transition matrix as an image.
#'
#'@param prGr - growth transition matrix
#'@param log - flag to plot ln-scale values
#'
#'@return a ggplot2 object
#'
#'@details None.
#'
#'@importFrom RColorBrewer brewer.pal
#'@importFrom reshape2 melt
#'
#'@import ggplot2
#'
#'@export
#'
plotSizeTransitionMatrix<-function(prGr,log=FALSE){

    colors <- brewer.pal(11,"Spectral")
    
    mPrGr<-melt(prGr);
    if (log) {mPrGr$value=log(mPrGr$value);}
    
    p<-ggplot(mPrGr,aes(from,to,z=value,fill=value)) +
        labs(x="pre-molt size (mm CW)",y="post-molt size (mm CW)",fill="pr(to|from)") + 
        scale_x_continuous(expand = c(0, 0)) + 
        scale_y_continuous(expand = c(0, 0)) +
        scale_fill_gradient2(low=colors[10],mid=colors[6],high=colors[2],
                             midpoint=0.5,limits=c(0,1)) +
        geom_raster() + 
        theme(plot.background = element_rect(fill = "grey90"), 
              legend.background = element_rect(fill = "grey90")) + 
        geom_contour() + 
        geom_line(data=melt(as.numeric(rownames(prGr))),mapping=aes(x=value,y=value))
    print(p)
}
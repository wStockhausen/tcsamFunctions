#'
#'@title Make bubble plots of a matrix or list of matrices
#'
#'@description Function to plot a matrix or list of matrices using a bubble plot.
#'
#' @param lst : list or matrix
#' @param x : column name for x axis
#' @param x.lbl : label for x axis
#' @param y : column name for y axis
#' @param y.lbl : label for y axis
#' @param z : column name for bubble sizes
#' @param z.lbl : label for sizes
#' @param fill : column name for bubble fill colours
#' @param fill.lbl : label for fills
#' @param scaleBy : value to scale bubbles by
#' @param normalize : flag (T/F) to normalize sums to 1
#' @param  xlims : x-axis limits (or NULL)
#' @param  ylims : y-axis limits (or NULL)
#' @param colors : color palette to use (default is \code{RColorBrewer::brewer.pal(11,"Spectral")})
#' @param alpha : transparency level
#' @param showPlot : flag to print plot
#' 
#'@importFrom reshape2 melt
#'@importFrom RColorBrewer brewer.pal
#'@import ggplot2
#'
#'@export
#'

# require('reshape2');
# require('ggplot2');
# require('RColorBrewer')
    
makeBubblePlots<-function(lst,
                         x='year',x.lbl=x,
                         y='size',y.lbl=y,
                         z='value',z.lbl=z,
                         fill='sex',fill.lbl=fill,
                         scaleBy=1,
                         normalize=FALSE,
                         xlims=NULL,
                         ylims=NULL,
                         colors=RColorBrewer::brewer.pal(11,"Spectral"),
                         alpha=0.5,
                         showPlot=TRUE){
    if (is.list(lst)){
        if (normalize) {
            for (el in names(lst)){
                m<-lst[[el]];
                cat('\n\nBefore normalization by age\n')
#                print(m);
                for (cn in colnames(m)) {
                    if (sum(m[,cn])>0) {m[,cn]<-m[,cn]/sum(m[,cn]);}
                }
                cat('\n\nAfter normalization by age\n')
#                print(m);
                lst[[el]]<-m;
            }
        }
        
        dfr<-melt(lst);
        nms<-names(dfr);
#        cat('names = ',nms,'\n')
        nms[which(names(dfr)=='L1')]<-fill;
#        cat('names = ',nms,'\n')
        names(dfr)<-nms;
#        print(dfr)
    }
    if (is.matrix(lst)){
        m<-lst;
        if (normalize) {
            for (cn in colnames(m)) {
                sm<-sum(m[,cn]);
                if (sm>0) {m[,cn]<-m[,cn]/sm;}
            }
        }        
        dfr<-melt(m);
        dfr[[fill]]<-fill;
    }
    
    mxz<-max(dfr[[z]],na.rm=TRUE);
    
#    cat(names(dfr),'\n')
#    cat('max z = ',mxz,'\n')
    
#    dfr[[z]]<-dfr[[z]]*(scaleBy/mxz)
    
    p<-ggplot(data=dfr,guide=FALSE) +
        labs(x=x.lbl,y=y.lbl,fill=fill.lbl) + 
        geom_point(aes_string(x=x,y=y,size=z,fill=fill),shape=21,alpha=alpha,na.rm=TRUE) + 
        scale_size_area(max_size=scaleBy*mxz) + 
        scale_x_continuous(expand = 0.05*c(1, 1),limits=xlims) + 
        scale_y_continuous(expand = 0.05*c(1, 1),limits=ylims) +
        theme(plot.background = element_rect(fill = "grey90"))

    if (showPlot) print(p);

    return(p);
}

# p1<-makeBubblePlots(sg,x='age',y='size',z='value',fill='maturity',scaleBy=20,normalize=FALSE)
# p2<-makeBubblePlots(sg,x='age',y='size',z='value',fill='maturity',scaleBy=20,normalize=TRUE)
# p3<-makeBubblePlots(sg$immature,x='age',y='size',z='value',fill='immature',scaleBy=20,normalize=FALSE)
# p4<-makeBubblePlots(sg$immature,x='age',y='size',z='value',fill='immature',scaleBy=20,normalize=TRUE,ylims=c(0,185))




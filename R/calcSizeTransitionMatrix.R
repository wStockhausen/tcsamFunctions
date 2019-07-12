#' 
#' @title Calculate a size transition matrix.
#' 
#' @description Function to calculate a size transition matrix.
#'
#' @param coeffs - list including function 'type' and parameter values
#' @param sizes - vector of sizes at which to calculate the sie transition probabilities
#' @param showPlot - flag (T/F) to plot the matrix
#' @param log - flag to use ln-scale on plot
#' @param colors - paleete to use to color the plot
#' 
#' @return 
#' 
#' @details If 
#' * 'type'='LinearGrowthIncrement', calls \code{calcSizeTransitionMatrix.LinearGrowthIncrement}.
#'
#' @export
#'
calcSizeTransitionMatrix<-function(coeffs=list(type='LinearGrowthIncrement',
                                               a=0.70000,
                                               b=0.883118,
                                               beta=0.75,
                                               maxDZ=50),
                                   sizes=seq(from=25,to=185,by=1),
                                   showPlot=FALSE,
                                   log=FALSE,
                                   colors=rainbow(1000)){
    if (coeffs$type=='LinearGrowthIncrement'){
        return(invisible(calcSizeTransitionMatrix.LinearGrowthIncrement(coeffs=coeffs,sizes=sizes,showPlot=showPlot,log=log,colors=colors)))
    }
    return(NULL);
}

#prGr<-calcSizeTransitionMatrix(colors=cm.colors(20),showPlot=TRUE);


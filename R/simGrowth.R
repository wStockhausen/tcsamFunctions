#'
#'@title Simulate individual growth trajectories based on a size transition matrix
#'
#'@description Function to simulate growth trajectories based on a size transition matrix.
#'
#'@param ny - number of years to simulate
#'@param zAtR - size at recruitment
#'@param M - annual rate of natural mortality
#'@param prGr - size transition matrix
#'@param grCoeffs - growth coefficients (if prGr is not given)
#'@param prM2M - probability of molt-to-maturity 
#'@param m2mCoeffs - molt-to-maturity coefficients (if prM2M is not given)
#'@param sizes - vector of size bins
#'@param showPlot - flag to plot simulation
#'@param normalize - flag to normalize annual size comps to 1
#'
#'@return list with matrices nAtZI and nAtZM representing immature and mature crab abundance 
#'as a function of elapsed time.
#'
#'@details Essentially calculates a cohort progression
#' 
#'@export
#'
simGrowth<-function(ny=21,zAtR=1,M=0.23,
                    prGr=NULL,grCoeffs=NULL,
                    prM2M=NULL,m2mCoeffs=NULL,
                    sizes=seq(from=25,to=185,by=1),
                    showPlot=TRUE,
                    normalize=FALSE){
    if (is.null(prGr)){
        if (!is.null(grCoeffs)){
            prGr<-calcSizeTransitionMatrix(coeffs=grCoeffs,sizes=sizes,showPlot=TRUE);
        } else {
            prGr<-calcSizeTransitionMatrix(sizes=sizes,showPlot=TRUE);
        }
    }
    if (is.null(prM2M)){
        if (!is.null(m2mCoeffs)){
            prM2M<-calcPrMoltToMaturity(coeffs=m2mCoeffs,sizes=sizes);
        } else {
            prM2M<-calcPrMoltToMaturity(sizes=sizes);
        }
    }
    
    nr<-nrow(prGr);
    print(ny);
    print(nr);
    nAtZI<-matrix(data=0,ncol=ny,nrow=nr,dimnames=list(size=sizes,age=0:(ny-1)));
    nAtZM<-matrix(data=0,ncol=ny,nrow=nr,dimnames=list(size=sizes,age=0:(ny-1)));
    if (length(zAtR)==nr) {nAtZI[,1]<-zAtR;} else {nAtZI[1,1]<-zAtR;}
    
    for (y in 2:ny){
        nAtZI[,y]<-prGr %*% nAtZI[,y-1];                 #growth without regard to terminal molt
        nAtZM[,y]<-(nAtZM[,y-1]+prM2M*nAtZI[,y])*exp(-M);#terminally-molting crab + mortality
        nAtZI[,y]<-(1-prM2M)*nAtZI[,y]*exp(-M);          #non-terminally-molting crab + mortality
        #cat('y = ',y,'\n')
        #cat('\t',t(nAtZI[,y]));
    }
    
    lst<-list(immature=nAtZI,mature=nAtZM);
    
    if (showPlot){
        p<-makeBubblePlots(lst,x='age',y='size',z='value',fill='maturity',scaleBy=20,normalize=normalize)
        print(p)
    }
    return(lst);
}

#sg<-simGrowth(M=0.23,sizes=seq(from=25,to=185,by=5),showPlot=TRUE,normalize=FALSE);

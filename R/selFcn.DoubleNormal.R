#'
#'@title Calculate a double normal selectivity function
#'
#'@description Function to calculate (and optionally plot the components of) a double normal selectivity function.
#'
#'@param x : age/size bins at which to compute values
#'@param params : vector of parameters (see Details)
#'@param showPlot : flag (T/F) to plot components of function
#'@param test : flag (T/F) to output list with detailed info (T) rather than just vector of sel function values (F)
#'
#'@details The \code{params} vector elements are:
#' \itemize{
#'  \item{1: location of ascending peak}
#'  \item{2: logistic scale width of plateau}
#'  \item{3: ascending slope  (NOT log scale; i.e. = exp(p\[3\]) if p\[3\] from SS)}
#'  \item{3: descending slope (NOT log scale; i.e. = exp(p\[4\]) if p\[4\] from SS)}
#'  \item{5: logistic scale value in initial bin; i.e. initial = 1./(1.+exp(-params\[5\]));}
#'  \item{6: logistic scale value in final bin;   i.e. final   = 1./(1.+exp(-params\[6\]));}
#' }
#'
#'@import graphics
#'
#'@export
#'
selFcn.DoubleNormal<-function(x,
                              params,
                              showPlot=FALSE,
                              test=FALSE) {
    #in ADMB code,
    #   x was len_bins_m (midpoints of size/age bins)
    #   binwidth was binwidth(nlength/2) (now assuming constant binwidth)
    #   minz was len_bin_m[1] (now assumed to be min(x)+binwidth/2)
    #   maxz was len_bin_m(nlength) (now assumed to be max(x)+binwidth/2)
    binwidth<-x[2]-x[1];
    minz<-min(x);
    maxz<-max(x);

    sel<-0*x;#selectivity curve
    asc<-0*x;#ascending limb
    dsc<-0*x;#descending limb
    asc_scl<-0*x;#scaled ascending limb
    dsc_scl<-0*x;#scaled descending limb
    join1<-0*x;#join1
    join2<-0*x;#join2

    peak1    <-params[1];
    peak2    <-params[1]+binwidth+ (0.99*maxz-params[1]-binwidth)/(1.+exp(-params[2]));
    upselex  <-params[3]; #was exp(params(3))
    downselex<-params[4]; #was exp(params(4))
    if(params[5]>-999) {
        point1<-1./(1.+exp(-params[5]));
        t1min<-exp(-((minz-params[1])^2)/upselex); # unscaled fcn at first bin
    }
    if(params[6]>-999) {
        point2<-1./(1.+exp(-params[6]));
        t2min<-exp(-((maxz-peak2)^2)/downselex);  #unscaled  fcn at last bin
    }
    for (j in 1:length(x)) {
        t1<-x[j]-peak1;
        t2<-x[j]-peak2;
        join1[j]<-1./(1.+exp(-(20./(1.+abs(t1)))*t1));
        join2[j]<-1./(1.+exp(-(20./(1.+abs(t2)))*t2));
        asc[j]<-exp(-(t1^2)/upselex);
        asc_scl[j]<-asc[j];
        if(params[5]>-999) {
            asc_scl[j]<-point1+(1.-point1)*(exp(-(t1^2)/upselex)-t1min)/(1.-t1min);
        }
        dsc[j]<-exp(-(t2^2)/downselex);
        dsc_scl[j]<-dsc[j];
        if(params[6]>-999) {
            dsc_scl[j]<-1+(point2-1.)*(exp(-(t2^2)/downselex)-1.)/(t2min-1.);
        }
        sel[j]<-asc_scl[j]*(1.-join1[j])+join1[j]*(1.-join2[j]+dsc_scl[j]*join2[j]);
    }

    if (showPlot){
        plot(x,sel,col="red");
        points(x,asc,col="blue",pch=24);
        lines(x,asc_scl,col="blue");
        lines(x,join1,col="violet");
        points(x,dsc,col="green",pch=25);
        lines(x,dsc_scl,col="green");
        lines(x,join2,col="yellow");
        points(x,sel,col="red");
        lines(x,sel,col="red");
    }
    
    if (!test) return(sel);
    
    #--output tibble
    tbl=tibble::tibble(x=x,sel=sel,
                       asc=asc,asc_scl=asc_scl,join1=join1,
                       dsc=dsc,dsc_scl=dsc_scl,join2=join2);
    return(list(params=params,
                binwidth=binwidth,minz=minz,maxz=maxz,
                peak1=peak1,peak2=peak2,
                upselex=upselex,downselex=downselex,
                point1=point1,t1min=t1min,
                point2=point2,t2min=t2min,
                tbl=tbl));
}
# if (FALSE){
#     #--test the above using Mike Byerly's parameters
#     x = 20:80;
#     params = c(54.0494, 15, exp(4.00043), exp(-15), -15, 15)
#     lst = selFcn.DoubleNormal(x,params,showPlot=TRUE,test=TRUE);
#     tbl = lst$tbl;
#     require(ggplot2);
#     #--unscaled and scaled ascending components and join
#     ggplot(tbl,aes(x=x))+
#         geom_line(aes(y=asc),    colour="red")   + geom_point(aes(y=asc),colour="red") +
#         geom_line(aes(y=asc_scl),colour="blue")  + geom_point(aes(y=asc_scl),colour="blue") +
#         geom_line(aes(y=join1),  colour="green") + geom_point(aes(y=join1),colour="green");
#     #--unscaled and scaled descending components and join
#     ggplot(tbl,aes(x=x))+
#         geom_line(aes(y=dsc),    colour="red")   + geom_point(aes(y=dsc),colour="red") +
#         geom_line(aes(y=dsc_scl),colour="blue")  + geom_point(aes(y=dsc_scl),colour="blue") +
#         geom_line(aes(y=join2),  colour="green") + geom_point(aes(y=join2),colour="green");
#     #--ascending, descending components and final selectivity curve
#     ggplot(tbl,aes(x=x))+
#         geom_line(aes(y=asc_scl),colour="red")   + geom_point(aes(y=asc_scl),colour="red") +
#         geom_line(aes(y=dsc_scl),colour="blue")  + geom_point(aes(y=dsc_scl),colour="blue") +
#         geom_line(aes(y=sel),    colour="green") + geom_point(aes(y=sel),    colour="green");
#     #--resulting selectivity curve
#     ggplot(tbl,aes(x=x))+
#         geom_line(aes(y=sel),    colour="green") + geom_point(aes(y=sel),    colour="green");
# }

#'
#'@title Calculate a double normal selectivity function with "4a" parameterization
#'
#'@description Function to calculate (and optionally plot the components of) a double normal selectivity function 
#'with "4a" parameterization.
#'
#'@param zBs : vector of sizes at which to compute values
#'@param params : vector of parameters
#'@param fsZ : max possible size on the "shelf"
#'@param showPlot : flag (T/F) to plot components of function
#'
#'@details Should replicate TCSAM02 dblnormal4a selectivity function.
#'
#'\itemize{
#' \item{params\[1\]: size at which ascending limb hits 1}
#' \item{params\[2\]: width of ascending limb}
#' \item{params\[3\]: scaled size at which descending limb departs from 1}
#' \item{params\[4\]: width of descending limb}
#'}
#'
#'@import ggplot2
#'@import reshape2
#'@import tibble
#'
#'@export
#'
selFcn.DoubleNormal4a<-function(zBs,
                                params,
                                fsZ,
                                showPlot=TRUE) {
    slp = 5.0;
    ascMnZ = params[1];#--size at which ascending limb hits 1
    ascWdZ = params[2];#--width of ascending limb
    sclInc = params[3];#--scaled size at which descending limb departs from 1
    dscWdZ = params[4];#--width of descending limb
    dscMnZ = (fsZ-ascMnZ)*sclInc + ascMnZ;
    ascN  = mfexp(-0.5*square((zBs-ascMnZ)/ascWdZ));
    ascJ  = 1.0/(1.0+mfexp(slp*(zBs-(ascMnZ))));
    dscN = mfexp(-0.5*square((zBs-dscMnZ)/dscWdZ));
    dscJ = 1.0/(1.0+mfexp(-slp*(zBs-(dscMnZ))));
    s = elem_prod(elem_prod(ascJ,ascN)+(1.0-ascJ), elem_prod(dscJ,dscN)+(1.0-dscJ));
    
    if (showPlot){
        dfr  = tibble::tibble(zBs=zBs,ascN=ascN,ascJ=ascJ,dscN=dscN,dscJ=dscJ,s=s);
        mdfr = reshape2::melt(dfr,id.vars="zBs",variable.name="variable",value.name="value");
        p = ggplot2::ggplot(mdfr,mapping=ggplot2::aes(x=zBs,y=value,colour=variable));
        p = p + ggplot2::geom_line();
        p = p + ggplot2::geom_vline(xintercept=ascMnZ,linetype=1,colour="grey",size=2,alpha=0.2);
        p = p + ggplot2::geom_vline(xintercept=dscMnZ,linetype=1,colour="grey",size=2,alpha=0.2);
        print(p);
    }
    return(s);
}
# require(wtsADMB);
# zBs = seq(27.5,182.5,5);
# params = c(140,23.6046517434,0.237607333911,29.9999999390);
# selFcn.DoubleNormal4a(zBs,params,182);

#'
#'@title Calculate a ascending normal selectivity function with "3" parameterization
#'
#'@description Function to calculate (and optionally plot the components of) an ascending normal selectivity function 
#'with "3" parameterization.
#'
#'@param zBs : vector of sizes at which to compute values
#'@param params : vector of parameters
#'@param fsZ    : max possible size at which selectivity reaches 1
#'@param showPlot : flag (T/F) to plot components of function
#'
#'@details Should replicate TCSAM02 ascnormal3 selectivity function.
#'
#'params are:\cr
#'  1: delta size from max possible (params\[4\])\cr
#'  2: ln-scale selectivity at size givens by params\[3\]\cr
#'  3: size at selectivity is given by params\[2\]\cr
#'  4: max possible size at which peak is reached
#'
#'@import ggplot2
#'@import reshape2
#'@import tibble
#'
#'@export
#'
selFcn.AscNormal3<-function(zBs,
                            params,
                            fsZ,
                            showPlot=TRUE) {
    slp = 5.0;
    ascZ1     = fsZ-params[1];#--size at which ascending limb hits 1
    ascSref   = params[2];    #--selectivity at ascZref
    ascZref   = params[3];    #--size at which selectivity reaches ascSref
    ascN      = mfexp(log(ascSref)*square((zBs-ascZ1)/(ascZref-ascZ1)));
    ascJ      = 1.0/(1.0+mfexp(slp*(zBs-(ascZ1))));
    s = elem_prod(ascJ,ascN)+(1.0-ascJ);
    
    if (showPlot){
        dfr  = tibble::tibble(zBs=zBs,ascN=ascN,ascJ=ascJ,s=s);
        mdfr = reshape2::melt(dfr,id.vars="zBs",variable.name="variable",value.name="value");
        p = ggplot2::ggplot(mdfr,mapping=ggplot2::aes(x=zBs,y=value,colour=variable));
        p = p + ggplot2::geom_line();
        p = p + ggplot2::geom_vline(xintercept=ascZ1,linetype=2,colour="grey",size=2,alpha=0.2);
        p = p + ggplot2::geom_vline(xintercept=ascZref,linetype=2,colour="grey",size=2,alpha=0.2);
        p = p + ggplot2::geom_hline(yintercept=ascSref,linetype=2,colour="grey",size=2,alpha=0.2);
        print(p);
    }
    return(s);
}
# zBs = seq(27.5,182.5,5);
# params = c(80.4968168760,exp(-1.67535227638),32);
# selFcn.AscNormal3(zBs,params,132);

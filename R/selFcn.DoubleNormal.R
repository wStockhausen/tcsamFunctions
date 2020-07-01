#'
#'@title Calculate a double normal selectivity function
#'
#'@description Function to calculate (and optionally plot the components of) a double normal selectivity function.
#'
#'@param x - age/size bins at which to compute values
#'@param params - list or vector of parameters
#'@param showPlot - flag (T/F) to plot components of function
#'
#'@details TBD!!
#'
#'@import graphics
#'
#'@export
#'
selFcn.DoubleNormal<-function(x,
                              params,
                              showPlot=FALSE) {
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

    peak2<-params[1]+binwidth+ (0.99*maxz-params[1]-binwidth)/(1.+exp(-params[2]));
    upselex<-params[3];   #was exp(params(3))
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
        t1<-x[j]-params[1];
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

    return(sel);
}
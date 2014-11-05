#'
#'@title Calculate a size transition matrix.
#'
#'@description Function to calculate a size transition matrix
#'
#'@param a - 
#'@param b - 
#'@param beta - gamma distribution scale parameter
#'@param maxDZ - max growth increment allowed (truncates potential growth)
#'@param coeffs - list with coefficients a, b, beta, maxDZ 
#'@param bins - size bins to include in matrix
#'@param showPlot - flag (T/F) to plot matrix
#'@param log - flag to plot ln-scale matrix (in conjunction w/ showPlot)
#'
#'@export
#'
calcSizeTransitionMatrix.LinearGrowthIncrement<-function(a=0.70000,
                                                         b=0.883118,
                                                         beta=0.75,
                                                         maxDZ=50,
                                                         coeffs=NULL,
                                                         sizes=seq(from=27.5,to=182.5,by=5),
                                                         showPlot=TRUE,
                                                         log=FALSE,
                                                         colors=rainbow(1000)){
    if (!is.null(coeffs)){
        if (!is.null(coeffs$a))     a<-coeffs$a;
        if (!is.null(coeffs$b))     b<-coeffs$b;
        if (!is.null(coeffs$beta))  beta<-coeffs$beta;
        if (!is.null(coeffs$maxDZ)) maxDZ<-coeffs$maxDZ;
        cat('a = ',a,'\nb = ',b,'\nbeta = ',beta,'\nmaxDZ = ',maxDZ,'\n')
    }
#         //compute growth transition matrix for this pc
#         prGr_zz.initialize();
#         dvar_vector mnZ = mfexp(grA)*pow(zBs,grB);//mean size after growth from zBs
#         dvar_vector alZ = (mnZ-zBs)/grBeta;//scaled mean growth increment from zBs
#         for (int z=1;z<nZBs;z++){//pre-molt growth bin
#             dvar_vector dZs =  zBs(z,nZBs) - zBs(z);//realized growth increments (note non-neg. growth only)
#             if (debug) cout<<"dZs: "<<dZs.indexmin()<<":"<<dZs.indexmax()<<endl;
#             dvar_vector prs = elem_prod(pow(dZs,alZ(z)-1.0),mfexp(-dZs/grBeta)); //pr(dZ|z)
#             if (debug) cout<<"prs: "<<prs.indexmin()<<":"<<prs.indexmax()<<endl;
#             if (prs.size()>10) prs(z+10,nZBs) = 0.0;//limit growth range TODO: this assumes bin size is 5 mm
#             if (debug) cout<<prs<<endl;
#             prs = prs/sum(prs);//normalize to sum to 1
#             if (debug) cout<<prs<<endl;
#             prGr_zz(z)(z,nZBs) = prs;
#         }
#         prGr_zz(nZBs,nZBs) = 1.0; //no growth from max size
#         prGr_czz(pc) = trans(prGr_zz);//transpose so rows are post-molt (i.e., "to") z's so n+ = prGr_zz*n
    zBs<-sizes;
    nZBs<-length(zBs);
    dZ<-zBs[2]-zBs[1]; #size bin increment
    mnZ<-exp(a)*zBs^b;
    a1Z<-(mnZ-zBs)/beta;
    prGr<-matrix(data=0,nrow=nZBs,ncol=nZBs)
    for (z in 1:(nZBs-1)){
        dZs<-zBs[z:nZBs]-zBs[z];#realized growth increment (non-negative only)
        idx<-dZs<=maxDZ;#limit max growth to maxDZ
        prs<-(dZs^(a1Z[z]-1.0))*exp(-dZs/beta)*idx; #unscaled, truncated pr(dZ|z)
        prs<-prs/sum(prs);#normalize to sum to 1
        prGr[z,z:nZBs]<-prs;
    }
    prGr[nZBs,nZBs]<-1;
    dimnames(prGr)<-list(from=as.character(zBs),to=as.character(zBs));

    #transpose matrix so "from" corresponds to columns, not rows
    prGr<-t(prGr);
    
    if (showPlot){
        plotSizeTransitionMatrix(prGr,log=log)
    }
    
    return(prGr);
}

#prGr<-calcSizeTransitionMatrix.LinearGrowthIncrement(colors=cm.colors(20),beta=3);

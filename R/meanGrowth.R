################################################################################
#' 
#' @title Retrieve a list of growth parameters
#' 
#' @description Function to retrieve a list of growth parameters.
#' 
#' @return list with elements 'assessment', 'tmb', and 'DonaldsonEtAl1981', each of which 
#' is a tibble with sex-specific growth parameters from the assessment model, the TMB
#' analysis, and the Donaldson et al., 1981 study.
#' 
#' @details Port of TCSAM02 size transition calculations from ADMB to R
#' 
#' @importFrom tibble tibble
#' 
#' @export
#' 
getGrowthParams<-function(){
    #--from assesssment model
    assessment=rbind(tibble::tibble(sex="male",  pGrA=32.4293064532,pGrB=164.709991529,zGrA=25,zGrB=125),
                     tibble::tibble(sex="female",pGrA=33.5510328222,pGrB=114.997323294,zGrA=25,zGrB=100));
    #--from TMB model
    tmb=rbind(tibble::tibble(sex="male",  pGrA=32.2467787,pGrB=100.0611878,zGrA=25,zGrB=80),
              tibble::tibble(sex="female",pGrA=32.6168227,pGrB=93.3340492,zGrA=25,zGrB=80));
    #--post_molt = c0 + c1*pre_molt
    DonaldsonEtAl1981=
        rbind(tibble::tibble(sex="male",  c0=0.02, c1=1.32, instars=2:10,  range=c(10,42)),
              tibble::tibble(sex="male",  c0=3.98, c1=1.19, instars=11:13, range=c(42,84)),
              tibble::tibble(sex="male",  c0=14.75,c1=1.07, instars=14:18, range=c(84,130)),
              tibble::tibble(sex="female",c0=0.02, c1=1.32, instars=2:10,  range=c(10,42)),
              tibble::tibble(sex="female",c0=3.98, c1=1.19, instars=11:12, range=c(42,80)),
              tibble::tibble(sex="female",c0=17.59,c1=0.96, instars=13,    range=c(70,90)));
    return(list(assessment=assessment,tmb=tmb,DonaldsonEtAl1981=DonaldsonEtAl1981));
}

# getInstarSizes.DonaldsonEtAl1981<-function(){
#     
# }

################################################################################
#' 
#' @title Calculate mean growth as in TCSAM02
#' 
#' @description Function to calculate mean growth as in TCSAM02.
#' 
#' @param zBs - size bins
#' @param params - list of growth parameters
#' 
#' @return vector
#' 
#' @details Port of TCSAM02 size transition calculations from ADMB to R
#' 
#' @importFrom tibble tibble
#' 
#' @export
#' 
calcMeanGrowth.tcsam02<-function(zBs,
                                 params){
    mnZs = params$grA*exp(log(params$grB/params$grA)/log(params$zGrB/params$zGrA)*log(zBs/params$zGrA));
    return(tibble::tibble(sex=params$sex,premolt=zBs,postmolt=mnZs));
}

################################################################################
#' 
#' @title Calculate mean growth as in Donaldson et al., 1981
#' 
#' @description Function to calculate mean growth as in Donaldson et al., 1981.
#' 
#' @param zBs - size bins
#' 
#' @return vector
#' 
#' @details See Donaldson et al., 1981.
#' 
#' @importFrom tibble tibble
#' 
#' @export
#' 
calcMeanGrowth.DonaldsonEtAl1981<-function(zBs){
    params<-getGrowthParams()$DonaldsonEtAl1981;
    dfr<-NULL;
    for (r in 1:nrow(params)){
        ps = params[rw,];
        idx = (ps$range[1]<=zBs)&(zBs<=ps$range[2]);
        if (sum(idx)>0){
            mnZs = ps$c0 +  ps$C1 * zBs[idx];
            dfr<-rbind(dfr,tibble::tibble(sex=ps$sex,premolt=zBs[idx],postmolt=mnZs));
        }
    }
    return(dfr);
}
################################################################################
#' 
#' @title Calculate the size transition matrix as in TCSAM02
#' 
#' @description Function to calculate the size transition matrix as in TCSAM02.
#' 
#' @param zBs - size bins
#' @param delZ - size bin width (mm CW)
#' @param params - list of growth parameters
#' 
#' @return matrix
#' 
#' @details Port of TCSAM02 size transition calculations from ADMB to R
#' 
#' @export
#' 
calcSizeTransitionMtrix.tcsam02<-function(zBs,
                                          delZ,
                                          params){
    mnZs = calcMeanGrowth.tcsam02(zBs,params)
    mnIs = mnZs - zBs;          #mean molt increments
    invBeta = 1.0/params$grBeta;#inverse scale (rate) for gamma density function
    alIs = mnIs*invBeta;        #gamma density alpha (location) parameters
    mnpIs = mnZs - (zBs-delZ/2);#mean molt increment (adjusted to start of size bin)
    alpIs = mnpIs*invBeta;      #gamma density alpha (location) parameters
     
    nZBs = length(zBs);
    zCs  = c(zBs-delZ/2,zBs[nZBs]+delZ/2);
    prGr_zz = matrix(0,nrow=nZBs,ncol=nZBs);
    for (z in 1:(nZBs-1)){
        cprs = vector(0,length=(nZBs-z+1));
        prs = vector(0,length=(nZBs-z+1));
        sclIs = (zCs[(z+1):(nZBs+1)]-zBs[z])*invBeta;#scaled increments at size bin cut points
        cprs[z] = pgamma(sclIs(z+1),alIs(z));
        prs[z]  = cprs[z];
        for (zp in (z+1):nZBs){
            cprs[zp] = cumd_gamma(sclIs(zp+1),alIs(z));
            prs[zp]  = cprs[zp]-cprs[zp-1];#--cumulative pr from zCs(zp) to zCs(zp+1)
        }
        prs[nZBs] = prs[nZBs] + 1.0 - cprs(nZBs);         #--treat final size bin as accumulator
        if (length(prs)>maxZBEx) prs[(z+maxZBEx):nZBs] = 0.0;#--limit growth range
        prs = prs/sum(prs);#--normalize to sum to 1
        if (debug) cat(prs,"\n");
        prGr_zz[z,z:nZBs] = prs;
    }#--z
    prGr_zz[nZBs,nZBs] = 1;#no growth from accumulator bin
    return(prGr_zz);
}


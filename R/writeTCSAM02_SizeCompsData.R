#' 
#' @title Write size composition data to a connection in TCSAM02 input format
#' 
#' @description Function to write size composition data to a connection in TCSAM02 input format.
#' 
#' @param con : connection to use (default=stdout)
#' @param dfrZCs : size comps dataframe to write to connection
#' @param dfrSSs : sample size dataframe corresponding to dfrZCs
#' @param cutpts : vector of cutpoints that were used to create the size comps
#' @param tail_compression : two-element vector of compresssion factors
#' @param optFit : objective function fitting option (e.g., "BY_X","BY_XM")
#' @param likeType : likelihood type ("MULTINOMIAL" or "DIRICHLET-MULTINOMIAL")
#' @param likeWgt : likelihood weight (multiplier)
#' @param unitsIn : input catch units ("ONES" or "MILLIONS")
#' @param unitsOut : output catch units ("MILLIONS" or "ONES)
#' 
#' @return Invisibly returns the connection \code{con} to allow piping.
#' 
#' @details The input values should be in ONES for abundance or KG for biomass data.
#' 
#' The input dataframes should have the following column names
#' \itemize{
#'   \item{year}
#'   \item{sex}
#'   \item{maturity}
#'   \item{shell condition}
#'   \item{size (dfrZCs only)}
#'   \item{value (dfrZCs) or ss (dfrSSs)}
#' }
#' 
#' @export
#' 
writeTCSAM02_SizeCompsData<-function(con=stdout(),
                                     dfrZCs=NULL,
                                     dfrSSs=NULL,
                                     cutpts=NULL,
                                     tail_compression=c(0,0),
                                     optFit=c("BY_X","BY_XM"),
                                     likeType=c("MULTINOMIAL","DIRICHLET-MULTINOMIAL"),
                                     likeWgt = 1,
                                     unitsIn=c("ONES","MILLIONS"),
                                     unitsOut=c("MILLIONS","ONES")){
    #--use only first element of these
    optFit   = toupper(optFit[1]);
    likeType = toupper(likeType[1]);
    likeWgt  = likeWgt[1];
    unitsIn  = toupper(unitsIn[1]);
    unitsOut = toupper(unitsOut[1]);
    
    scl = getScaleForAbundance(unitsIn,unitsOut);
    
    #--column names in dataframe
    yr  = "year";
    sx  = "sex";
    mt  = "maturity";
    sc  = "shell condition";
    sz  = "size";
    vlc = "value";#--value column in dfrZCs
    ss  = "ss";   #--value column in dfrSSs
    
    tmpZCs = dfrZCs;
    tmpSSs = dfrSSs;
    bins<-(cutpts[2:length(cutpts)]+cutpts[1:(length(cutpts)-1)])/2;
    uYs<-sort(unique(tmpZCs[[yr]]));
    uFCs<-unique(tmpZCs[,c(sx,mt,sc)]);
    #cat("uFCs:\n")
    #print(uFCs);
    #cat("dfrSS:\n")
    #print(dfrSS);
    cat("SIZE_FREQUENCY_DATA  #required keyword\n",file=con);
    cat(optFit,          "            #objective function fitting option\n",file=con);
    cat(likeType,        "            #likelihood type\n",file=con);
    cat(likeWgt,         "            #likelihood weight\n",file=con);
    cat(tail_compression,"            #tail compression factors\n",file=con);
    cat(length(uYs),     "            #number of years of data\n",file=con);
    cat(unitsOut,        "            #units\n",file=con);
    cat(length(cutpts)," #number of size bin cutpoints\n",file=con);
    cat("#size bin cutpts (mm CW)\n",file=con);
    cat(cutpts,"\n",file=con);
    cat("#--------------\n",file=con);
    cat(nrow(uFCs),"    #number of factor combinations\n",file=con);
    for (iFC in 1:nrow(uFCs)){
      fc<-uFCs[iFC,];
      #cat("uFC[",iFC,",]:\n");
      #print(fc);
      cat(toupper(subForTCSAM(fc[[sx]],"ALL_SEX")),
          toupper(subForTCSAM(fc[[mt]],"ALL_MATURITY")),
          toupper(subForTCSAM(fc[[sc]],"ALL_SHELL")),"\n",file=con);
          cat("#use?  DMi  year    ss    ",bins,"\n",file=con);
          for (y in uYs){
            ids<-(tmpSSs[[yr]]==y)&
                  (tmpSSs[[sx]]==fc[[sx]])&
                  (tmpSSs[[mt]]==fc[[mt]])&
                  (tmpSSs[[sc]]==fc[[sc]]);
            idz<-(tmpZCs[[yr]]==y)&
                  (tmpZCs[[sx]]==fc[[sx]])&
                  (tmpZCs[[mt]]==fc[[mt]])&
                  (tmpZCs[[sc]]==fc[[sc]]);
            rw<-paste(tmpZCs[idz,vlc,drop=TRUE]*scl,collapse=" ");
            cat(1,0,y,tmpSSs[ids,ss,drop=TRUE],rw,"\n",sep="    ",file=con);#--default use?=1, Dirichlet Multinomial index=0
          }#--y
    }#--iFC
    
    return(invisible(con));
}


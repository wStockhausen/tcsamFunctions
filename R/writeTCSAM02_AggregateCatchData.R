#' 
#' @title Write aggregate catch data to a connection in TCSAM02 input format
#' 
#' @description Function to write aggregate catch data to a connection in TCSAM02 input format.
#' 
#' @param con - connection to use (default=stdout)
#' @param dfr - aggregate catch dataframe to write to connection
#' @param type - catch type ("ABUNDANCE" or "BIOMASS")
#' @param optFit - objective function fitting option (e.g., "BY_X","BY_XM")
#' @param likeType - likelihood type
#' @param likeWgt - likelihood weight (multiplier)
#' @param unitsIn - input catch units
#' @param unitsOut - output catch units
#' 
#' @return Invisibly returns the connection \code{con} to allow piping.
#' 
#' @details The input dataframe should have the following column names
#' \itemize{
#'   \item{year}
#'   \item{sex}
#'   \item{maturity}
#'   \item{shell condition}
#'   \item{value - value column}
#'   \item{cv - cv column}
#' }
#' 
#' @export
#' 
writeTCSAM02_AggregateCatchData<-function(con=stdout(),
                                          dfr=NULL,
                                          type=c("ABUNDANCE","BIOMASS"),
                                          optFit=c("BY_X","BY_XM"),
                                          likeType=c("NORM2","NORMAL","LOGNORMAL"),
                                          likeWgt = 1,
                                          unitsIn=ifelse(toupper(type)=="ABUNDANCE",
                                                       c("ONES","MILLIONS"),
                                                       c("KG","THOUSANDS_MT","MILLIONS_LBS","MT","LBS")),
                                          unitsOut=ifelse(toupper(type)=="ABUNDANCE",
                                                       c("MILLIONS","ONES"),
                                                       c("THOUSANDS_MT","MILLIONS_LBS","MT","KG","LBS"))){
    #--use only first element of these
    type     = toupper(type[1]);
    optFit   = toupper(optFit[1]);
    likeType = toupper(likeType[1]);
    likeWgt  = likeWgt[1];
    unitsIn  = toupper(unitsIn[1]);
    unitsOut = toupper(unitsOut[1]);
    
    if (type=="ABUNDANCE") scl = getScaleForAbundance(unitsIn,unitsOut);
    if (type=="BIOMASS")   scl = getScaleForBiomass(unitsIn,unitsOut);
    
    #--column names in dataframe
    yr  = "year";
    sx  = "sex";
    mt  = "maturity";
    sc  = "shell condition";
    vlc = "value";#--value column
    cvc = "cv";   #--cv column

    tmp = dfr;
    uYs<-sort(unique(tmp[[yr]]));
    uFCs<-unique(tmp[,c(sx,mt,sc)]);
    cat(paste0("AGGREGATE_",type),"     #required keyword\n",file=con);
    cat(optFit,     "            #objective function fitting option\n",file=con);
    cat(likeType,   "            #likelihood type\n",file=con);
    cat(likeWgt,    "            #likelihood weight\n",file=con);
    cat(length(uYs),"            #number of years\n",file=con,sep='');
    cat(unitsOut,   "            #units, catch abundance\n",file=con);
    cat(nrow(uFCs), "            #number of factor combinations\n",file=con);
    for (iFC in 1:nrow(uFCs)){
      fc<-uFCs[iFC,];
      cat(toupper(subForTCSAM(fc[[sx]],"ALL_SEX")),
          toupper(subForTCSAM(fc[[mt]],"ALL_MATURITY")),
          toupper(subForTCSAM(fc[[sc]],"ALL_SHELL")),"\n",file=con);
      cat("#use?  year    value    cv\n",file=con);
      for (y in uYs){
        idv<-(tmp[[yr]]==y)&
              (tmp[[sx]]==fc[[sx]])&
              (tmp[[mt]]==fc[[mt]])&
              (tmp[[sc]]==fc[[sc]]);
        vl<-tmp[idv,vlc,drop=TRUE];  #--unscaled catch values (ones )
        cv<-tmp[idv,cvc,drop=TRUE];  #--effective cv
        cat(1,y,scl*vl,cv,"\n",sep="    ",file=con);#--default: use?=1
      }#--y
    }#--fc
    
    return(invisible(con));
}

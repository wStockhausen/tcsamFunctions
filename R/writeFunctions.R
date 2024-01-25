#' 
#' @title Make substitutions for "undetermined"
#' 
#' @description Function to make substitutions for "undetermined".
#' 
#' @param x : character vector on which to make substitutions
#' @param str : string to substitute
#' 
#' @return modified character vector
#' 
#' @details None
#' 
subForTCSAM<-function(x,str){
xp <- ifelse(tolower(x)=="undetermined",str,x);
xp <- gsub(" ","_",xp,fixed=TRUE);
return(xp);
}

################################################################################
#' 
#' @title Get multiplicative scale factor from input abundance units to output units
#' 
#' @description Function to get multiplicative scale factor from input abundance units to output units.
#' 
#' @param unitsIn : units for input data
#' @param unitsOut : units for output data
#' 
#' @return multiplicative scale: output = scale * input.
#' 
#' @details None.
#' 
#' @export
#' 
getScaleForAbundance<-function(unitsIn,unitsOut){
  #--determine scaling from input units to output units
  getScale<-function(units){
      MILLIONS<-1000000;     #scale to convert to millions
      if (units=="MILLIONS") {
        scl <- 1.0/MILLIONS;
      } else if (units=="ONES") {
        scl <- 1.0;
      } else {
        msg<-paste0("\nERROR in tcsamFunction::getScaleForAbundance.",
                    "\ninput value for units ('",units,"') is invalid.",
                    "\nValid values are 'ONES','MILLIONS'.\n");
        stop(msg);
      }
      return(scl);
  }
  scaleIn = getScale(unitsIn);
  scaleOut = getScale(unitsOut);
  sclA = scaleOut/scaleIn;
  return(sclA);
}
################################################################################
#' 
#' @title Get multiplicative scale factor from input biomass units to output units
#' 
#' @description Function to get multiplicative scale factor from input biomass units to output units.
#' 
#' @param unitsIn : units for input data
#' @param unitsOut : units for output data
#' 
#' @return multiplicative scale: output = scale * input.
#' 
#' @details None.
#' 
#' @export
#' 
getScaleForBiomass<-function(unitsIn,unitsOut){
  #--determine scaling from input units to output units
  getScale<-function(units){
      MILLIONS<-1000000;     #scale to convert to millions
      LBStoKG <- 0.45359237; #multiplicative factor to get kg from lbs
      #--multiplicative factor from KG to units
      if (units=="THOUSANDS_MT") {
        scl <- 1.0/MILLIONS;
      } else if (units=="MILLIONS_LBS") {
        scl <- 1.0/(MILLIONS*LBStoKG);
      } else if (units=="KG") {
        scl <- 1.0;
      } else {
        msg<-paste0("\nERROR in tcsamFunctions::getScaleForBiomass",
                    "\ninput value for units ('",units,"') is invalid.",
                    "\nValid values are 'THOUSANDS_MT','MILLIONS_LBS','KG'.\n");
        stop(msg);
      }
  };
  scaleIn = getScale(unitsIn);
  scaleOut = getScale(unitsOut);
  sclB = scaleOut/scaleIn;
  return(sclB);
}
################################################################################
#' 
#' @title Create an input list for aggregate catch data
#' 
#' @description Function to create an input list for aggregate catch data.
#' 
#' @param type : type of aggregate catch data ("ABUNDANCE" or "BIOMASS")
#' @param dfr : the dataframe
#' @param cv : the default cv (NULL to use cv's from \code{dfr})
#' @param minErr : minimum assumed error in catch data (in 1's or kg, depending on \code{type})
#' @param optFit : objective function fitting option (e.g., "BY_X","BY_XM")
#' @param likeType : likelihood type ("NORM2", "NORMAL" or "LOGNORMAL")
#' @param likeWgt : likelihood multiplier
#' @param unitsIn : units for input data (possibilities depend on \code{type})
#' @param unitsOut : units for output data (possibilities depend on \code{type})
#' 
#' @return a list (see Details).
#' 
#' @details See below:
#' 
#' Output list has elements:
#' \itemize{
#' \item{dfr - dataframe}
#' \item{cv - default cv or dataframe with 'year', 'cv', and 'minErr' as columns (e.g., for fishery data)}
#' \item{minErr - minimum assumed error in catch data (in 1's)}
#' \item{optFit - objective function fitting option}
#' \item{likeType - likelihood type}
#' \item{likeWgt - likelihood multiplier}
#' \item{unitsIn - input units}
#' \item{unitsOut - output units}
#' }
#' 
#' Possible units are:
#' \itemize{
#'   \item{type = "ABUNDANCE": "ONES" or "MILLIONS"}
#'   \item{type = "BIOMASS": 'THOUSANDS_MT', 'MILLIONS_LBS', or 'KG'}
#' }
#' 
#' @export
#' 
inputList_AggregateCatchData<-function(type=c("ABUNDANCE","BIOMASS"),
                                       dfr=NULL,
                                       cv=0.05,
                                       minErr=100,
                                       optFit=c("BY_X","BY_XM"),
                                       likeType=c("NORM2","NORMAL","LOGNORMAL"),
                                       likeWgt=1,
                                       unitsIn =ifelse(toupper(type[1])=="ABUNDANCE",
                                                       c("ONES","MILLIONS"),
                                                       c('KG','THOUSANDS_MT', 'MILLIONS_LBS')),
                                       unitsOut=ifelse(toupper(type[1])=="ABUNDANCE",
                                                       c("MILLIONS","ONES"),
                                                       c('THOUSANDS_MT','KG','MILLIONS_LBS'))
                                   ){
    if(is.numeric(cv)) {
        cv     = cv[1];
        minErr = minErr[1];
    } else {
        #--alternative is a data.frame
        minErr = NULL;
    }
    lst=list(dfr=dfr,
            cv=cv,
            minErr=minErr,
            optFit=toupper(optFit[1]),
            likeType=toupper(likeType[1]),
            likeWgt=likeWgt[1],
            unitsIn =toupper(unitsIn[1]),
            unitsOut=toupper(unitsOut[1]));
    return(lst);
}
################################################################################
#' 
#' @title Create an input list for size comps data
#' 
#' @description Function to create an input list for size comps data.
#' 
#' @param dfrZCs : size comps dataframe
#' @param dfrSSs : sample sizes dataframe
#' @param cutpts : vector of cutpoints for size comps
#' @param tail_compression : 2-element vector giving tail compression factors
#' @param optFit : objective function fitting option (e.g., "BY_X","BY_XM")
#' @param likeType : likelihood type ("NORM2", "NORMAL" or "LOGNORMAL")
#' @param likeWgt -likelihood multiplier
#' @param unitsIn : units for input data ("ONES" or "MILLIONS")
#' @param unitsOut : units for output data ("ONES" or "MILLIONS")
#' 
#' @return a list (see Details).
#' 
#' @details See below:
#' 
#' Output list has elements:
#' \itemize{
#' \item{dfrZCs - size comps dataframe}
#' \item{dfrSSs - sample sizes dataframe}
#' \item{cutpts - cutpoints for size comps}
#' \item{tail_compression - 2-element vector giving tail compression factors}
#' \item{optFit - objective function fitting option}
#' \item{likeType - likelihood type}
#' \item{likeWgt - likelihood multiplier}
#' \item{unitsIn - input units}
#' \item{unitsOut - output units}
#' }
#' 
#' @export
#' 
inputList_SizeCompsData<-function(dfrZCs=NULL,
                                  dfrSSs=NULL,
                                  cutpts=NULL,
                                  tail_compression=c(0,0),
                                  optFit=c("BY_X","BY_XM"),
                                  likeType=c("MULTINOMIAL","DIRICHLET-MULTINOMIAL"),
                                  likeWgt=1,
                                  unitsIn =c("ONES","MILLIONS"),
                                  unitsOut=c("MILLIONS","ONES")
                               ){
    lst=list(dfrZCs=dfrZCs,
            dfrSSs=dfrSSs,
            cutpts=cutpts,
            tail_compression=tail_compression,
            optFit=toupper(optFit[1]),
            likeType=toupper(likeType[1]),
            likeWgt=likeWgt[1],
            unitsIn =toupper(unitsIn[1]),
            unitsOut=toupper(unitsOut[1]));
    return(lst);
}
################################################################################
#' 
#' @title Create an input list for effort data
#' 
#' @description Function to create an input list for effort data.
#' 
#' @param dfr : effort data dataframe
#' @param avgInterval : averaging interval (e.g., \code{[1992:-1]})
#' @param likeType : likelihood type ("NORM2","NORMAL",or "LOGNORMAL")
#' @param likeWgt : likelihood weight (multiplier)
#' @param unitsIn : input effort units ("ONES","MILLIONS")
#' @param unitsOut : output effort units ("ONES","MILLIONS")
#' 
#' @return a list (see Details).
#' 
#' @details See below:
#' 
#' The input dataframe \code{dfr} should have the following column names
#' \itemize{
#'   \item{year}
#'   \item{effort}
#' }
#' 
#' The output list has elements:
#' \itemize{
#'   \item{dfr - dataframe with columns 'year' and 'effort'}
#'   \item{avgInterval - averaging interval given in TCSAM02 time block format}
#'   \item{likeType - likelihood type}
#'   \item{likeWgt - likelihood weight}
#'   \item{unitsIn - units for input effort}
#'   \item{unitsOut - units for output effort}
#' }
#' 
#' @importFrom dplyr select
#' 
#' @export
#' 
inputList_EffortData<-function(dfr=NULL,
                               avgInterval="[1992:-1]",
                               likeType=c("NORM2","NORMAL","LOGNORMAL"),
                               likeWgt=1,
                               unitsIn =c("ONES","MILLIONS"),
                               unitsOut=c("MILLIONS","ONES")){
    lst=list(dfr=dfr %>% dplyr::select(year,effort),
             avgInterval=avgInterval,
             likeType=toupper(likeType[1]),
             likeWgt=likeWgt[1],
             unitsIn =toupper(unitsIn[1]),
             unitsOut=toupper(unitsOut[1]));
    return(lst);
}
################################################################################



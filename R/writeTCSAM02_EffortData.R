#' 
#' @title Write effort data to a connection in TCSAM02 input format
#' 
#' @description Function to write effort data to a connection in TCSAM02 input format.
#' 
#' @param con : connection to use (default=stdout)
#' @param dfr : effort data dataframe to write to connection
#' @param closed : vector of closure years (or NULL)
#' @param avgInterval : averaging interval (e.g., "\[1992:-1\]")
#' @param likeType : likelihood type ("NORM2","NORMAL",or "LOGNORMAL")
#' @param likeWgt : likelihood weight (multiplier)
#' @param unitsIn : input effort units ("ONES","MILLIONS")
#' @param unitsOut : output effort units ("ONES","MILLIONS")
#' 
#' @return Invisibly returns the connection \code{con} to allow piping.
#' 
#' @details The input values should be in ONES.
#' 
#' The input dataframe should have the following column names
#' \itemize{
#'   \item{year}
#'   \item{effort}
#' }
#' 
#' @importFrom dplyr filter
#' 
#' @import magrittr
#' 
#' @export
#' 
writeTCSAM02_EffortData<-function(con=stdout(),
                                  dfr=NULL,
                                  closed=NULL,
                                  avgInterval="[1992:-1]",
                                  likeType=c("NORM2","NORMAL","LOGNORMAL"),
                                  likeWgt = 1,
                                  unitsIn=c("ONES","MILLIONS"),
                                  unitsOut=c("ONES","MILLIONS")){
    #--use only first element of these
    likeType = toupper(likeType[1]);
    likeWgt  = likeWgt[1];
    unitsIn  = toupper(unitsIn[1]);
    unitsOut = toupper(unitsOut[1]);
    
    scl = getScaleForAbundance(unitsIn,unitsOut);
    
    #--column names in dataframe
    yr  <- "year";
    vlc = "effort";#--value column
    
    tmp<-dfr;
    if (!is.null(closed)) tmp %<>% dplyr::filter(!(year %in% closed));
    uYs<-sort(unique(tmp[[yr]]));
    cat("EFFORT_DATA    #required keyword\n",file=con);
    cat(avgInterval,"  #interval over which to average effort/fishing mortality\n",file=con);
    cat(likeType,   "      #likelihood type\n",file=con);
    cat(likeWgt,    "      #likelihood weight\n",file=con);
    cat(unitsOut,   "      #potlift units\n",file=con);
    cat(paste0(length(uYs),"    #number of years of directed effort data\n"),file=con);
    cat("#year	potlifts\n",file=con);
    for (y in uYs) {
      cat(y,scl*tmp$effort[tmp[[yr]]==y],"\n",file=con);
    }
    
    return(invisible(con));
}



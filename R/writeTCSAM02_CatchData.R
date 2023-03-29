#'
#' @title Write catch data type (retained catch, total catch, discard catch, index catch) to a connection
#'
#' @description Function to write catcht data type (retained catch, total catch, discard catch, index catch) to a connection.
#'
#' @param con : connection (default = stdout())
#' @param type : data type ("RETAINED","TOTAL","DISCARD" or "INDEX")
#' @param closed : vector of years when fishery was closed (if fleet is a fishery)
#' @param lstAbd : list with abundance data information (see \link{inputList_AggregateCatchData})
#' @param lstBio : list with biomass data information (see \link{inputList_AggregateCatchData})
#' @param lstZCs : list with size comps data information (see \link{inputList_SizeCompsData})
#'
#' @return Invisibly returns the connection to facilitate piping.
#'
#' @details 
#' 
#' See \link{inputList_AggregateCatchData} for list structure of \code{lstAbd} and \code{lstBio}.
#' See \link{inputList_SizeCompsData} for list structure of \code{lstZCs}
#' 
#' @import dplyr
#' @import magrittr
#'
#' @export
#'
writeTCSAM02_CatchData<-function(con=stdout(),
                                 type=c("RETAINED","TOTAL","DISCARD","INDEX"),
                                 closed=NULL,
                                 lstAbd=inputList_AggregateCatchData("ABUNDANCE"),
                                 lstBio=inputList_AggregateCatchData("BIOMASS"),
                                 lstZCs=inputList_SizeCompsData()
                           ){

  #--write flags for various quantities
  writeAbd<-!is.null(lstAbd) || !is.null(lstAbd$dfr);
  writeBio<-!is.null(lstBio) || !is.null(lstBio$dfr);
  writeZCs<-!is.null(lstZCs) || !is.null(lstZCs$dfrZCs);

  cat("CATCH_DATA     #required keyword\n",file=con);
  cat(writeAbd,"           #has aggregate catch abundance (numbers)\n",file=con);
  cat(writeBio,"           #has aggregate catch biomass (weight)\n",file=con);
  cat(writeZCs,"           #has size frequency data\n",file=con);

  #--catch abundance
  cat("#------------AGGREGATE CATCH DATA (NUMBERS)------------#\n",file=con);
  if (!writeAbd){
    cat("#------------no data-----------\n",file=con);
  } else {
    #--determine cv's and write data to connection
    dfr = lstAbd$dfr;
    if (!is.null(closed)) dfr %<>% dplyr::filter(!(year %in% closed));
    if (!is.null(lstAbd$cv)) {
        errScl = getScaleForAbundance("ONES",lstAbd$unitsIn);#--need minErr on same scale as values
        dfr %<>% dplyr::rowwise() %>% dplyr::mutate(cv=max(lstAbd$cv,errScl*lstAbd$minErr/value));#--effective cv
    }
    writeTCSAM02_AggregateCatchData(con=con,
                                    dfr=dfr,
                                    type="ABUNDANCE",
                                    optFit=lstAbd$optFit,
                                    likeType=lstAbd$likeType,
                                    likeWgt=lstAbd$likeWgt,
                                    unitsIn=lstAbd$unitsIn,
                                    unitsOut=lstAbd$unitsOut);
    rm(dfr)
  }#--writeAbd

  #--catch biomass
  cat("#------------AGGREGATE CATCH DATA (BIOMASS)------------#\n",file=con);
  if (!writeBio){
    cat("#------------no data-----------\n",file=con);
  } else {
    #--determine cv's and write data to connection
    dfr = lstBio$dfr;
    if (!is.null(closed)) dfr %<>% dplyr::filter(!(year %in% closed));
    if (!is.null(lstBio$cv)) {
        errScl = getScaleForBiomass("KG",lstBio$unitsIn);#--need minErr on same scale as values
        dfr %<>% dplyr::rowwise() %>% dplyr::mutate(cv=max(lstBio$cv,errScl*lstBio$minErr/value));#--effective cv
    }
    writeTCSAM02_AggregateCatchData(con=con,
                                    dfr=dfr,
                                    type="BIOMASS",
                                    optFit=lstBio$optFit,
                                    likeType=lstBio$likeType,
                                    likeWgt=lstBio$likeWgt,
                                    unitsIn=lstBio$unitsIn,
                                    unitsOut=lstBio$unitsOut);
    rm(dfr)
  }#--writeB

  #--size compositions
  cat("#------------NUMBERS-AT-SIZE DATA-----------\n",file=con);
  if (!writeZCs){
    cat("#------------no data-----------\n",file=con);
  } else {
    dfrZCs <- lstZCs$dfrZCs;
    if (!is.null(closed)) dfrZCs %<>% dplyr::filter(!(year %in% closed));
    dfrSSs <- lstZCs$dfrSSs;
    if (!is.null(closed)) dfrSSs %<>% dplyr::filter(!(year %in% closed));
    writeTCSAM02_SizeCompsData(con=con,
                               dfrZCs=dfrZCs,
                               dfrSSs=dfrSSs,
                               cutpts=lstZCs$cutpts,
                               tail_compression=lstZCs$tail_compression,
                               optFit=lstZCs$optFit,
                               likeType=lstZCs$likeType,
                               likeWgt=lstZCs$likeWgt,
                               unitsIn=lstZCs$unitsIn,
                               unitsOut=lstZCs$unitsOut);
  }#--writeZ
  return(invisible(con));
}

#'
#' @title Write fleet data (fishery or survey) to a TCSAM input file
#'
#' @description Function to write fleet data to a TCSAM input file.
#'
#' @param con - connection to write to
#' @param fleet - TCSAM fleet name
#' @param type - fleet type ('FISHERY' or 'SURVEY')
#' @param closed - vector of years when fishery was closed, if fleet is a fishery
#' @param lstRC - list with retained catch information (see [details])
#' @param lstTC - list with total catch information (see [details])
#' @param lstDC - list with discard catch information (see [details])
#' @param lstIC - list with index catch information (see [details])
#' @param lstEff - list with effort data information
#'
#' @return Invisibly returns the connection to facilitate piping.
#'
#' @details 
#' 
#' The input lists \code{lstRC}, \code{lstTC}, \code{lstDC}, and \code{lstIC}
#' must have the following structure, if not NULL:
#' \itemize{
#'   \item{lstAbd - NULL, or the \code{lstAbd} input to [writeTCSAM02_CatchData()] (see [inputList_AggregateCatchData()])}
#'   \item{lstBio - NULL, or the \code{lstAbd} input to[ writeTCSAM02_CatchData()] (see [inputList_AggregateCatchData()])}
#'   \item{lstZCs - NULL, or the \code{lstAbd} input to [writeTCSAM02_CatchData()] (see [inputList_SizeCompsData()])}
#' }
#' 
#' \code{lstEff} contains input values passed to [writeTCSAM02_EffortData()] and must have 
#' the following structure ([inputList_EffortData()]), if not NULL:
#' \itemize{
#'   \item{dfr - dataframe with columns 'year' and 'effort'}
#'   \item{avgInterval - averaging interval given in TCSAM02 time block format}
#'   \item{likeType - likelihood type}
#'   \item{likeWgt - likelihood weight}
#'   \item{unitsIn - units for input effort}
#'   \item{unitsOut - units for output effort}
#' }
#' 
#' @seealso [inputList_AggregateCatchData()], [inputList_SizeCompsData()], [inputList_EffortData()].
#' 
#' @export
#'
writeInputFile_FleetData<-function(con=stdout(),
                                   fleet=NULL,
                                   type=c("FISHERY","SURVEY"),
                                   closed=NULL,
                                   lstRC=NULL,
                                   lstTC=NULL,
                                   lstDC=NULL,
                                   lstIC=NULL,
                                   lstEff=NULL
                                  ){

    type = toupper(type[1]);
    
    #--write flags for various quantities
    hasRC = !is.null(lstRC);
    hasTC = !is.null(lstTC);
    hasDC = !is.null(lstDC);
    hasIC = !is.null(lstIC);
    hasEff = !is.null(lstEff);

    cat("#---------------------------\n",file=con);
    cat("#--TCSAM02 fleet data file--\n",file=con);
    cat("#---------------------------\n",file=con);
    cat(type,   "       #required keyword\n",file=con);
    cat(fleet,  "       #fleet name\n",file=con,sep='');
    cat(hasIC,  "       #has index catch data?\n",file=con);
    cat(hasRC,  "       #has retained catch data?\n",file=con);
    cat(hasDC,  "       #has observed discard catch data\n",file=con);
    cat(hasTC,  "       #has observed total catch data\n",file=con);
    cat(hasEff, "       #has effort data?\n",file=con);
    cat("#------------INDEX CATCH DATA------------	\n",file=con);
    if (!hasIC){
        cat("#---no data---\n",file=con);
    } else {
        con = writeTCSAM02_CatchData(con=con,
                                     type="INDEX",
                                     closed=NULL,        #--closures do not pertain
                                     lstAbd=lstIC$lstAbd,
                                     lstBio=lstIC$lstBio,
                                     lstZCs=lstIC$lstZCs);
    }#--hasIC
    
    cat("#------------RETAINED CATCH DATA------------#\n",file=con);
    if (!hasRC){
      cat("#---no data---\n",file=con);
    } else {
        con = writeTCSAM02_CatchData(con=con,
                                     type="RETAINED",
                                     closed=closed,
                                     lstAbd=lstRC$lstAbd,
                                     lstBio=lstRC$lstBio,
                                     lstZCs=lstRC$lstZCs);
    }#--hasRC

    cat("#------------DISCARD CATCH DATA------------#\n",file=con);
    if (!hasDC){
      cat("#---no data---\n",file=con);
    } else {
        con = writeTCSAM02_CatchData(con=con,
                                     type="DISCARD",
                                     closed=closed,
                                     lstAbd=lstDC$lstAbd,
                                     lstBio=lstDC$lstBio,
                                     lstZCs=lstDC$lstZCs);
    }#--hasDC

    cat("#------------TOTAL CATCH DATA------------#\n",file=con);
    if (!hasTC){
      cat("#---no data---\n",file=con);
    } else {
        con = writeTCSAM02_CatchData(con=con,
                                     type="TOTAL",
                                     closed=closed,
                                     lstAbd=lstTC$lstAbd,
                                     lstBio=lstTC$lstBio,
                                     lstZCs=lstTC$lstZCs);
    }#--hasTC

    cat("#------------EFFORT DATA------------#\n",file=con);
    if(!hasEff){
      cat("#---no data---\n",file=con);
    } else {
        con = writeTCSAM_EffortData(con=con,
                                    closed=closed,
                                    dfr=lstEff$dfr,
                                    avgInterval=lstEff$avgInterval,
                                    likeType=lstEff$likeType,
                                    likeWgt=lstEff$likeWgt,
                                    unitsIn=lstEff$unitsIn,
                                    unitsOut=lstEff$unitsOut);
    } #--hasEff

    return(invisible(con));
}

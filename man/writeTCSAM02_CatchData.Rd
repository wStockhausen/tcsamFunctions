% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/writeTCSAM02_CatchData.R
\name{writeTCSAM02_CatchData}
\alias{writeTCSAM02_CatchData}
\title{Write catch data type (retained catch, total catch, discard catch, index catch) to a connection}
\usage{
writeTCSAM02_CatchData(
  con = stdout(),
  type = c("RETAINED", "TOTAL", "DISCARD", "INDEX"),
  closed = NULL,
  lstAbd = inputList_AggregateCatchData("ABUNDANCE"),
  lstBio = inputList_AggregateCatchData("BIOMASS"),
  lstZCs = inputList_SizeCompsData()
)
}
\arguments{
\item{con}{: connection (default = stdout())}

\item{type}{: data type ("RETAINED","TOTAL","DISCARD" or "INDEX")}

\item{closed}{: vector of years when fishery was closed (if fleet is a fishery)}

\item{lstAbd}{: list with abundance data information (see \link{inputList_AggregateCatchData})}

\item{lstBio}{: list with biomass data information (see \link{inputList_AggregateCatchData})}

\item{lstZCs}{: list with size comps data information (see \link{inputList_SizeCompsData})}
}
\value{
Invisibly returns the connection to facilitate piping.
}
\description{
Function to write catcht data type (retained catch, total catch, discard catch, index catch) to a connection.
}
\details{
See \link{inputList_AggregateCatchData} for list structure of \code{lstAbd} and \code{lstBio}.
See \link{inputList_SizeCompsData} for list structure of \code{lstZCs}
}

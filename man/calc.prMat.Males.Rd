% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calc.prMat.Males.R
\name{calc.prMat.Males}
\alias{calc.prMat.Males}
\title{Function to calculate Pr(mature|z) for new shell males}
\usage{
calc.prMat.Males(z, shell_condition, z50 = 102.32849563, slp = 0.06440451)
}
\arguments{
\item{z}{: size (mm CW)}

\item{shell_condition}{: \sQuote{NEW_SHELL} or \sQuote{OLD_SHELL}}

\item{z50}{: size at 50\% maturity}

\item{slp}{: slope for logistic curve (\eqn{mm^{-1}})}
}
\value{
a vector of Pr(mature|z) corresponding to the input sizes and shell condition classification.
}
\description{
Logistic function for Pr(mature|size) for new shell males
based on Rugolo and Turnock's fitted logistic curve. Parameters are based on my fit to their curve.
}
\details{
Default parameter values are based on my fit to Rugolo and Turnock's fitted logistic curve. Old shell
crab are assumed to be terminally-molted, so \eqn{Pr(mature|z)=1} for all sizes \eqn{z}.
}

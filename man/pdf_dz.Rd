% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calcGammaPDFs.R
\name{pdf_dz}
\alias{pdf_dz}
\title{Evaluate a gamma pdf for \emph{molt increment} parameterized by mean and scale}
\usage{
pdf_dz(dz, mnDZ, beta)
}
\arguments{
\item{dz}{: the molt increment}

\item{mnDZ}{: the mean molt increment}

\item{beta}{: the gamma distribution scale parameter}
}
\value{
value of the corresponding value of the gamma distribution
}
\description{
Function to evaluate the pdf for \emph{molt increment} parameterized by the mean and scale.
}
\details{
Uses \link{pdf_gamma} to calculate the pdf, where
\eqn{\alpha} (the shape parameter) is given by
\deqn{\alpha = \frac{mnDZ}{\theta}}
and \eqn{\theta} (the scale parameter) is \code{beta}
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calcGammaPDFs.R
\name{plot_pdf}
\alias{plot_pdf}
\title{Plot a gamma pdf for post-molt size parameterized by mean and scale}
\usage{
plot_pdf(z_, mnZ, beta)
}
\arguments{
\item{z_}{\itemize{
\item premolt size
}}

\item{mnZ}{\itemize{
\item mean post-molt size
}}

\item{beta}{\itemize{
\item scale parameter
}}
}
\value{
ggplot2 plot object showing pdf and mean
}
\description{
Function to plot a gamma pdf for post-molt size parameterized by the mean and scale.
}
\details{
Expands \code{z_} by \code{0:100} to plot the pdf.
Overplots output from \link{pdf_z}, \link{pdf_dz}, and \link{pdf_gamma} to
demonstrate equivalency.

The underlying gamma pdf (\link{pdf_gamma}) is defined as
\deqn{pdf(x) = \frac{1}{\Gamma(\alpha) \theta^\alpha} x^(\alpha-1) e^{-\frac{x}{\theta}}}
where \code{x} is the \emph{molt increment} (i.e., post-molt size = \code{z_+x}),
\eqn{\theta =} \code{beta} is the scale parameter, and
\eqn{\alpha} is the shape parameter such that
\deqn{\alpha = \frac{\mu}{\theta}}
where \eqn{\mu} is the \emph{mean molt increment} (\code{mnZ-z_}).
}

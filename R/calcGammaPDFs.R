
#' 
#' @title Evaluate a gamma pdf for *post-molt size* given pre-molt size, parameterized by the mean post-molt size and scale
#' 
#' @description Function to evaluate the pdf for *post-molt size* given pre-molt size, parameterized by the 
#' mean post-molt size and scale.
#' 
#' @param zp : the post-molt size
#' @param z_ : the pre-molt size
#' @param mnZ : the mean post-molt size
#' @param beta : the gamma distribution scale parameter
#' 
#' @return value of the corresponding value of the gamma distribution
#' 
#' @details Converts the post-molt size to molt increment and mean molt increment 
#' and uses \link{pdf_dz} to calculate the corresponding value of the pdf.
#' 
#' @export
#' @md
#' 
pdf_z<-function(zp,z_,mnZ,beta){
 dz = zp-z_;
 mnDZ = mnZ-z_;
 return(pdf_dz(dz,mnDZ,beta));
}


#' 
#' @title Evaluate a gamma pdf for *molt increment* parameterized by mean and scale
#' 
#' @description Function to evaluate the pdf for *molt increment* parameterized by the mean and scale.
#' 
#' @param dz : the molt increment
#' @param mnDZ : the mean molt increment
#' @param beta : the gamma distribution scale parameter
#' 
#' @return value of the corresponding value of the gamma distribution
#' 
#' @details Uses \link{pdf_gamma} to calculate the pdf, where 
#' \eqn{\alpha} (the shape parameter) is given by 
#'   \deqn{\alpha = \frac{mnDZ}{\theta}}
#' and \eqn{\theta} (the scale parameter) is \code{beta}
#' @export
#' @md
#' 
pdf_dz<-function(dz,mnDZ,beta){
 alpha = mnDZ/beta;
 return(pdf_gamma(dz,alpha,beta))
}

#' 
#' @title Evaluate a gamma pdf
#' 
#' @description Function to evaluate a gamma pdf.
#' 
#' @param x - vector of values at which to evaluate the pdf
#' @param alpha - the gamma pdf shape parameter
#' @param theta - the gamma pdf scale parameter
#' 
#' @return vector of values of the gamma pdf evaluated at x
#' 
#' @details The underlying gamma pdf (\link{pdf_gamma}) is defined as
#'    \deqn{pdf(x) = \frac{1}{\Gamma(\alpha) \theta^\alpha} x^(\alpha-1) e^{-\frac{x}{\theta}}}
#' where \eqn{\alpha} (the shape parameter) is related to the mean \eqn{\mu} by
#'   \deqn{\alpha = \frac{\mu}{\theta}}
#' 
#' @export
#' @md
#' 
pdf_gamma<-function(x,alpha,theta){
 nll = lgamma(alpha) + alpha*log(theta)-(alpha-1)*log(x)+x/theta;
 return(exp(-nll));
}

#' 
#' @title Plot a gamma pdf for post-molt size parameterized by mean and scale
#' 
#' @description Function to plot a gamma pdf for post-molt size parameterized by the mean and scale.
#' 
#' @param z_ - premolt size
#' @param mnZ - mean post-molt size
#' @param beta - scale parameter
#' 
#' @return ggplot2 plot object showing pdf and mean
#' 
#' @details Expands \code{z_} by \code{0:100} to plot the pdf. 
#' Overplots output from \link{pdf_z}, \link{pdf_dz}, and \link{pdf_gamma} to
#' demonstrate equivalency.
#' 
#' The underlying gamma pdf (\link{pdf_gamma}) is defined as
#'    \deqn{pdf(x) = \frac{1}{\Gamma(\alpha) \theta^\alpha} x^(\alpha-1) e^{-\frac{x}{\theta}}}
#' where \code{x} is the *molt increment* (i.e., post-molt size = \code{z_+x}),  
#' \eqn{\theta =} \code{beta} is the scale parameter, and 
#' \eqn{\alpha} is the shape parameter such that
#'   \deqn{\alpha = \frac{\mu}{\theta}}
#' where \eqn{\mu} is the *mean molt increment* (\code{mnZ-z_}).
#' 
#' @import ggplot2
#' @importFrom tibble tibble
#' @importFrom tidyr pivot_longer
#' @importFrom wtsPlots getStdTheme
#' @export
#' @md
#' 
plot_pdf<-function(z_,mnZ,beta){
  zp  = z_+0:100;
  pz  = pdf_z(zp,z_,mnZ,beta);
  pdz = pdf_dz(zp-z_,mnZ-z_,beta);
  pg  = pdf_gamma(zp-z_,(mnZ-z_)/beta,beta);
  dfr = tibble::tibble(zp=zp,z_=z_,pz=pz,pdz=pdz,pg=pg) |> 
         tidyr::pivot_longer(cols=pz:pg);
  p = ggplot(dfr,aes(x=zp,y=value,colour=name)) + 
      geom_line() +
      geom_vline(xintercept=mnZ,colour="red") + 
      wtsPlots::getStdTheme();
  return(p);
}

#--test: plot_pdf(z_=23.6,mnZ=30.4874,beta=0.492);


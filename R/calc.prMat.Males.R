#'
#'@title Function to calculate Pr(mature|z) for new shell males
#'
#'@description Logistic function for Pr(mature|size) for new shell males
#'   based on Rugolo and Turnock's fitted logistic curve. Parameters are based on my fit to their curve.
#'
#'@param z : size (mm CW)
#'@param shell_condition : \sQuote{NEW_SHELL} or \sQuote{OLD_SHELL}
#'@param z50 : size at 50% maturity
#'@param slp : slope for logistic curve (\eqn{mm^{-1}})
#'
#'@return a vector of Pr(mature|z) corresponding to the input sizes and shell condition classification.
#'
#'@details Default parameter values are based on my fit to Rugolo and Turnock's fitted logistic curve. Old shell 
#'crab are assumed to be terminally-molted, so \eqn{Pr(mature|z)=1} for all sizes \eqn{z}.
#'
#'@export
#'
calc.prMat.Males<-function(z,
                           shell_condition,
                           z50=102.32849563,
                           slp=0.06440451){
    prM<-(1.0/(1.0+exp(-slp*(z-z50))))*(toupper(shell_condition)=="NEW_SHELL");
    prM<-prM+1.0*(toupper(shell_condition)=="OLD_SHELL");
    return(prM);
}

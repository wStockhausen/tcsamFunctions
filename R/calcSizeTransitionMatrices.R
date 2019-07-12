#'
#'@title Calculate size transition matrices for Tanner crab.
#'
#'@description Function to calculate sex-specific size transition matrices for Tanner crab.
#'
#'@param coeffs - list of sex-specific lists of growth parameters (a, b, beta)
#'@param sizes - vector of sizes at which to evaluate
#'
#'@return list with sex-specific size transition matrices
#'
#'@details Uses \code{calcSizeTransitionMatrix} to calculate sex-specific matrices
#'
#'@export
#'
calcSizeTransitionMatrices<-function(coeffs=list(  MALE=list(a=0.42577,b=0.971389,beta=0.75,type='LinearGrowthIncrement'),
                                                 FEMALE=list(a=0.70000,b=0.883118,beta=0.75,type='LinearGrowthIncrement')),
                                     sizes=seq(from=27.5,to=182.5,by=5)){
   nms<-names(coeffs);
   res<-list();
   for (nm in nms){
       prGr<-calcSizeTransitionMatrix(coeffs[[nm]],sizes=sizes);
       res[[nm]]<-prGr;
   }
   
   return(invisible(prGr));
}
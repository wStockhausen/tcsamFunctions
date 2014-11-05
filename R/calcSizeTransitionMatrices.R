#'
#'@title Calculate size transition matrices for Tanner crab.
#'
#'@description Function to calculate sex-specific size transition matrices for Tanner crab.
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
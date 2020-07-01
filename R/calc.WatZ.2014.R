#'
#'@title Convert Tanner crab size to weight based on the 2014 assessment (Rugolo and Turnock)
#'
#'@description Power-law functions for Tanner crab weight-at-size (g) by sex and maturity
#'   based on the 2014 assessment (Rugolo and Turnock's parameters).
#'   
#'@param z : vector of sizes (mm CW)
#'@param sex : 'MALE' or 'FEMALE'
#'@param maturity : 'IMMATURE' or 'MATURE'
#'@param male : list of regression coefficients for immature, mature crab (see details)
#'@param female : list of regression coefficients for immature, mature crab (see details)
#'
#'@return a vector of weights corresponding to the input sizes
#'
#'@details 'male' and 'female' are lists with named elements 'a' and 'b', whose values
#'reflect the parameters in the eq: w = a*z^b
#'
#'@export
#'
#-----------------------------------------------------------
calc.WatZ.2014<-function(z,sex,maturity,
                        male=list(immature=list(a=0.00016,b=3.136),
                                    mature=list(a=0.00016,b=3.136)),
                        female=list(immature=list(a=0.00064,b=2.794),
                                      mature=list(a=0.00034,b=2.956))){
    idx.m<-toupper(sex)==  'MALE';
    idx.f<-toupper(sex)=='FEMALE';
    idx.imm<-toupper(maturity)=='IMMATURE';
    idx.mat<-toupper(maturity)==  'MATURE';
    
    wgt<-(male$immature$a*z^male$immature$b)*(idx.m&idx.imm);
    wgt<-wgt+(male$mature$a*z^male$mature$b)*(idx.m&idx.mat);
    wgt<-wgt+(female$immature$a*z^female$immature$b)*(idx.f&idx.imm);
    wgt<-wgt+(female$mature$a*z^female$mature$b)    *(idx.f&idx.mat);
    return(wgt);
}

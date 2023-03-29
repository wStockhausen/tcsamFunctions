#'
#'@title Convert Tanner crab size to weight using the current (2015) relationships
#'
#'@description Power-law functions for Tanner crab weight-at-size (g) by sex and maturity
#'   based on 2015 NMFS survey.
#'   
#'@param z : vector of sizes (mm CW)
#'@param sex : vector of sexes ('MALE' or 'FEMALE')
#'@param maturity : vector of maturity states ('IMMATURE', 'MATURE', or 'UNDETERMINED')
#'@param male : list of regression coefficients for immature, mature crab (see details)
#'@param female : list of regression coefficients for immature, mature crab (see details)
#'
#'@return a vector of weights in g corresponding to the input sizes, sexes, and maturity states
#'
#'@details 'male' and 'female' are lists with named elements 'a' and 'b', whose values
#'reflect the parameters in the eq:  \eqn{w = a \cdot z^b}
#'
#'@export
#'
#-----------------------------------------------------------
calc.WatZ<-function(z,sex,maturity,
                    male=list(immature    =list(a=0.00027,b=3.022134),
                              mature      =list(a=0.00027,b=3.022134),
                              undetermined=list(a=0.00027,b=3.022134)),
                    female=list(immature=list(a=0.000562,b=2.816928),
                                mature  =list(a=0.000441,b=2.898686))){
    idx.m<-toupper(sex)=='MALE';
    idx.f<-toupper(sex)=='FEMALE';
    idx.imm<-toupper(maturity)=='IMMATURE';
    idx.mat<-toupper(maturity)=='MATURE';
    idx.und<-toupper(maturity)=='UNDETERMINED';
    
    wgt<-(male$immature$a    *z^male$immature$b    )*(idx.m&idx.imm)+
         (male$mature$a      *z^male$mature$b      )*(idx.m&idx.mat)+
         (male$undetermined$a*z^male$undetermined$b)*(idx.m&idx.und)+
         (female$immature$a  *z^female$immature$b  )*(idx.f&idx.imm)+
         (female$mature$a    *z^female$mature$b    )*(idx.f&idx.mat);
    
    return(wgt);
}

#'
#'@title Calculate the probability at crab size z for the molt to maturity, pr(mature|z).
#'
#'@description Function to calculate the probability at crab size z for the 
#'molt to maturity, pr(mature|z). Note that it is up to the user to define whether
#'"size" refers to pre-molt or post-molt size. If pre-molt, then pr(mature|z) is the
#'probability of molting to maturity (terminal molt) on the NEXT molt whereas if size is post-molt
#'then pr(mature|z) is the probability that a crab that just molted to size z also matured.
#'
#'@param coeffs - a named list identifying function type and coefficient values
#'@param sizes - sizes at which to evaluate the function
#'
#'@details coeffs is a named list with the 1st element 'type' identifying the function type (e.g. 'logistic')
#'and remaining elements as function coefficients (e.g., slope, z50)
#'
#'@export
#'
calcPrMoltToMaturity<-function(coeffs=list(type='logistic',z50=130,slope=0.1),
                               sizes=seq(from=25,to=185,by=1)){
    
    if (coeffs$type=='logistic'){
        pr<-1/(1+exp(-coeffs$slope*(sizes-coeffs$z50)));
    }
    
    return(pr);
}
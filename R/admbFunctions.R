#'
#' @title Implemetation of ADMB function "mfexp"
#'
#' @description Function that implements ADMB function "mfexp"
#'
#' @param x - vector
#'
#' @return exp(x)
#'
#' @details R version of ADMB mfexp.
#'
#' @export
#'
mfexp<-function(x){
    return(exp(x));
}

#'
#' @title Implemetation of ADMB function "square"
#'
#' @description Function that implements ADMB function "square"
#'
#' @param x - vector
#'
#' @return x*x
#'
#' @details R version of ADMB square.
#'
#' @export
#'
square<-function(x){
    return(x*x);
}

#'
#' @title Implemetation of ADMB function "elem_prod"
#'
#' @description Function that implements ADMB function "elem_prod"
#'
#' @param x - vector
#' @param y - vector
#'
#' @return elem_prod(x,y)
#'
#' @details R version of ADMB elem_prod. x and y must be the same length, or
#' one or both must have only 1 element.
#'
#' @export
#'
elem_prod<-function(x,y){
    z = x*y;
    return(z);
}


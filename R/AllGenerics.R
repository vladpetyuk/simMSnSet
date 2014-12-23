

#' Simulates components of quantitative LC-MS data
#' @param object an instance of simulatable classs. It should contain 
#'          some parameter and underlying statistical distribution that
#'          generate random samplings.         
#' @return The intention is that it returns the instance of the same class,
#'          but now with filled slots that needs to be simulated.
#' @export
setGeneric("simulate",
           function(object, ...)
               standardGeneric("simulate"))


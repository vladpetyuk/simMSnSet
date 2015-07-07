#' @include AllGenerics.R
#' @include proteins.R
#' @include phenotypes.R
NULL


#' Constructor for Experiment class object
#' @param proteins an instance of \link{Proteins-class}
#' @param phenotypes an instance of \link{Phenotypes-class}
#' @return \link{Experiment-class} object
#' @name Experiment-constructor
#' @export
Experiment <- function(proteins, phenotypes){
    new("Experiment",
        proteins=proteins,
        phenotypes=phenotypes,
        simulated=FALSE)
}




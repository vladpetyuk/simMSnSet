# Definition of phenotypes class and corresponding contructors
# Somehow need to be conformed with pData of eSet (MSnSet)
# In essense Phenotypes is affine to phenoData

#' @include AllGenerics.R
NULL


#' Constructor for Phenotypes class object
#' @param originalPhenotypes character vector with original 
#'          (intended) phenotype assignments
#' @param proportionOutlying proportion of outlying samples. 
#'          This will be used in simulation as parameter of Bernoulli distribution.
#' @return \link{Phenotypes-class} object
#' @name Phenotypes-constructor
#' @export
Phenotypes <- function(originalPhenotypes=rep(LETTERS[1:3],each=4),
                       proportionOutlying = 0.1,
                       meanFoldChange=1,
                       pipettingAccuracy=0.15)
{
    #
    sample.dup <- unique(originalPhenotypes[duplicated(originalPhenotypes)])
    sampleNames <- make.unique(c(sample.dup, originalPhenotypes))[-seq_along(sample.dup)]
    #
    phenotypesObj <- new("Phenotypes", 
                         n=length(originalPhenotypes),
                         originalPhenotypes=originalPhenotypes, 
                         proportionOutlying=proportionOutlying,
                         meanFoldChange=meanFoldChange,
                         pipettingAccuracy=pipettingAccuracy,
                         sampleNames=sampleNames,
                         simulated=FALSE)
}





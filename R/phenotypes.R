# Definition of phenotypes class and corresponding contructors
# Somehow need to be conformed with pData of eSet (MSnSet)
# In essense Phenotypes is affine to phenoData

#' @include AllGenerics.R
NULL


#' Phenotypes class.
#'
#' Details about Phenotypes class.
#' 
#' Class holding information about the samples
#' 
#' @slot n number of samples
#' @slot proportionOutlying proportion of outlying or 
#'          erroneusly assigned samples.  It is not uncommong that in
#'          large-scale experiments up to ~5% samples appear as outliers.
#' @slot original character vector with original or intended 
#'          phenotype (group) assignment
#' @slot simulated character vector with phenotype assigments after pretending
#'          that some of the samples are actually outliers
#' @slot sampleNames character vector with sample names
#' 
#' @examples
#' # 10 samples (5 controls + 5 cases)
#' # Note, the number of outlying samples is somewhat extreme - 50%.
#' # This is for demonstration purpose only.
#' ph <- Phenotypes(original=c('Ctrl', 'Ctrl', 'Ctrl', 'Ctrl', 'Ctrl',
#'                             'Case', 'Case', 'Case', 'Case', 'Case'),
#'                  proportionOutlying=0.5)
#' ph <- simulate(ph, seed=123)
#' show(ph)
#'  
#' @exportClass Phenotypes
#' 
setClass(Class="Phenotypes",
         representation(
             n="integer",
             proportionOutlying="numeric",
             original="character",
             simulated="character",
             sampleNames="character")
)



#' Constructor for Phenotypes class object
#' @param original character vector with original 
#'          (intended) phenotype assignments
#' @param proportionOutlying proportion of outlying samples. 
#'          This will be used in simulation as parameter of Bernoulli distribution.
#' @return \link{Phenotypes-class} object
#' @name Phenotypes-constructor
#' @export
Phenotypes <- function(original=rep(LETTERS[1:3],each=4),
                       proportionOutlying = 0.1)
{
    #
    sample.dup <- unique(original[duplicated(original)])
    sampleNames <- make.unique(c(sample.dup, original))[-seq_along(sample.dup)]
    #
    phenotypesObj <- new("Phenotypes", 
                         n=length(original),
                         original=original, 
                         proportionOutlying=proportionOutlying,
                         sampleNames=sampleNames)
}




#' @describeIn Phenotypes Selects which samples are going to be outliers given
#'              the value of the \code{proportionOutlying} slot.
#' @param seed an integer for \code{set.seed}
#' 
setMethod(
    "simulate",
    signature(object="Phenotypes"),
    definition=function(object, seed=NULL)
    {
        probs <- rep(object@proportionOutlying, object@n)
        outlierNames <- paste('Out', seq_len(object@n), sep='')
        if(!is.null(seed)) set.seed(seed) # for the sake for reproducibility
        idx <- sapply(probs, rbinom, n=1, size=1) == 1
        simulated <- object@original
        simulated[idx] <- outlierNames[idx]
        object@simulated <- simulated
        return(object)
    }
)



#' @describeIn Phenotypes Show method for a Phenotype object
#' 
setMethod(
    "show",
    signature(object="Phenotypes"),
    definition=function(object){
        cat("An object of class \"Phenotypes\"\n")
        cat("Number of samples:", object@n, "\n")
        cat("Original phenotypes:", object@original, "\n")
        cat("Proportion of outliers:", object@proportionOutlying, "\n")
        cat("Actual phenotypes:", 
            if(identical(object@simulated, character(0)))
                "not simulated"
            else
                object@simulated, "\n")
    }
)




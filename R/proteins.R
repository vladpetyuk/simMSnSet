# Definition of proteins class and corresponding contructors
# Somehow need to be conformed with fData of eSet (MSnSet)
# In essense Proteins is affine to featureData

#' @include AllGenerics.R
NULL



#' Proteins class.
#'
#' Details about Proteins class.
#' 
#' Class holding the information about the proteins and peptides
#' 
#' @slot nProt number of proteins
#' @slot nPep number of peptides
#' @slot lambda parameter of Poisson distribution controlling the number of 
#'          peptides per protein during the simulation
#' @slot meanIoinizationIntensity a numeric value of mean peptide intensity
#' @slot pepCounts integer vector with number of peptides per each protein
#' @slot protNames character vector with generated protein names
#' @slot pepNames character vector with generacted peptide names
#' @slot ionizationEff numeric vector with simulated ionization efficiencies 
#'          for the peptides
#' @slot features data.frame containing peptide and protein names
#' 
#' @examples
#' pr <- Proteins(nProt=10)
#' pr <- simulate(pr)
#' show(pr)
#'  
#' @exportClass Proteins
#' 
setClass(Class="Proteins",
         representation(
             nProt="integer",
             nPep="integer",
             lambda="numeric",
             meanIoinizationIntensity="numeric",
             pepCounts="integer",
             protNames="character",
             pepNames="character",
             ionizationEff="numeric",
             features="data.frame")
)



#' Constructor for Proteins class object
#' @param nProt number of proteins
#' @param lambda parameter of Poisson distribution 
#'          controlling number of peptides per protein
#' @param meanIoinizationIntensity numeric value of an average 
#'          peptide ionization intensity
#' @param pepCounts (optional) integer vector with number of peptides 
#'          per each protein. If provided then \code{lambda} is ignored.
#' @return \link{Proteins-class} object
#' @name Proteins-constructor
#' @export
Proteins <- function(nProt=20, 
                     lambda=2, 
                     meanIoinizationIntensity=1e6, 
                     pepCounts=NULL, ...)
{
    #
    protNames <- paste("prot", seq_len(nProt), sep='')
    #
    protObj <- new("Proteins", 
                   nProt=as.integer(nProt),
                   protNames=protNames,
                   lambda=lambda,
                   meanIoinizationIntensity=meanIoinizationIntensity)
    #
    # in case peptide counts for each protein are provided
    if(!is.null(pepCounts)){
        pepNames <- paste(rep(protNames, pepCounts), 
                          unlist(sapply(pepCounts, seq_len)), sep='_pep')
        features <- data.frame(prot=rep(protNames, pepCounts), 
                               pep=pepNames)
        protObj@pepCounts <- pepCounts
        protObj@pepNames <- pepNames
        protObj@features <- features
    }
    return(protObj)
}



#' @describeIn Proteins Simulates the number of peptides per protein 
#'              (if not provided) and ionization efficiencies of each peptide.
#' @param seed an integer for \code{set.seed}
#' 
setMethod(
    "simulate",
    signature(object="Proteins"),
    definition=function(object, seed=NULL)
    {
        if(!is.null(seed)) set.seed(seed) # for the sake for reproducibility
        pepCounts <- rpois(object@nProt, object@lambda) + 1
        nPep <- sum(pepCounts)
        pepNames <- paste(rep(object@protNames, pepCounts), 
                          unlist(sapply(pepCounts, seq_len)), sep='_pep')
        features <- data.frame(prot=rep(object@protNames, pepCounts), 
                               pep=pepNames)
        #
        ionizationEff <- rlnorm(nPep,
                                meanlog = log(object@meanIoinizationIntensity))
        #
        object@pepCounts <- as.integer(pepCounts)
        object@nPep <- as.integer(nPep)
        object@pepNames <- pepNames
        object@features <- features
        object@ionizationEff <- ionizationEff
        return(object)
    }
)



#' @describeIn Proteins Show method for a Proteins object
#' 
setMethod(
    "show",
    signature(object="Proteins"),
    definition=function(object){
        cat("An object of class \"Proteins\"\n")
        cat("Number of proteins:", object@nProt, "\n")
        cat("Number of peptides:", 
            if(identical(object@nPep, integer(0)))
                "not simulated"
            else
                object@nPep, "\n")
    }
)



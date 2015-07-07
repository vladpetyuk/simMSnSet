# Definition of proteins class and corresponding contructors
# Somehow need to be conformed with fData of eSet (MSnSet)
# In essense Proteins is affine to featureData

#' @include AllGenerics.R
NULL


#' Constructor for Proteins class object
#' @param nProt number of proteins
#' @param lambda parameter of Poisson distribution 
#'          controlling number of peptides per protein
#' @param meanIoinizationIntensity numeric value of an average 
#'          peptide ionization intensity
#' @param pepCounts (optional) integer vector with number of peptides 
#'          per each protein. If provided then \code{lambda} is ignored.
#' @param noise measurement noise. So far just a single value for all
#'          peptides/proteins.
#' @return \link{Proteins-class} object
#' @name Proteins-constructor
#' @export
Proteins <- function(nProt=20, 
                     propChanging=0.2,
                     fdrPeptide=0.05,
                     fdrMatching=0.01,
                     lambda=2, 
                     meanIoinizationIntensity=1e6,
                     missThreshQuantile=0.05,
                     missThreshSharpness=10,
                     noise=1,
                     pepCounts=NULL, ...)
{
    #
    protNames <- paste("prot", seq_len(nProt), sep='')
    #
    protObj <- new("Proteins", 
                   nProt=as.integer(nProt),
                   propChanging=propChanging,
                   fdrPeptide=fdrPeptide,
                   fdrMatching=fdrMatching,
                   protNames=protNames,
                   lambda=lambda,
                   noise=noise,
                   meanIoinizationIntensity=meanIoinizationIntensity,
                   missThreshQuantile=missThreshQuantile,
                   missThreshSharpness=missThreshSharpness, ...)
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




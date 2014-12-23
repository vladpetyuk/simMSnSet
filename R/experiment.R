#' @include AllGenerics.R
#' @include proteins.R
#' @include phenotypes.R
NULL



#' Experiment class.
#'
#' Details about Experiment class.
#' 
#' The central class of simMSnSet package
#' 
#' @slot phenotypes Phenotypes and Samples
#' @slot proteins Peptides and Proteins
#' @slot intensities Simulated intensities
#' @slot simulated A logical flag indicating if the experiment has been simulated or not
#' 
#' @examples
#' # 10 proteins measured across 6 samples (3 controls + 3 cases)
#' pr <- Proteins(nProt=10)
#' ph <- Phenotypes(original=c('Control', 'Control', 'Control', 
#'                             'Case', 'Case', 'Case'))
#' ex <- Experiment(pr, ph)
#' ex <- simulate(ex)
#' image(ex)
#'  
#' @exportClass Experiment
#' 
setClass(Class="Experiment",
         representation(
             phenotypes="Phenotypes",
             proteins="Proteins",
             intensities="matrix",
             simulated="logical")
)



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



#' @describeIn Experiment Simulates the entire experiment
#' @param seed an integer for \code{set.seed}
#' 
setMethod(
    "simulate",
    signature(object="Experiment"),
    definition=function(object, seed=NULL)
    {
        if(!is.null(seed)) set.seed(seed) # for the sake for reproducibility
        object@phenotypes <- simulate(object@phenotypes)
        object@proteins <- simulate(object@proteins)
        #
        intensities <- matrix(rlnorm(object@phenotypes@n * 
                                     object@proteins@nPep),
                              ncol=object@phenotypes@n)    
        colnames(intensities) <- object@phenotypes@sampleNames
        rownames(intensities) <- object@proteins@pepNames
        names(dimnames(intensities)) <- c("peptides","samples")
        object@intensities <- intensities
        object@simulated <- TRUE
        return(object)
    }
)



#' @describeIn Experiment Show method for an Experiment object
#' 
setMethod(
    "show",
    signature(object="Experiment"),
    definition=function(object)
    {
        cat("___Phenotypes___\n")
        show(object@phenotypes)
        cat("___Proteins___\n")
        show(object@proteins)
    }
)



#' @describeIn Experiment Showing simulated experiment as an image
#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot geom_raster scale_fill_gradient2 aes
#' @exportMethod image
setMethod(
    "image",
    "Experiment",
#     signature(object="Experiment"),
    function(x, ...)
    {
        stopifnot(x@simulated)
        int <- x@intensities
        # log2 transform, zero-centering
        int <- log2(int)
        int <- sweep(int, 1, rowMeans(int), '-')
        xm <- melt(int)
        xm$peptides <- ordered(xm$peptides, levels=rev(x@proteins@pepNames))
        #
        ggplot(xm, aes(x=samples, y=peptides, fill=value)) +
            geom_raster() +
            scale_fill_gradient2()
    }
)


'YYYY-MM-DD' 

# ' @describeIn Experiment Coerce simulated experiment into 
# '              \link[MSnbase]{MSnSet-class} object
# ' @exportMethod as
#' @importFrom MSnbase MSnSet experimentData<-
#' @importClassesFrom MSnbase MIAPE
setAs(from="Experiment", 
      to="MSnSet", 
      def=function(from){
          
                if(!from@simulated){
                    warning("The Expriment has not been simulated!. Returning NA")
                    return(NA)
                }
                
                # MIAPE
                # this will contain simulation parameters
                miape <- new("MIAPE",
                             title="Simulation of MSnSet object",
                             name=Sys.info()['effective_user'],
                             samples=list(sampleNames=from@phenotypes@sampleNames),
                             dateStamp=format(Sys.Date(), "%Y-%m-%d"),
                             other=list(lambda=from@proteins@lambda))
                fd <- from@proteins@features
                rownames(fd) <- fd$pep
                pd <- data.frame(phenotype=from@phenotypes@original, 
                                 row.names=from@phenotypes@sampleNames)
                obj <- MSnSet(exprs=from@intensities, 
                              fData=fd,
                              pData=pd)
                experimentData(obj) <- miape
                if(validObject(obj)){
                    return(obj)
                }else{
                    warning("Invalid MIAPE info! Returning NA")
                    return(NA)
                }
          }
)




# library("MSnbase")
# M <- matrix(rnorm(12), 4)
# pd <- data.frame(otherpdata = letters[1:3])
# fd <- data.frame(otherfdata = letters[1:4])
# x0 <- MSnSet(M, fd, pd)
# show(x0)



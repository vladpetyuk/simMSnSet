

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




#' @describeIn Phenotypes Show method for a Phenotype object
#' 
setMethod(
    "show",
    signature(object="Phenotypes"),
    definition=function(object){
        cat("An object of class \"Phenotypes\"\n")
        cat("Number of samples:", object@n, "\n")
        cat("Original phenotypes:", object@originalPhenotypes, "\n")
        cat("Proportion of outliers:", object@proportionOutlying, "\n")
        cat("Actual phenotypes:", 
            if(!object@simulated) "not simulated"
            else object@simulatedPhenotypes, "\n")
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



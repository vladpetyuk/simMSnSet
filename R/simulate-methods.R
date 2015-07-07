# The core. Simulation methods.

#' @describeIn Proteins Simulates the number of peptides per protein 
#'              (if not provided) and ionization efficiencies of each peptide.
#' @param object \link{Proteins-class} object
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
        #
        probs <- rep(object@propChanging, object@nProt)
        object@changingProtIdx <- sapply(probs, rbinom, n=1, size=1) == 1
        return(object)
    }
)




#' @describeIn Phenotypes Selects which samples are going to be outliers given
#'              the value of the \code{proportionOutlying} slot.
#' @param object \link{Phenotypes-class} object
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
        simulatedPhenotypes <- object@originalPhenotypes
        simulatedPhenotypes[idx] <- outlierNames[idx]
        object@simulatedPhenotypes <- simulatedPhenotypes
        object@sampleBiases <- rlnorm(object@n, 
                                      meanlog = 0, 
                                      sdlog = object@pipettingAccuracy)
        object@simulated <- TRUE
        return(object)
    }
)





#' @describeIn Experiment Simulates the entire experiment
#' @param seed an integer for \code{set.seed}
#' @param object \link{Experiment-class} object
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
        # Simulating intensities as matrix.
        # The key component of the "Experiment".
        # Y <- yIonization * yGroupMeans * yNoise * ySysErr * (1-idxOutlier) + 
        #      yOutlyingMeasurements * idxOutlier
        #
        simIonization <- matrix(rep(object@proteins@ionizationEff, 
                                    object@phenotypes@n), 
                                ncol=object@phenotypes@n, byrow=F)
        #
        simNoise <- matrix(rlnorm(object@phenotypes@n * object@proteins@nPep, 
                                  sdlog = object@proteins@noise),
                              ncol=object@phenotypes@n)
        #
        simSampleBiases <- matrix(rep(object@phenotypes@sampleBiases, 
                                      object@proteins@nPep), 
                                  ncol=object@phenotypes@n, byrow=TRUE)
        #
        simOutliers <- matrix(rbinom(object@proteins@nPep * object@phenotypes@n, 1, 
                                     prob=object@proteins@fdrMatching),
                              ncol=object@phenotypes@n)
        #
        simOutlyingMeasurements <- matrix(
            rlnorm(object@proteins@nPep*object@phenotypes@n, 
                   meanlog = log(object@proteins@meanIoinizationIntensity), 
                   sdlog = object@proteins@noise), # Dispertion of outliers. This may be tweaked later.
            nrow=object@proteins@nPep, 
            ncol=object@phenotypes@n)
        #
        # PHENOTYPE_ORIGINAL -- ph2@original
        # PHENOTYPE -- ph2@simulated
        # N_PROT -- pr2@nProt
        # PROPORTION_CHANGING_PROTEINS -- pr2@propChanging
        # FDR_PEPTIDE -- pr2@fdrPeptide
        simMeans <- .generate_group_means(object@phenotypes@originalPhenotypes,
                                              object@phenotypes@simulatedPhenotypes,
                                              object@proteins@nProt,
                                              object@proteins@pepCounts,
                                              object@proteins@features,
                                              object@proteins@propChanging,
                                              object@proteins@changingProtIdx,
                                              object@phenotypes@meanFoldChange,
                                              object@proteins@fdrPeptide)
        #
        simFinal <- simIonization * simNoise * simMeans * simSampleBiases * 
                    (1 - simOutliers) + simOutlyingMeasurements * simOutliers
        #
        if(identical(object@proteins@missThreshAbsolute, numeric()))
            object@proteins@missThreshAbsolute <- 
                quantile(simFinal, object@proteins@missThreshQuantile)
        thresh <- object@proteins@missThreshAbsolute
        sharp <- object@proteins@missThreshSharpness
        pMiss <- 1 - 1/(1+(thresh/simFinal)^sharp)
        idxMissing <- apply(pMiss, c(1,2), rbinom, n=1, size=1) == 1
        simMissing <- matrix(1, ncol=object@phenotypes@n, nrow=object@proteins@nPep)
        simMissing[idxMissing] <- NA
        #
        simFinal <- simFinal * simMissing
        # 1
#         colnames(simIonization) <- object@phenotypes@sampleNames
#         rownames(simIonization) <- object@proteins@pepNames
#         names(dimnames(simIonization)) <- c("peptides","samples")
#         object@simIonization <- simIonization
        #
        # 2
        colnames(simNoise) <- object@phenotypes@sampleNames
        rownames(simNoise) <- object@proteins@pepNames
        names(dimnames(simNoise)) <- c("peptides","samples")
        object@simNoise <- simNoise
        #
        # 3
        colnames(simMeans) <- object@phenotypes@sampleNames
        rownames(simMeans) <- object@proteins@pepNames
        names(dimnames(simMeans)) <- c("peptides","samples")
        object@simMeans <- simMeans
        #
        # 4
        colnames(simSampleBiases) <- object@phenotypes@sampleNames
        rownames(simSampleBiases) <- object@proteins@pepNames
        names(dimnames(simSampleBiases)) <- c("peptides","samples")
        object@simSampleBiases <- simSampleBiases
        #
        # 5
        colnames(simOutliers) <- object@phenotypes@sampleNames
        rownames(simOutliers) <- object@proteins@pepNames
        names(dimnames(simOutliers)) <- c("peptides","samples")
        object@simOutliers <- simOutliers
        #
        # 6
        colnames(simMissing) <- object@phenotypes@sampleNames
        rownames(simMissing) <- object@proteins@pepNames
        names(dimnames(simMissing)) <- c("peptides","samples")
        object@simMissing <- simMissing
        #
        # Z
        colnames(simFinal) <- object@phenotypes@sampleNames
        rownames(simFinal) <- object@proteins@pepNames
        names(dimnames(simFinal)) <- c("peptides","samples")
        object@simFinal <- simFinal
        #
        object@simulated <- TRUE
        return(object)
    }
)





# Notes:
# The function is taken from old simulation code and left as is.
# inputs:
# PHENOTYPE_ORIGINAL -- ph2@original
# PHENOTYPE -- ph2@simulated
# N_PROT -- pr2@nProt
# pepCounts -- pr2@pepCounts
# features -- object@proteins@features
# PROPORTION_CHANGING_PROTEINS -- pr2@propChanging
# changingProtIdx -- pr2@changingProtIdx
# FDR_PEPTIDE -- pr2@fdrPeptide
# return:
# yGroupMeans -- pure mean peptide abundances fold changes
.generate_group_means <- function(PHENOTYPE_ORIGINAL,
                                  PHENOTYPE,
                                  N_PROT,
                                  pepCounts,
                                  features,
                                  PROPORTION_CHANGING_PROTEINS,
                                  changingProtIdx,
                                  meanFoldChange,
                                  FDR_PEPTIDE)
{
    # Generating group means
    #------------------------------------------------------------------------
    NG <- length(unique(PHENOTYPE))
    mean.prot.folds <- matrix(1, nrow=N_PROT, ncol=NG)
    # a more straight way to generate a percentage of TRUE idx
    # changingProtIdx <- sample(c(rep(TRUE,round(nProt*PROPORTIONCHANGING)),
    #                             rep(FALSE,nProt - round(nProt*PROPORTIONCHANGING))))
    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    simChanges <- t(replicate(N_PROT, 
                              c(1, rlnorm(NG-1, meanlog = 0, sdlog = meanFoldChange))))
    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    mean.prot.folds[changingProtIdx,] <- mean.prot.folds[changingProtIdx,] *
        simChanges[changingProtIdx,]
    # now expand from proteins to peptides
    # do I need to add technical peptide noise?
    mean.pep.folds <- mean.prot.folds[rep(seq_along(pepCounts), pepCounts),]
    #------------------------------------------------------------------------
    # post handling for reporting
    colnames(mean.prot.folds) <- unique(PHENOTYPE)
    mean.prot.folds.df <- as.data.frame(mean.prot.folds)
    for(i in setdiff(unique(PHENOTYPE_ORIGINAL), unique(PHENOTYPE))){
        mean.prot.folds.df[[i]] <- NA
    }
    mean.prot.folds <- as.matrix(mean.prot.folds.df[,unique(PHENOTYPE_ORIGINAL)])
    #------------------------------------------------------------------------
    
    
    
    # Now generate some false peptide identificaitons
    #------------------------------------------------------------------------
    # # weighing probabilities by protein representation
    # probs <- FDR_PEPTIDE*nPeps[rep(seq_along(nPeps), nPeps)]/
    #     mean(nPeps[rep(seq_along(nPeps), nPeps)])
    probs <- rep(FDR_PEPTIDE, nrow(features))
    idxOutlyingPeptides <- sapply(probs, rbinom, n=1,size=1) == 1
    # generate ALTERNATIVE group means
    mean.outlying.pep.folds <- matrix(1, nrow=nrow(features), ncol=NG)
    # switch to Bernoulli distro
    probs <- rep(PROPORTION_CHANGING_PROTEINS, nrow(features))
    changingPepIdx <- sapply(probs, rbinom, n=1, size=1) == 1
    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    simChanges <- t(replicate(nrow(features), c(1, rlnorm(NG-1, meanlog = 0, sdlog = 1))))
    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    mean.outlying.pep.folds[changingPepIdx,] <- 
        mean.outlying.pep.folds[changingPepIdx,] * simChanges[changingPepIdx,]
    # now expand from proteins to peptides
    # do I need to add technical peptide noise?
    mean.pep.folds[idxOutlyingPeptides,] <- mean.outlying.pep.folds[idxOutlyingPeptides,]
    #------------------------------------------------------------------------
    
    
    
    # The pre-final milestone.  Generating clean expression matrix.
    #------------------------------------------------------------------------
    # mm - model matrix
    mm <- model.matrix(~ PHENOTYPE + 0)
    # loops through peptides
    yGroupMeans <- sapply(seq_len(nrow(features)), function(i){
        multipliers <- t(t(mm) * mean.pep.folds[i,])
        rowSums(multipliers * 1)})
    yGroupMeans <- t(yGroupMeans)
    #------------------------------------------------------------------------
    return(yGroupMeans)
}



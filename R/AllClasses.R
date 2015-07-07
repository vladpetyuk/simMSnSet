

#' Phenotypes class.
#'
#' Details about Phenotypes class.
#' 
#' Class holding information about the samples
#' 
#' @slot n number of samples
#' @slot proportionOutlying proportion of outlying or 
#'          erroneusly assigned samples.  It is not uncommong that in
#'          large-scale experiments up to ~5\% samples appear as outliers.
#' @slot originalPhenotypes character vector with original or intended 
#'          phenotype (group) assignment
#' @slot simulatedPhenotypes character vector with phenotype assigments 
#'          after pretending that some of the samples are actually outliers
#' @slot sampleNames character vector with sample names
#' @slot meanFoldChange numeric value controlling the applitude of average fold
#'          of change between the conditions. Used as \code{sdlog} in
#'          \code{rlnorm} function.
#' @slot pipettingAccuracy a parameter for generating sample systematic biases
#' @slot sampleBiases normalization scaling factors
#' @slot simulated logical indicated if the object has been simulated
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
             originalPhenotypes="character",
             simulatedPhenotypes="character",
             meanFoldChange="numeric",
             pipettingAccuracy="numeric",
             sampleBiases="numeric",
             simulated="logical",
             sampleNames="character")
)

# * checking for code/documentation mismatches ... WARNING
# S4 class codoc mismatches from documentation object 'Phenotypes-class':
#     Slots for class 'Phenotypes'
# Code: n originalPhenotypes pipettingAccuracy proportionOutlying
# sampleBiases sampleNames simulated simulatedPhenotypes
# Docs: n proportionOutlying




#' Proteins class.
#'
#' Details about Proteins class.
#' 
#' Class holding the information about the proteins and peptides
#' 
#' @slot nProt number of proteins
#' @slot fdrPeptide FDR of peptide identifications
#' @slot fdrMatching FDR of peptide to spectrum or LC-MS feature matching
#' @slot propChanging proportion of proteins that are changed between the 
#'          different phenotype groups.  At this point it is independent of
#'          phenotype groups.  It shoud, however, depend on phenotype.
#' @slot changingProtIdx index indicating which proteins are changing
#' @slot nPep number of peptides
#' @slot lambda parameter of Poisson distribution controlling the number of 
#'          peptides per protein during the simulation
#' @slot meanIoinizationIntensity a numeric value of mean peptide intensity
#' @slot noise measurement noise. So far one value for all. May be in the
#'          future convert to signal/noise ratio.
#' @slot pepCounts integer vector with number of peptides per each protein
#' @slot protNames character vector with generated protein names
#' @slot pepNames character vector with generacted peptide names
#' @slot ionizationEff numeric vector with simulated ionization efficiencies 
#'          for the peptides
#' @slot missThreshAbsolute soft threshold in absolute intensity. 
#'          It has precedence over missThreshQuantile.
#' @slot missThreshQuantile quantile to determine soft threshold value.
#' @slot missThreshSharpness parameter controlling the sharpness of the
#'           thresholding step function. The suggested interval for the
#'           sharpness values is [2,100]. At 2 it approximates linear 
#'           (not necessarily the best setting). At 100 it is almost a 
#'           step-function. The suggested values are >5.
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
             fdrPeptide="numeric",
             fdrMatching="numeric",
             propChanging="numeric",
             changingProtIdx="logical",
             nPep="integer",
             lambda="numeric",
             meanIoinizationIntensity="numeric",
             noise="numeric",
             pepCounts="integer",
             protNames="character",
             pepNames="character",
             ionizationEff="numeric",
             missThreshAbsolute="numeric",
             missThreshQuantile="numeric",
             missThreshSharpness="numeric",
             features="data.frame")
)






#' Experiment class.
#'
#' Details about Experiment class.
#' 
#' The central class of simMSnSet package
#' 
#' @slot phenotypes Phenotypes and Samples
#' @slot proteins Peptides and Proteins
#' @slot simFinal simulated final intensities
#' @slot simNoise simulated noise of intensity measurements
#' @slot simMeans simulated actual fold of changes for each peptide/sample
#' @slot simSampleBiases simulated systematic biases due to pipetting errors
#' @slot simOutliers indexes of simulated outlying samples
#' @slot simMissing index of missing measurements
#' @slot simulated A logical flag indicating if the experiment has been simulated or not
#' 
#' @examples
#' # 10 proteins measured across 6 samples (3 controls + 3 cases)
#' pr <- Proteins(nProt=10)
#' ph <- Phenotypes(originalPhenotypes=c('Control', 'Control', 'Control', 
#'                             'Case', 'Case', 'Case'))
#' ex <- Experiment(pr, ph)
#' ex <- simulate(ex)
#' image(ex)
#'  
#' @exportClass Experiment
#' @family Experiment-class
#' 
setClass(Class="Experiment",
         representation(
             phenotypes="Phenotypes",
             proteins="Proteins",
             simFinal="matrix",
             simNoise="matrix",
             simMeans="matrix",
             simSampleBiases="matrix",
             simOutliers="matrix",
             simMissing="matrix",
             simulated="logical")
)

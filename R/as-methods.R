
# ' As("Foo", "SpatialPointsDataFrame")
# '
# ' @name as
# ' @family Foo
# '
# ' @importClassesFrom sp SpatialPointsDataFrame

#' As("Experiment", "MSnSet")
#' @name as
#' @family Experiment-class
# ' @describeIn Experiment Coerce simulated experiment into 
# '              \link[MSnbase]{MSnSet-class} object
#' @importFrom MSnbase MSnSet experimentData<-
#' @importClassesFrom MSnbase MIAPE
# ' @exportMethod as
#' @param from \link{Experiment-class} object

# ' @importFrom MSnbase MSnSet experimentData<-
# ' @importClassesFrom MSnbase MIAPE
# ' @describeIn Experiment Coerce simulated experiment into 
# '              \link[MSnbase]{MSnSet-class} object
# ' @exportMethod as
# '

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
          pd <- data.frame(phenotype=from@phenotypes@originalPhenotypes, 
                           row.names=from@phenotypes@sampleNames)
          obj <- MSnSet(exprs=from@simFinal, 
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


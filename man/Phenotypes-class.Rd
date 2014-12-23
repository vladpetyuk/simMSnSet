% Generated by roxygen2 (4.0.2): do not edit by hand
\docType{class}
\name{Phenotypes-class}
\alias{Phenotypes-class}
\alias{show,Phenotypes-method}
\alias{simulate,Phenotypes-method}
\title{Phenotypes class.}
\usage{
\S4method{simulate}{Phenotypes}(object, seed = NULL)

\S4method{show}{Phenotypes}(object)
}
\arguments{
\item{seed}{an integer for \code{set.seed}}
}
\description{
Details about Phenotypes class.
}
\details{
Class holding information about the samples
}
\section{Methods (by generic)}{
\itemize{
\item \code{simulate}: Selects which samples are going to be outliers given
the value of the \code{proportionOutlying} slot.

\item \code{show}: Show method for a Phenotype object
}}
\section{Slots}{

\describe{
\item{\code{n}}{number of samples}

\item{\code{proportionOutlying}}{proportion of outlying or
erroneusly assigned samples.  It is not uncommong that in
large-scale experiments up to ~5% samples appear as outliers.}

\item{\code{original}}{character vector with original or intended
phenotype (group) assignment}

\item{\code{simulated}}{character vector with phenotype assigments after pretending
that some of the samples are actually outliers}

\item{\code{sampleNames}}{character vector with sample names}
}}
\examples{
# 10 samples (5 controls + 5 cases)
# Note, the number of outlying samples is somewhat extreme - 50\%.
# This is for demonstration purpose only.
ph <- Phenotypes(original=c('Ctrl', 'Ctrl', 'Ctrl', 'Ctrl', 'Ctrl',
                            'Case', 'Case', 'Case', 'Case', 'Case'),
                 proportionOutlying=0.5)
ph <- simulate(ph, seed=123)
show(ph)
}

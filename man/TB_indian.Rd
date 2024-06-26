% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{TB_indian}
\alias{TB_indian}
\title{An example TB dataset with Indian population data.}
\format{
SummarizedExperiment
}
\usage{
TB_indian
}
\description{
An example dataset containing the gene expression and metadata in a
SummarizedExperiment object for an Indian population. Active TB contamination
of the 44 subjects is denoted for each as a "1"(active) or "0"
(latent/not present), and can be accessed via \code{TB_indian$label}. The
SummarizedExperiment object contains 2 assays (counts and log(counts)),
and the column names give the unique subject identification number along
with the subject's gender.
}
\details{
This dataset was published as part of a study to assess performance of
published TB signatures in a South Indian population (Leong et. al 2018).
RNA sequencing was performed on whole blood PAX gene samples collected
from 28 TB patients and 16 latent TB infected (LTBI) subjects enrolled
as part of an ongoing household contact study. Whole blood RNA-Seq analysis
was performed on all 44 participants.
}
\examples{
data("TB_indian")
}
\references{
Leong S., Zhao Y., et. al. (2018). Existing blood transcriptional classifiers
accurately discriminate active tuberculosis from latent infection in
individuals from south India. \emph{Tuberculosis} \strong{109}, 41-51.
doi: \url{https://doi.org/10.1016/j.tube.2018.01.00210.1016/j.tube.2018.01.002}.
}
\keyword{datasets}

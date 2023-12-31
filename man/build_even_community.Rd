% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/build_even_community.R
\name{build_even_community}
\alias{build_even_community}
\title{Build even community matrix}
\usage{
build_even_community(
  n.taxa = 5,
  n.samples = 5,
  n.reads = 5000,
  taxa.dist = "normal",
  taxa.sd = 0
)
}
\arguments{
\item{n.taxa}{Positive numeric vector of length 1. The number of taxa (columns) to generate. Default = 5.}

\item{n.samples}{Positive numeric vector of length 1. The number of samples (rows) to generate. Default = 5.}

\item{n.reads}{Positive numeric vector of length 1. The number of observations within each sample. This is analogous to per-sample sequencing depth.}

\item{taxa.dist}{Currently, only the 'normal' distribution is accepted.}

\item{taxa.sd}{The standard deviation for the randomization function that generates taxa abundances. Higher values generate more zeros for each taxon.}
}
\value{
Taxon abundance matrix with samples as rows and taxa as columns. class='matrix'
}
\description{
Builds a community matrix with samples as rows and taxa as columns.
Taxa abundances are drawn from a normal distribution for each sample (this is currently the only distribution supported).
}
\examples{
comm <- build_even_community(n.taxa = 100,n.samples = 44,n.reads = 3000, taxa.sd = 30)

}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{goodGenes}
\alias{goodGenes}
\alias{goodGenes,GeneAberExpr-method}
\title{goodGenes processes the Ice cream HMM example}
\usage{
goodGenes(data, max)
}
\value{
bool vector of rows with less than 40% over max columns.
}
\description{
check in the sparse matrix if at least 40% of the row values are less than max
}

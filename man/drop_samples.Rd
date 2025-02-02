% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_manip.R
\name{drop_samples}
\alias{drop_samples}
\title{Drop samples}
\usage{
snpdata = drop_samples(snpdata, samples.to.be.dropped)
}
\arguments{
\item{snpdata}{a SNPdata object}

\item{samples.to.be.dropped}{a vector of samples to be dropped}
}
\value{
a SNPdata object where the specified samples have been removed
}
\description{
remove a set of samples from the SNPdata object
}

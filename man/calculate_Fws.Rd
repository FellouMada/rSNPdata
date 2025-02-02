% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_manip.R
\name{calculate_Fws}
\alias{calculate_Fws}
\title{Calculate the complexity of infection in every sample (Fws)}
\usage{
snpdata =  calculate_Fws(snpdata)
}
\arguments{
\item{snpdata}{a SNPdata object}
}
\value{
a SNPdata object with an additional column in the meta table
\enumerate{
  \item Fws: within host genetic diversity value
  \item COI: the complexity of infection: 1 for Fws>0.95, 2 for Fws<=0.95
}
}
\description{
Fws is the within host genetic diversity
}

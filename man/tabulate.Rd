% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tabulation.R
\name{tabulate}
\alias{tabulate}
\title{Tabulate MCMC Output}
\usage{
tabulate(results, burnin = 0.1)
}
\arguments{
\item{results}{MCMC output as generated by run_mcmc.}

\item{burnin}{We discard the first (burnin * 100) percent of the results}
}
\value{
A data frame.
}
\description{
This function creates a table of the posterior probabilities associated with each possible transmission link.
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/recbase.R
\name{recbase}
\alias{recbase}
\title{Baseline estimation and plot}
\usage{
recbase(fit, formula, recdata, tau = NULL, np = 40)
}
\arguments{
\item{fit}{the output from the rectime() function.}

\item{formula}{the same formula used in the rectime() function.}

\item{recdata}{the rectime data object.}

\item{tau}{the pre-specified time point such that the recurrent event process could be
potentially observed beyond tau. Usually it is the length of the study. If left unspecified,
the default is the maximum of the time variable.}

\item{np}{number of time points where the baseline rate function is evaluated. The time
points are equally spaced between 0 to tau.}
}
\description{
Estimate baseline rates at equally spaced time points using the same method selected for
model estimation.
}

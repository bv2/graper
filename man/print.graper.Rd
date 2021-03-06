% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/print.graper.R
\name{print.graper}
\alias{print.graper}
\title{Print a graper object}
\usage{
\method{print}{graper}(x, ...)
}
\arguments{
\item{x}{fitted graper model as obtained from  \code{\link{graper}}}

\item{...}{additional print arguments}
}
\value{
Print output.
}
\description{
Function to print a fitted graper model.
}
\examples{
# create data
dat <- makeExampleData()
# fit the graper model
fit <- graper(dat$X, dat$y, dat$annot)
# print a summary of the fitted model
print(fit)
}

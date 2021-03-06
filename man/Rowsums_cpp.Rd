% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{Rowsums_cpp}
\alias{Rowsums_cpp}
\title{Rowsums_cpp}
\usage{
Rowsums_cpp(P, SK)
}
\arguments{
\item{P}{is a matrix}

\item{SK}{is a data sketch}
}
\value{
a vector of length \eqn{2m}
}
\description{
Implementation of a matrix rows' summation and a sketch vector subtraction
}
\details{
The sums over the rows of the matrix \eqn{P \in R^{2m \time n}} is a sketch of the decomposition,
\emph{i.e.} \eqn{SK(D\cdot A)}. It is defined as \eqn{P = [cos(W\cdot D \cdot A);sin(W\cdot D \cdot A)]}, where \eqn{W} is a frequency matrix, \eqn{D} is
a dictionary and \eqn{A} is a code matrix. This function computes a vector of the difference between a
data sketch and the decomposition sketch: \eqn{SK - SK(D\cdot A)}.
}
\seealso{
\code{\link{Gradient_D_cpp_parallel}}
}
\keyword{internal}

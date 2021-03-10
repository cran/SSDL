#' @title SSDL-package
#' 
#' @description R package \code{SSDL} implements the Sketched Stochastic Dictionary Learning method that builds a dictionary from 
#' large-scale data collection by operating on a compressed data version referred to as a data sketch.
#' The [chickn](https://CRAN.R-project.org/package=chickn) package is used to carry out the data compression. 
#' SSDL package is designed to handle voluminous data encoded as a matrix, which cannot be loaded in memory. 
#' To do this, SSDL package relies on the Filebacked Big Matrix class of [bigstatsr](https://github.com/privefl/bigstatsr) package, 
#' which allows to access and manipulate matrix-like data stored in files on disk. 
#' @seealso \code{\link{CDL}}
#' @docType package
#' @author Olga Permiakova, Thomas Burger
#' @import pracma
#' @importFrom Rcpp evalCpp
#' @importFrom Rdpack reprompt
#' @useDynLib SSDL, .registration=TRUE
#' @name SSDL
NULL  
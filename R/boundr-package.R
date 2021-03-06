#' The 'boundr' package.
#'
#' @description This R package is used to build and analyze structural graphical models as proposed by Pearl (2009, chapter 8). This approach allows for bounded estimation of counterfactuals in an experiment where (i) all variables are discrete and (ii) not all variables are directly manipulated by the randomized assignment mechanism. This is done by modeling the correlated latent background variables determining all observable variables. Inference is conducted by constructing a Bayesian statistical model.
#'
#' @docType package
#' @name boundr-package
#' @aliases boundr
#' @useDynLib boundr, .registration = TRUE
#' @import methods
#' @import Rcpp
#' @import rlang
#' @import purrr
#' @import dplyr
#' @import stringr
#' @import forcats
#' @import tidyr
#' @import tibble
#' @importFrom magrittr %>% %<>%
#' @importFrom rstan sampling vb
#' @importFrom loo loo loo.array
#'
#' @references
#' Chickering, D. M., & Pearl, J. (1996). A clinician’s tool for analyzing non-compliance. Proceedings of the National Conference on Artificial Intelligence, 2, 1269–1276.
#' Pearl, J. (2009). Causality: Models, Reasoning and Inference (Second Ed). Cambridge: Cambridge University Press.
#' Stan Development Team (2019). RStan: the R interface to Stan. R package version 2.19.2. https://mc-stan.org
#'
NULL

#' evclass: A package for evidential classification
#'
#' The evclass package currently contains functions for two evidential classifiers: the evidential
#' K-nearest neighbor (EK-NN) rule (Denoeux, 1995; Zouhal and Denoeux, 1998) and the evidential
#' neural network (Denoeux, 2000). In contrast with classical statistical classifiers, evidential
#' classifier quantify the uncertainty of the classification using Dempster-Shafer mass functions.
#'
#' The main functions are: \code{\link{EkNNinit}}, \code{\link{EkNNfit}} and \code{\link{EkNNval}}
#' for the initialization, training and evaluation of the EK-NN classifier,
#' \code{\link{proDSinit}}, \code{\link{proDSfit}} and \code{\link{proDSval}} for the
#' evidential neural network classifier, and \code{\link{decision}} for decision-making.
#'
#' @docType package
#' @name evclass
#'
#' @seealso \code{\link{EkNNinit}}, \code{\link{EkNNfit}},
#'\code{\link{EkNNval}}, \code{\link{proDSinit}}, \code{\link{proDSfit}}, \code{\link{proDSval}}.
#'
#' @references
#'T. Denoeux. A k-nearest neighbor classification rule based on Dempster-Shafer
#'theory. IEEE Transactions on Systems, Man and Cybernetics, 25(05):804--813, 1995.
#'
#'T. Denoeux. Analysis of evidence-theoretic decision rules for pattern
#'classification. Pattern Recognition, 30(7):1095--1107, 1997.
#'
#'T. Denoeux. A neural network classifier based on Dempster-Shafer theory.
#'IEEE Trans. on Systems, Man and Cybernetics A, 30(2):131--150, 2000.
#'
#' L. M. Zouhal and T. Denoeux. An evidence-theoretic k-NN rule with parameter
#' optimization. IEEE Transactions on Systems, Man and Cybernetics Part C,
#' 28(2):263--271,1998.
#'
#'Available from \url{https://www.hds.utc.fr/~tdenoeux}.
#'
NULL

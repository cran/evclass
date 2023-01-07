#' evclass: A package for evidential classification
#'
#' The evclass package currently contains functions for three evidential classifiers: the evidential
#' K-nearest neighbor (EK-NN) rule (Denoeux, 1995; Zouhal and Denoeux, 1998), the evidential
#' neural network (Denoeux, 2000) and the RBF classifier with weight-of-evidence interpretation
#' (Denoeux, 2019; Huang et al., 2022), as well as methods to compute output mass functions from
#' trained logistic regression or multilayer classifiers as described in (Denoeux, 2019). In contrast
#' with classical statistical classifiers, evidential classifiers quantify the uncertainty of the
#' classification using Dempster-Shafer mass functions.
#'
#' The main functions are: \code{\link{EkNNinit}}, \code{\link{EkNNfit}} and \code{\link{EkNNval}}
#' for the initialization, training and evaluation of the EK-NN classifier;
#' \code{\link{proDSinit}}, \code{\link{proDSfit}} and \code{\link{proDSval}} for the
#' evidential neural network classifier; \code{\link{decision}} for decision-making;
#' \code{\link{RBFinit}}, \code{\link{RBFfit}} and \code{\link{RBFval}} for the RBF classifier;
#' \code{\link{calcAB}} and \code{\link{calcm}} for computing output mass functions from trained
#' logistic regression or multilayer classifiers.
#'
#' @docType package
#' @name evclass
#'
#' @seealso \code{\link{EkNNinit}}, \code{\link{EkNNfit}},
#'\code{\link{EkNNval}}, \code{\link{proDSinit}}, \code{\link{proDSfit}}, \code{\link{proDSval}},
#'\code{\link{RBFinit}}, \code{\link{RBFfit}} and \code{\link{RBFval}}, \code{\link{decision}},
#'\code{\link{calcAB}}, \code{\link{calcm}}.
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
#' T. Denoeux. Logistic Regression, Neural Networks and Dempster-Shafer Theory: a New Perspective.
#'Knowledge-Based Systems, Vol. 176, Pages 54–67, 2019.
#'
#' L., S. Ruan, P. Decazes and T. Denoeux. Lymphoma segmentation from 3D PET-CT images using a
#'  deep evidential network. International Journal of Approximate Reasoning, Vol. 149, Pages 39-60,
#'  2022.
#'
NULL

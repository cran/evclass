#' Glass dataset
#'
#' This data set contains the description of 214 fragments of glass originally
#' collected for a study in the context of criminal investigation. Each fragment has a measured
#' reflectivity index and chemical composition (weight percent of Na, Mg, Al, Si, K, Ca, Ba and Fe).
#' As suggested by Ripley (1994), 29 instances were discarded, and the remaining 185
#' were re-grouped in four classes: window float glass (70), window non-float glass (76), vehicle window
#' glass (17) and other (22). The data set was split randomly in a training set of size 89 and a test
#' set of size 96.
#'
#' @docType data
#'
#' @usage data(glass)
#'
#' @format A list with two elements:
#' \describe{
#' \item{x}{The 185 x 9 object-attribute matrix.}
#' \item{y}{A 185-vector containing the class labels.}
#' }
#'
#' @keywords datasets
#'
#' @references
#'
#' P. M.  Murphy and D. W. Aha.  UCI Reposition of machine learning databases.
#' [Machine readable data repository]. University of California, Departement of
#' Information and Computer Science, Irvine, CA.
#'
#' B.D.Ripley, Flexible nonlinear approaches to classification, in "From Statistics
#' to Neural Networks", V. Cherkassly, J. H. Friedman, and H. Wechsler, Eds.,
#' Berlin, Germany: Springer-Verlag, 1994, pp. 105--126.
#'
#' T. Denoeux. A neural network classifier based on Dempster-Shafer theory.
#'IEEE Trans. on Systems, Man and Cybernetics A, 30(2):131--150, 2000.
#'
#'
#' @examples
#' data(glass)
#' table(glass$y)
"glass"

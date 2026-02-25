#' @keywords internal
#'
#' @useDynLib scatR, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom methods as
#' @importFrom Matrix sparseMatrix Diagonal crossprod
#' @importFrom stats model.matrix model.frame terms formula dist median na.pass
#' @importFrom graphics plot par lines polygon
#' @importFrom grDevices adjustcolor
"_PACKAGE"

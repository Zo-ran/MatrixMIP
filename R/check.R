#' @noRd
check_input <- function(Initial_X, A, B, K, Tr, Scalar, Bscalar, Gpp, ...) {
  plist <- list(...)

  #check matrix
  if (!is.matrix(Scalar))
    stop("Scalar must be a matrix")
  if (!is.matrix(Bscalar))
    stop("Bscalar must be a matrix")
  check_matrix <- function(X, name) {
    if (!is.matrix(X)) {
      stop(paste(name, "must be a matrix"))
    } else if (dim(X)[1] != dim(X)[2]) {
      stop(paste("the matrix", name, "must be square"))
    }
  }
  check_matrix(A, "A")
  check_matrix(K, "K")
  check_matrix(Tr, "Tr")

  #check dimension
  if (!is.vector(B))
    stop("B must be a vector")
  if (!is.vector(Initial_X))
    stop("Initial_X must be a vector")
  if (!is.vector(Gpp))
    stop("Gpp must be a vector")
  if (dim(A)[2] != dim(K)[2] || dim(A)[2] != dim(Tr)[2] || dim(A)[2] != dim(Scalar)[2] || dim(A)[2] != dim(Bscalar)[2])
    stop("A, B, K, Tr, Scalar, Bscalar must have same number of columns")
  if (dim(A)[2] != length(Initial_X) || dim(A)[2] != length(B))
    stop("the length of vector Initial_X, B must be the same as the number of columns in matrix A")


  if (length(plist) > 0) {
    if (plist[[1]] < 0)
      stop("iterations must be a positive integer")

    if (plist[[2]] != "month" && plist[[2]] != "day")
      stop("timestep must be \"month\" or \"day\"")

    if (plist[[2]] == "month") {
      if (dim(Scalar)[1] < 12)
        stop("not enough samples for sasu; if timestep == \"month\", the number of rows of Scalar must be at least 12")
      if (dim(Bscalar)[1] < 12)
        stop("not enough samples for sasu; if timestep == \"month\", the number of rows of Bscalar must be at least 12")
      if (length(Gpp) < 12)
        stop("not enough samples for sasu; if timestep == \"month\", the length of Gpp must be at least 12")
    } else {
      if (dim(Scalar)[1] < 365)
        stop("not enough samples for sasu; if timestep == \"day\", the number of rows of Scalar must be at least 365")
      if (dim(Bscalar)[1] < 365)
        stop("not enough samples for sasu; if timestep == \"day\", the number of rows of Bscalar must be at least 365")
      if (length(Gpp) < 365)
        stop("not enough samples for sasu; if timestep == \"day\", the length of Gpp must be at least 365")
    }
  }
}

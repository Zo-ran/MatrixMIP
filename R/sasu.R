#' @title Semi-Analytical Spin-Up

#' @description  A semi-analytical solution to accelerate spin-up

#' @param iterations number of iterations of the model

#' @param timestep a character representing the time interval between each sample; timestep = "month" or "day"

#' @param Initial_X a vector of the initial size of each C pool

#' @param A a matrix indicating coefficients of transfer among C pools

#' @param B a vector indicating indicating plant C allocation; Only plant C pools will not have a value of 0

#' @param K a diagonal matrix indicating C leaving individual pools through mortality or decomposition

#' @param Tr a matrix indicating vertical transfer coefficients

#' @param Scalar a matrix indicating environmental scalars, which is used to modify the matrix K

#' @param Bscalar a matrix indicating C input through photosynthesis, which is used to modify the matrix B

#' @param Gpp a vector of daily or monthly gross primary production(Gpp)

#' @return a vector of the initial size of each C pool

#' @examples
#' pool_n <- 8   #the number of pools
#'
#' Initial_X <- rep(0, pool_n)
#'
#' A <- diag(pool_n)
#' A[4, 1] <- -1, A[4, 2] <- -1, A[4, 3] <- -0.001
#' A[5, 3] <- -0.999
#' A[6, 4] <- -0.45, A[6, 5] <- -0.275, A[6, 7] <- -0.42, A[6, 8] <- -0.45
#' A[7, 5] <- -0.275, A[7, 6] <- -0.296
#' A[8, 6] <- -0.004, A[8, 7] <- -0.03
#'
#' B <- c(0.16875, 0.16875, 0.1125, 0, 0, 0, 0, 0)
#'
#' K <- diag(c(0.00274, 0.00274, 0.0000684, 0.00913, 0.000472, 0.00684, 0.0000548, 0.00000137))
#'
#' Tr <- matrix(0, pool_n, pool_n)
#'
#'

#' @export

sasu <- function(iterations, timestep, Initial_X, A, B, K, Tr, Scalar, Bscalar, Gpp) {
  #### input check
  #check timestep
  if (timestep != "month" && timestep != "day") {
    stop("timestep must be \"month\" or \"day\"")
  }

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
    stop("A, B, K, Tr, Scalar, Bscalar must have same numsber of columns")
  if (dim(A)[2] != length(Initial_X) || dim(A)[2] != length(B))
    stop("the length of vector Initial_X, B must be the same as sthe number of columns in matrix A")
  if (dim(Scalar)[1] != dim(Bscalar)[1] || dim(Scalar)[1] != length(Gpp))
    stop("the length of vector Gpp and the number of columns in Scalar, Bscalar should be the same")


  #### some variables
  AVG_B <- B / sum(B)
  date_length <- dim(Scalar)[1] # the length of date
  pool_n <- dim(A)[2]           # the number of pool
  X <- matrix(0, nrow=date_length, ncol=pool_n) # the return value
  datelength4sasu <- if (timestep == "day") 365 else 12
  I4simu <- Gpp * sum(B) # calculate NPP (I4simu/I4sasu) from Gpp, NPP as model input
  I4sasu <- Gpp[1:datelength4sasu] * sum(B)


  #### spin-up using SASU
  X[1, ] <- Initial_X
  scal.mean <- colMeans(Scalar[1:datelength4sasu, ])
  scal.sasu <- diag(x = 0, pool_n, pool_n)
  diag(scal.sasu) <- scal.mean
  matrixAK_sasu <- A %*% (scal.sasu %*% K) + Tr
  matrixBI_sasu <- colMeans(Bscalar1[1:datelength4sasu, ]) * mat.B * mean(I4sasu)


  # temporal C storage capacity
  X[1, ] <- solve(matrixAK_sasu) %*% matrixBI_sasu

  for(m in 1:iterations) {
    for(k in 1:datelength4sasu) {
      scal.temp <- diag(x = 0, pn, pn); diag(scal.temp) <- as.numeric(Scalar[k, ])
      matrixAK <- mat.A %*% (scal.temp %*% K) + Tr
      if (k == 1)
        l <- k
      else
        l <- k - 1
      X[k, ] <- diag(pn) %*% X[l, ] + (Bscalar[k, ] * AVG_B) * I4simu[k] - matrixAK %*% X[l, ]
      if (k == datelength4sasu)
        X[1, ] = X[k, ]
    }
  }
  return(X[1, ])
}

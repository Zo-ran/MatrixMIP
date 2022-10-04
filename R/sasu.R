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
#' pool_n <- 15   #the number of pools
#'
#' Initial_X <- rep(0, pool_n)
#'
#' A <- diag(pool_n)
#' A[6, 1] <- -0.75
#' A[7, 1] <- -0.25
#' A[8, 2] <- -0.75
#' A[9, 2] <- -0.25
#' A[10, 3] <- -1
#' A[11, 4] <- -1
#' A[12, 5] <- -1
#' A[13, 6] <- -0.4125; A[13, 7] <- -0.45; A[13, 8] <- -0.3375; A[13, 9] <- -0.45; A[13, 10:12] <- -0.125; A[13, 14] <- -0.42
#' A[14, 6] <- -0.175; A[14, 8] <- -0.175; A[14, 10:12] <- -0.375; A[14, 13] <- -0.1664; A[14, 15] <- -0.45
#' A[15, 13] <- -0.004; A[15, 14] <- -0.03
#'
#' B <- c(0.1725, 0.15, 0.0225, 0.14, 0.015, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
#'
#' K <- diag(c(0.09, 0.05, 0.01, 0.0008, 0.001, 0.153857, 1.2, 0.190296, 1.5, 0.0555556, 0.166667, 0.138889, 0.5865, 0.0162857, 0.000557143))
#'
#' Tr <- matrix(0, pool_n, pool_n)
#'
#' datelength4sasu <- 12
#' Scalar <- matrix(0, nrow = datelength4sasu, ncol = pool_n)
#' Scalar[, 1:5] <- 1
#' Scalar[, 6:15] <- c(0.001125247, 0.026934919, 0.00386343, 0.036797823, 0.056063114, 0.23379085, 0.168964848, 0.259775597, 0.019196622, 0.061492205, 0.013498124, 0.004996017)
#'
#' Bscalar <- matrix(0, nrow = datelength4sasu, ncol = pool_n)
#' Bscalar[ , ] <- 1
#'
#' Gpp <- c(1.5074, 12.8787, 33.9259, 78.7088, 142.6238, 171.0413, 194.6502, 203.5172, 145.2821, 96.9997, 27.6457, 7.9289)
#'
#' initial <- sasu(100, "month", Initial_X, A, B, K, Tr, Scalar, Bscalar, Gpp)
#' initial

#' @export

sasu <- function(iterations, timestep, Initial_X, A, B, K, Tr, Scalar, Bscalar, Gpp) {
  #### input check
  check_input(Initial_X, A, B, K, Tr, Scalar, Bscalar, Gpp, iterations, timestep)

  #### some variables
  AVG_B <- B / sum(B)
  date_length <- dim(Scalar)[1] # the length of date
  pool_n <- dim(A)[2]           # the number of pools
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
  matrixBI_sasu <- colMeans(Bscalar[1:datelength4sasu, ]) * AVG_B * mean(I4sasu)


  # temporal C storage capacity
  X[1, ] <- solve(matrixAK_sasu) %*% matrixBI_sasu

  for (m in 1:iterations) {
    for (k in 1:datelength4sasu) {
      scal.temp <- diag(x = 0, pool_n, pool_n); diag(scal.temp) <- as.numeric(Scalar[k, ])
      matrixAK <- A %*% (scal.temp %*% K) + Tr
      if (k == 1)
        l <- k
      else
        l <- k - 1
      X[k, ] <- diag(pool_n) %*% X[l, ] + (Bscalar[k, ] * AVG_B) * I4simu[k] - matrixAK %*% X[l, ]
      if (k == datelength4sasu)
        X[1, ] = X[k, ]
    }
  }
  return(X[1, ])
}

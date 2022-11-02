#' @title Forward Run

#' @description predict C pool sizes and C fluxes

#' @param Initial_X a vector of the initial size of each C pool, which can be obtained by "sasu"

#' @param A a matrix indicating coefficients of transfer among C pools

#' @param B a vector indicating indicating plant C allocation; Only plant C pools will not have a value of 0

#' @param K a diagonal matrix indicating C leaving individual pools through mortality or decomposition

#' @param Tr a matrix indicating vertical transfer coefficients

#' @param Scalar a matrix indicating environmental scalars, which is used to modify the matrix K

#' @param Bscalar a matrix indicating C input through photosynthesis, which is used to modify the matrix B

#' @param Gpp a vector of daily or monthly gross primary production(Gpp)

#' @param soil_depth a vector of depth of soil

#' @param OP order of plant C pools; for example, if C pools 1, 3 and 4 are plant C pools, 'OP' should be c(1, 3, 4)

#' @param OL order of litter C pools; see 'OP' for details

#' @param OS order of soil C pools; see 'OP' for details

#' @return a list of C pool size, net ecosystem exchange of plant C pools, net ecosystem exchange of litter C pools, net ecosystem exchange of soil C pools, carbon input to plant, litter respiration, carbon input to litter, soil respiration and carbon input to soil
#' You can use the return value as follows:
#' result <- forward_run(...); result$C_pool_size; result$NEC_plant; result$NEC_litter; result$NEC_soil; result$Ip; result$Rlit; result$Ilit; result$Rsol; result$Isol           

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
#' 
#' soil_depth <- rep(0, pool_n)
#' OP <- 1:5; OL <- 6:12; OS <- 13:15
#' 
#' res <- forward_run(initial, A, B, K, Tr, Scalar, Bscalar, Gpp, soil_depth, OP, OL, OS)
#' res$C_pool_size

#' @export

forward_run <- function(Initial_X, A, B, K, Tr, Scalar, Bscalar, Gpp, soil_depth, OP = NULL, OL  = NULL, OS  = NULL) {
  #### input check
  check_input(Initial_X, A, B, K, Tr, Scalar, Bscalar, Gpp)

  pool_n <- dim(A)[2]    #the number of pool
  date_length <- min(dim(Scalar)[1], dim(Bscalar)[1], dim(Gpp)[1])
  AVG_B <- B / sum(B)
  I4simu <- Gpp * sum(B) # calculate NPP (I4simu) from Gpp, NPP as model input

  check_order(OP, "'OP'", pool_n)
  check_order(OL, "'OL'", pool_n)
  check_order(OS, "'OS'", pool_n)
  if (!is.vector(soil_depth))
    stop(paste("'soil_depth'", "must be a vector"))

  X <- matrix(0, nrow = date_length, ncol = pool_n)            # C pool size
  X[1, ] <- Initial_X
  Tn <- matrix(0, nrow = date_length, ncol = pool_n)           # Tn, residence time of individual pool in network
  Xc <- matrix(0, nrow = date_length, ncol = pool_n)           # Xc, ecosystem storage capacity
  Xp <- matrix(0, nrow = date_length, ncol = pool_n)           # Xp, ecosystem storage potential
  NEC <- matrix(0, nrow = date_length, ncol = pool_n)          # NEC, net ecosystem exchange


  Rlit <- if (is.null(OL)) NULL else rep(0, date_length)
  Rsol <- if (is.null(OS)) NULL else rep(0, date_length)
  Ilit <- if (is.null(OL)) NULL else rep(0, date_length)
  Isol <- if (is.null(OS)) NULL else rep(0, date_length)
  Ip <- if (is.null(OP)) NULL else rep(0, date_length)

  for (f in 1:date_length) {
    scal.temp <- Scalar[f, ]
    matrixAK <- A %*% (diag(scal.temp) %*% K) + Tr
    g <- if (f == 1) f else f - 1

    # update of carbon pool sizes
    X[f, ] <- diag(pool_n) %*% X[g, ] + (Bscalar[f, ] * AVG_B) * I4simu[f] - matrixAK %*% X[g, ]

    diag(matrixAK)[diag(matrixAK) < 1e-15] <- 1e-15
    Tch <- solve(matrixAK)      # Tch, chasing time
    Tn[f, ] <- Tch %*% (Bscalar[f, ] * AVG_B)
    Xc[f, ] <- Tn[f, ] * I4simu[f]
    Xp[f, ] <- Xc[f, ] - X[f, ]
    NEC[f, ] <- matrixAK %*% Xp[f, ]

    # calculate C fluxes
    cal_respiration <- function(order_vector) {
      if (length(order_vector) >= 2){
        return (sum(colSums(A[, order_vector]) * scal.temp[order_vector] * diag(K)[order_vector] * X[f, order_vector] * soil_depth[order_vector]))
      } else {
        return (sum(sum(A[, order_vector]) * scal.temp[order_vector] * diag(K)[order_vector] * X[f, order_vector] * soil_depth[order_vector]))
      }
    }

    cal_input <- function(order_vector) {
      return (sum(diag(A)[order_vector] * scal.temp[order_vector] * diag(K)[order_vector] * X[f, order_vector] * soil_depth[order_vector]))
    }

    if (!is.null(OS)) {
      Rsol[f] <- cal_respiration(OS)
      Isol[f] <- cal_input(OS)
    }
    if (!is.null(OL)) {
      Rlit[f] <- cal_respiration(OL)
      Ilit[f] <- cal_input(OL)
    }
    if (!is.null(OP))
      Ip[f] <- sum(Bscalar[f, OP] * AVG_B[OP]) * I4simu[f]
  }

  NEC <- NEC %*% diag(soil_depth)
  NEC_plant <- if (length(OP) >= 2) rowSums(NEC[, OP]) else NEC[, OP]
  NEC_litter <- if (length(OL) >= 2) rowSums(NEC[, OL]) else NEC[, OL]
  NEC_soil <- if (length(OS) >= 2) rowSums(NEC[, OS]) else NEC[, OS]

  res <- list(X, NEC_plant, NEC_litter, NEC_soil, Ip, Rlit, Ilit, Rsol, Isol)
  names(res) <- c("C_pool_size", "NEC_plant", "NEC_litter", "NEC_soil", "Ip", "Rlit", "Ilit", "Rsol", "Isol")
  return (res)
}

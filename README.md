# MatrixMIP

## Introduction

The aim of matrixMIP(Model Intercomparison Project) is to help you build land C cycle models in a in a unified matrix form. 

This model is based on the fact that although hundreds of models have been developed to represent
land C dynamics , the current generation of land C models inside ESMs all use multiple pools to represent various land C compartments and transfers among them. 

This common structure makes it possible to unify the land C models in a matrix form, by accommodating any number of pools, and by folding all C cycling processes into the terms of the matrix equation related to C input, C allocation into different plant organs, C turnover rate and its environmental modifier, and C transfers among pools.

## Installation

To get the current development version from github:

```R
# install.packages("devtools")
devtools::install_github("Zo-ran/MatrixMIP")
library("matrixMIP")
```

## Usage

```R
library(matrixMIP)

pool_n <- 15   #the number of pools

Initial_X <- rep(0, pool_n)

A <- diag(pool_n)
A[6, 1] <- -0.75
A[7, 1] <- -0.25
A[8, 2] <- -0.75
A[9, 2] <- -0.25
A[10, 3] <- -1
A[11, 4] <- -1
A[12, 5] <- -1
A[13, 6] <- -0.4125; A[13, 7] <- -0.45; A[13, 8] <- -0.3375; A[13, 9] <- -0.45; A[13, 10:12] <- -0.125; A[13, 14] <- -0.42
A[14, 6] <- -0.175; A[14, 8] <- -0.175; A[14, 10:12] <- -0.375; A[14, 13] <- -0.1664; A[14, 15] <- -0.45
A[15, 13] <- -0.004; A[15, 14] <- -0.03

B <- c(0.1725, 0.15, 0.0225, 0.14, 0.015, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)

K <- diag(c(0.09, 0.05, 0.01, 0.0008, 0.001, 0.153857, 1.2, 0.190296, 1.5, 0.0555556, 0.166667, 0.138889, 0.5865, 0.0162857, 0.000557143))

Tr <- matrix(0, pool_n, pool_n)

datelength4sasu <- 12
Scalar <- matrix(0, nrow = datelength4sasu, ncol = pool_n)
Scalar[, 1:5] <- 1
Scalar[, 6:15] <- c(0.001125247, 0.026934919, 0.00386343, 0.036797823, 0.056063114, 0.23379085, 0.168964848, 0.259775597, 0.019196622, 0.061492205, 0.013498124, 0.004996017)

Bscalar <- matrix(0, nrow = datelength4sasu, ncol = pool_n)
Bscalar[ , ] <- 1

Gpp <- c(1.5074, 12.8787, 33.9259, 78.7088, 142.6238, 171.0413, 194.6502, 203.5172, 145.2821, 96.9997, 27.6457, 7.9289)

initial <- sasu(100, "month", Initial_X, A, B, K, Tr, Scalar, Bscalar, Gpp)
initial
```


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


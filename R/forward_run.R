#' @title Forward Run

#' @description predict C pool sizes and C fluxes

#' @param Initial_X a vector of the initial size of each C pool, which can be obtained by "sasu"

#' @return a list

#' @examples sasu("month", x0, mat.A, mat.B0, mat.K, mat.Tr, scalar1, Bscalar1, I4simu0)

#' @export

forward_run <- function(Initial_X, A, B, K, Tr, Scalar, Bscalar, Gpp) {
  #### input check
  check_input(Initial_X, A, B, K, Tr, Scalar, Bscalar, Gpp)


  for(f in 1:datelength){
    scal.temp <- diag(x = 0, pn, pn); diag(scal.temp)=as.numeric(scalar1[f,])
    matrixAK = mat.A%*%(scal.temp%*%mat.K)+mat.Tr
    if(f==1) {g = f} else {g = f-1}
    # update of carbon pool sizes
    X[f,] = diag(pn)%*%X[g,]+(Bscalar1[f,]*mat.B)*I4simu[f]-matrixAK%*%X[g,]
    # Tch, chasing time
    diag(matrixAK)[diag(matrixAK)<1e-15] = 1e-15
    Tch = solve(matrixAK)
    # Tn, residence time of individual pool in network
    Tn[f,] = Tch%*%(Bscalar1[k,]*mat.B )
    # Xc, ecosystem storage capacity
    Xc[f,] = Tn[f,]*I4simu[f]
    # Xp, ecosystem storage potential
    Xp[f,] = Xc[f,] - X[f,]
    # NEC, net ecosystem exchange
    NEC[f,] = matrixAK%*%Xp[f,]
    # another way to calculate net ecosystem exchange
    # NEC2[f,] =  mat.B*I4simu[f]-matrixAK%*%X[g,]

    cflux$GPP[f] = I4simu0[f]
    # Ra, autotrophic respiration
    cflux$Ra[f] = (1-sum(Bscalar1[k,]*mat.B))*I4simu[f]


    # calculate C fluxes
    for(i in 1:length(nlst)) {

        # i = 1
        ordtemp = as.numeric(unlist(ordlst[i]))
        if(nlst[i]>=2) {
        cflux[f,i+6] = sum(colSums(mat.A2[,ordtemp])*diag(scal.temp)[ordtemp]*diag(mat.K)[ordtemp]*X[f,ordtemp]*diag(depthmatrix)[ordtemp])
        } else if(nlst[i]==1) {
        cflux[f,i+6] = sum(sum(mat.A2[,ordtemp])*diag(scal.temp)[ordtemp]*diag(mat.K)[ordtemp]*X[f,ordtemp]*diag(depthmatrix)[ordtemp])
        }
        if(nlst[i]>=1) { cflux[f,i+13] = sum(diag(mat.A2)[ordtemp]*diag(scal.temp)[ordtemp]*diag(mat.K)[ordtemp]*X[f,ordtemp]*diag(depthmatrix)[ordtemp]) }
    }

    for(j in 1:length(nlstP)) {
        ordtemp2 = unlist(ordlstP[j])
        cflux[f,j+20]=sum(Bscalar1[k,ordtemp2]*mat.B[ordtemp2])*I4simu[f]
    }
  }
  s
}

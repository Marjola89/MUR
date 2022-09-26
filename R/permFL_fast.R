#' Compute TFCE derived p-values using Freedman-Lane procedure faster!
#'
#' Given a NxV imaging matrix Y (N = number of subjects, V = number of vertices in the ventricular mesh), a NxC model matrix X  (N = number of subjects, C = number of variables + intercept term)
#' and the numbers of the column variables to extract, this function computes for each variable specified in extract a couple of arrays: in the first array the TFCE derived p-values map on the mesh and in the second array a binary variable expressing whether that vertex has reached significance using a whole-image threshold as specified in the TFCE original paper. The output is a matrix with a number of
#' columns equal to twice the lenght of extract and a number of rows equal to the number of vertices.
#' @param X is the design matrix. Number of rows = number of subjects in the study, number of columns = number of vertices in the atlas. Numerical varable must be normalized to 0-mean and unit-standard deviation. Categorical variables must be coded using dummy coding. The first column should contain the intercept (all 1s).
#' @param Y is the imaging matrix. Number of rows = N. Number of columns = V.
#' @param extract is an array expressing which covariates in X you want to extract.
#' @param A A V-dimensional vector containing the area associated with a vertex, usually its Voronoi area.
#' @param NNmatrix  Nx2 matrix containing the mesh edges. Important: to speed up the execution please avoid repetitions like (A,B) and (B,A).
#' @param nPermutations number of permutations in the permutation test, default is 1000.
#' @param E is the TFCE parameter, by default fixed to 0.5.
#' @param H is the TFCE parameter, by default fixed to 2.
#' @import igraph
#' @import mutools3D
#' @import parallel
#' @import doParallel
#' @import foreach
#' @import Rcpp
#' @import RcppEigen
#' @import RcppArmadillo
#' @import plyr
#' @import float
#' @export
#' @examples TFCEresults = permFL_fast(X, Y, extract, A, NNmatrix, nPermutations = 1000, E=0.5, H=2)

permFL_fast <- function(X, Y, extract, A, NNmatrix, nPermutations, E = 0.5, H = 2){

  registerDoParallel(detectCores())
<<<<<<< HEAD

=======
  
>>>>>>> 582365443a3ba1c27c2b7f749de0f7c76237114f
  set.seed(1234)
  # set seed for reproducibility

  Z <- X[, -extract]
  # compute Z (nuisance matrix)

  Zpinv <- Z %*% solve(t(Z) %*% Z) %*% t(Z)
  # Ypr<-eigenMapMatMult(Zpinv,Y) # faster multiplication!
  Ypr <- Zpinv %*% Y
  Yp <- Y - Ypr

  Rz <- fl(Yp) # float for precision matrix

  # parallelization for iF=1:10 in a loop
  nPerm <- nPermutations / 10

  permresP<-vector(mode = "list", length = nPerm) # a list of nPerm length
  for (iR in 1:nPerm){
    resP <- foreach(iF = 1:10, .combine=rbind)%dopar%{
      Yper <- Rz@Data[sample(1:nrow(Rz@Data)),]
      # Y permuted for the Freedman and Lane procedure

      resMUR <- murq(X, Yper, extract)
      rm(Yper)
      computed <- matrix(0, ncol = ncol(Y), length(extract))

      computed <- TFCE(h = round(resMUR[, 2], 2), A = A, NNmatrix = NNmatrix, E = E, H = H)
      #compute TFCE
      return(computed)
    }
    permresP[[iR]] <- resP

  }

  resP <- ldply(permresP, rbind)
  # for each element of a list, apply function then combine results into the resP data frame
  closeAllConnections()

  significance <- matrix(0, ncol = length(extract)*2, nrow = ncol(Y))
  # TFCE derived p-values

  results <- matrix(0, ncol = 3*length(extract), nrow = ncol(Y))
  results <- murq(X, Y, extract)
  tfceScores <- list()

  # compute the residual matrix of Z
  iEx <- 1
  tfceScores[[iEx]] <- TFCE(results[, 2+(iEx-1)*3], A = A, NNmatrix = NNmatrix, E = E, H = H)

  # list of TFCE scores to analyse
  TFCEmatrix <- resP[seq(1, nrow(resP), by=length(extract)),]


  minimum = sort(apply(TFCEmatrix, 1, min))
  if (length(which(minimum < 0) > 0)) {
    thrMin = minimum[ceiling(0.05 * nrow(TFCEmatrix))]
  } else {
    thrMin = 0
  }

  maximum = sort(apply(TFCEmatrix, 1, max))
  if (length(which(maximum > 0) > 0)) {
    thrMax = maximum[floor(0.95 * nrow(TFCEmatrix))]
  }else{
    thrMax = 0
  }

  for(a in 1:ncol(Y)){
    if(tfceScores[[iEx]][a] >= 0){
      significance[a, 1+(iEx-1)*2]  <- length(which(TFCEmatrix[, a]> tfceScores[[iEx]][a]))/nPermutations
      if(tfceScores[[iEx]][a] > thrMax) significance[a, 2+(iEx-1)*2] = 1
    }

    if(tfceScores[[iEx]][a] < 0){
      significance[a, 1+(iEx-1)*2]  <- length(which(TFCEmatrix[, a]< tfceScores[[iEx]][a]))/nPermutations
      if(tfceScores[[iEx]][a] < thrMin) significance[a, 2+(iEx-1)*2] = 1
    }

    if(significance[a, 1+(iEx-1)*2] == 0) significance[a, 1+(iEx-1)*2] <- 1 / nPermutations # minimum pvalue achievable.
  }



  length(which(significance[,1] < 0.05))
  TFCEresults = list("pvalues" = significance, "TFCEmatrix" = TFCEmatrix, "tfceScores" = tfceScores)

  #else TFCEresults = significance

  rm(significance)
  rm(TFCEmatrix)
  rm(tfceScores)

  return(TFCEresults)
}

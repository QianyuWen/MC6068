##' Low rank matrix completion
##'
##' This package provide SOFT-IMPUTE algorithm for the nuclear norm regularized least-squares problem that scales to large problems. Iteratively, this algorithm replaces the missing elements with those obtained from a soft-thresholded SVD. At every iteration SOFT-IMPUTE decreases the value of the objective function towards its minimum, and at the same time gets closer to the set of optimal solutions of the problem. Reference: Rahul Mazumder, Trevor Hastie, and Robert Tibshirani, 2010, Journal of Machine Learning Research 11(2010)2287-2322.
##' 
##' @title MC6068
##' @param x Incomplete matrix of your interest
##' @param iniMatrix Matrix with initial values
##' @param error Error for convergence
##' @param lambda Regularization parameter controlling the nuclear norm of the minimizer
##' @param silence Choice to mute iteration. T: mute; F: track each iteration
##' @return Complete matrix
##' @author Ye Yuan, Qianyu Wen
##' @export
##' @examples AdmmMC(x=incompleteMatrix, iniMatrix = NULL, error = 1e-16, lambda = 0.5, silence = F)
##' 

#' @export
AdmmMC <- function(x, iniMatrix = NULL, error = 1e-16, lambda = 0.5, silence = F){
  
  nu2norm <- function(p){(norm(p, type = "F"))^2}
  
  diffchange <- function(p, q){
    nu2norm(p-q)/nu2norm(q)
  }
  
  softThreshold <- function(yy, lambda){
    sign(yy) * pmax(abs(yy) - lambda, 0)
  }
  
  f <- function(z,w,lambda){nu2norm(z-w) + lambda*sum(svd(z)$d)}
  
  if(is.null(iniMatrix)){
    init <- (min(x, na.rm = T) + max(x, na.rm = T))/2
    z <- matrix(init, nrow = nrow(x), ncol = ncol(x))
  } else {
    z <- iniMatrix
  }
  z_old <- matrix(0, nrow = nrow(x), ncol = ncol(x))
  zsvd <- svd(z)
  u_z <- zsvd$u
  d_z <- zsvd$d
  v_z <- zsvd$v
  step <- 0
  
  while (diffchange(z,z_old)>error){
    w <- x
    z_old <- z
    for (i in 1:nrow(w)) {
      for (j in 1:ncol(w)) {
        w[i,j] <- ifelse(is.na(w[i,j]), z[i,j], w[i,j])
      }
    }
    
    wsvd <- svd(w)
    
    for (i in 1:length(d_z)) {
      yy <- t(u_z[,i]) %*% w %*% v_z[,i]
      d_z[i] <- softThreshold(yy = yy, lambda = lambda)
      u_z[,i] <- wsvd$u[,i]
      v_z[,i] <- wsvd$v[,i]
    }
    
    z <- u_z %*% diag(d_z) %*% t(v_z)
    step <- step + 1
    if(silence==F){
      cat("Step:", step,"\n", "OFV:", f(z=z, w=w ,lambda=lambda),"\n")
    }
  }
  z[z<0] <- 0
  z
}
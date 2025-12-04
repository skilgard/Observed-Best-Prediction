library(Matrix)
library(matlib)

## X is design matrix for fixed effects
## consists of predictor variables 
X <- cbind(c(rep(1,2), rep(0,7)),
           c(rep(0,2), rep(1,3), rep(0, 4)),
           c(rep(0,5), rep(1,4)))

## Z is a design matrix for random effect
## determined by study 
Z <- cbind(c(1, rep(0, 8)),
           c(rep(0,2), 1, rep(0,6)),
           c(rep(0,5), rep(1,2), rep(0,2)),
           c(0,1,0,1,1,0,0,1,1))

## observed response 
y <- c(110, 100, 110, 100, 100, 110, 110, 100, 100)

theta <- c(1, 1)

Ke <- diag(9)
Kb <- diag(4)

# Kb is the fixed matrix used to calculate D. It is calculated from the data
# Ke is the fixed matrix used to calculate R. It is calculated from the data
# b-hat is random effect vector
Calc_EBLUP <- function(Kb, Ke, y, X, Z, theta, method = c("QR", "proj"))
{
  covmat <- Calc_var_comp(Kb, Ke, y, X, Z, theta, method = method)
  D <- covmat[[1]]
  R <- covmat[[2]]
  beta <- Calc_beta(Kb, Ke, y, X, Z, theta, method = method)

  b_hat <- D %*% t(Z) %*% inv(Z %*% D %*% t(Z) + R) %*% (y - X %*% beta)

  return(b_hat)
}

## function to calculate V-hat
## estimated variance matrix for vector y
Calc_V <- function(Kb, Ke, y, X, Z, theta, method = c("QR", "proj")) {
  covmat <- Calc_var_comp(Kb, Ke, y, X, Z, theta, method = method)
  D <- covmat[[1]]
  R <- covmat[[2]]
  V <- Z %*% D %*% t(Z) + R
  return(V)
}

## beta_hat is the estimated fixed effect parameter 
Calc_beta <- function(Kb, Ke, y, X, Z, theta, method = c("QR", "proj"))
{
  V <- Calc_V(Kb, Ke, y, X, Z, theta, method = method)
  beta <- inv(t(X) %*% inv(V) %*% X) %*% t(X) %*% inv(V) %*% y
  return(beta) 
}

## this function returns the REML estimators for sigma_b^2 and sigma_e^2
Calc_var_comp <- function(Kb, Ke, y, X, Z, theta, num_iter = 1e3, tol = 1e-5, method = c("QR", "proj"))
{
  
    
    ei <- eigen(Ke) # Will Ke always be positive. If so, a symmetric = TRUE parameter can be added
    V <- ei$vectors
    Ke.inv.sqrt <- V %*% diag(1 / sqrt(ei$values)) %*% t(V)

    y.new <- Ke.inv.sqrt %*% y
    X.new <- Ke.inv.sqrt %*% X
    Z.new <- Ke.inv.sqrt %*% Z


    K <- Z.new %*% Kb %*% t(Z.new)
    cov_list <- LMM_REML(y = y.new, X=X.new, theta = theta, K = K, num_iter = num_iter, tol = tol, method = "QR")
    vc_reml <- cov_list$Theta[cov_list$num_iter,]  # extract the REML estimators for the variance components from the cov_list

    D <- vc_reml[1] * Kb
    R <- vc_reml[2] * Ke
    return_list <- list( D = D, R = R)

    return(return_list)
}


# trace function
tr <- function(A)
{
  return(sum(diag(A)))
}

##Here the function has a K parameter, but I now have a Kb and a Ke parameter. 
## I don't know what I should modify in this function 
## using the matricies you gave me, I'm pretty sure they need to be different sizes
LMM_REML <- function(y, X, theta, K, num_iter = 1e3, tol = 1e-5, method = c("QR", "proj"))
{
  n <- length(y)
  q <- rankMatrix(X)
  
  id_mat <- diag(rep(1,n-q))
  
  
  # matrix for parameter updates in EM algorithm
  Theta <- matrix(0, num_iter, 2)
  Theta[1,] <- theta
  
  for(i in 2:num_iter)
  {
    sigma_a.i <- Theta[i-1,1]
    sigma_e.i <- Theta[i-1,2]
    gamma.i   <- sigma_a.i/sigma_e.i
    
    if(method == "QR")
    {
      # error contrast matrix
      qrX <- qr(X)
      y.tilde <- qr.qty(qrX,y)[(rankMatrix(X)+1):n]
      A.t <-t(qr.Q(qrX,complete=TRUE)[,(rankMatrix(X)+1):n])
      A <- t(A.t)
      
      # marginal covariance matrix for y
      V_gamma <- A.t %*% K %*% A
      Sigma_gamma <- gamma.i * V_gamma + id_mat
      Sigma_gamma.inv <- solve(Sigma_gamma)
      
      # update parameters
      sigma_a.i1 <- sigma_a.i + 
        1/(n-q) * gamma.i^2 * t(y.tilde) %*% Sigma_gamma.inv %*% V_gamma %*% Sigma_gamma.inv %*% y.tilde -
        1/(n-q) * sigma_a.i^2/sigma_e.i * tr(Sigma_gamma.inv %*% V_gamma)
      
      
      sigma_e.i1 <- 1/(n-q) * (sigma_a.i * tr(V_gamma) - sigma_a.i^2/sigma_e.i * tr(V_gamma %*% Sigma_gamma.inv %*% V_gamma)
                               + crossprod(t(id_mat - gamma.i * Sigma_gamma.inv %*% V_gamma) %*% y.tilde))
      
    }
    
    if(method == "proj")
    {
      In <- diag(1,n)
      B_gamma <- gamma.i * K + In
      B_gamma.inv <- solve(B_gamma)
      P_gamma <- B_gamma.inv - B_gamma.inv %*% X %*% 
                 solve(t(X) %*% B_gamma.inv %*% X) %*% t(X) %*% B_gamma.inv
      P_X <- X %*% solve(crossprod(X)) %*% t(X)
      
      sigma_a.i1 <- sigma_a.i +
                    1/(n-q) * gamma.i^2 * t(y) %*% P_gamma %*% K %*% P_gamma %*% y -
                    1/(n-q) * sigma_a.i^2/sigma_e.i * tr(P_gamma %*% K)
      
      sigma_e.i1 <- 1/(n-q) * (sigma_a.i * tr(K %*% (In - P_X)) - sigma_a.i^2/sigma_e.i * tr(P_gamma %*% K %*% (In-P_X) %*% K)
                               + t(y) %*% (In - gamma.i * P_gamma %*% K) %*% (In-P_X) %*% 
                                 (In - gamma.i * K %*% P_gamma) %*% y)
      
    }
    
    Theta[i,1] <- sigma_a.i1
    Theta[i,2] <- sigma_e.i1
    
    cat(Theta[i,], "\n")
    
    if(crossprod(Theta[i,]-Theta[i-1,]) < tol)
    {
      break
    }
    
  }
  
  return_list <- list(Theta = Theta, num_iter = i)
  return(return_list)
}


Calc_var_comp(Kb, Ke, y, X, Z, theta)

Calc_V(Kb, Ke, y, X, Z, theta)

Calc_beta(Kb, Ke, y, X, Z, theta)

Calc_EBLUP(Kb, Ke, y, X, Z, theta)

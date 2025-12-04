library(MASS)
source("MSPE_and_Inverse_Fn.R")

# This function takes a matrix and finds the inverse 
# if it's positive semi definite, it uses the Cholensky method 
# if not, it uses ginv from MASS
TC_Inverse <- function(A) {
  tryCatch( 
    {
      return(Inverse(A))
    },
    error=function(e) {
      message('Matrix not positive semi-definite')
      return(ginv(A))
    },
    error=function(e) {
      message('ginv() did not work')
      return(solve(A))
    }
  )
}

# Capitol letter (eg M) is a matrix
# Capitalized word/phrase (eg Phi) is an expression with matricies 
# lowercase word/phase (eg sigmaE or gamma) is a scalar

# random error covariance = sigmaE * Q
# random effects covariance = sigmaV * K

## test the equation for sigma from the paper against 0.4 (the true variance component)


## because the derivatives of MSPE are so long and repetitive, I've split them up for
## clarity and efficiency
## sigmaE is the covariance comp for the random error
## sigmaV is the covariance comp for the random effect 
## gamma is the signal to noise ratio (cov comp for standard error/ cov comp for random error)
## w is the mu from the book
Calc_SignaltoNoise <- function(w, tol = 1e-7, max_rep = 2500, m, t, mu, gamma = 1, K, Q, X, methodSigmaE){
  
  GammaVector <- c() # each new gamma appends to the end, running record of all the gammas 
  Phi_2der_minEigen <- c()
  DerMSPE_1_values <- c()
  DerMSPE_2_values <- c()
  
  
  n <- nrow(X) # used to create identity matrix to update sigmaE 
  
  for(i in 1:max_rep) {
    
    # updates Phi and all of its derivative
    RecurrentInverse_1 <- TC_Inverse(gamma * K + Q) # this reduced the number of inverses calculated
    Phi <- RecurrentInverse_1 %*% Q %*% Q %*% RecurrentInverse_1
    Phi_1der <- -1 * (RecurrentInverse_1 %*% K %*% Phi + Phi %*% K %*% RecurrentInverse_1)
    
    
    # updates beta function
    RecurrentInverse_2 <- TC_Inverse(t(X) %*% Phi %*% X)
    Beta <- RecurrentInverse_2 %*% t(X) %*% Phi %*% mu
    
    # MODIFY: choose one and comment out the other 
    # Equation 2.5 from the paper, as of 7/6
    # updates sigmaE
    V <- diag(n) + gamma * K
    V_inv <- TC_Inverse(V)
    P <- V_inv - V_inv %*% X %*% TC_Inverse(t(X) %*% V_inv %*% X) %*% t(X) %*% V_inv
    
    if(methodSigmaE == "EQ") {
      sigmaE <- (1 /(n - qr(X)$rank)) * t(mu) %*% P %*% mu
      sigmaE <- as.numeric(sigmaE)
    }
    
    if(methodSigmaE == "GIVEN")
      sigmaE <- 0.4 # actual sigmaE
    
    RecurrentExpression <- mu - X %*% Beta
    
    DerMSPE_1 <- (2 / n) * sigmaE %*% tr(K %*% Phi) - (1 / n) * t(RecurrentExpression) %*% RecurrentInverse_1 %*%
      K %*% RecurrentInverse_1 %*% R %*% R %*% RecurrentInverse_1 %*% RecurrentExpression
    
    DerMSPE_1_values <- append(DerMSPE_1_values, as.numeric(DerMSPE_1))
    
    DerMSPE_2 <- (2 / n) * sigmaE * tr(K %*% Phi_1der) + (1 / n) * t(RecurrentExpression) %*% RecurrentInverse_1 %*%
      K %*% RecurrentInverse_1 %*% K %*% Phi %*% RecurrentExpression + (1 / n) * t(RecurrentExpression) %*%
      Phi %*% K %*% RecurrentInverse_1 %*% K %*% RecurrentInverse_1 %*% RecurrentExpression -
      (1 / n) * t(RecurrentExpression) %*% Phi_1der %*% K %*% RecurrentInverse_1 %*% RecurrentExpression -
      (1 / n) * t(RecurrentExpression) %*% RecurrentInverse_1 %*% K %*% Phi_1der %*% RecurrentExpression 
    
    DerMSPE_2_values <- append(DerMSPE_1_values, as.numeric(DerMSPE_2))
    print(DerMSPE_2_values)
    
    
    
    gamma <- gamma - (t * (1/n) * DerMSPE_1 - gamma^(-1)) / (t * (1 / n) * DerMSPE_2 + gamma^(-2)) # update gamma
    gamma <- as.numeric(gamma)
    GammaVector <- append(GammaVector, gamma)
    
    # MSPEtest <- sigmaE * tr(Q) - sigmaE * tr(Q %*% RecurrentInverse_1 %*% Q) + t(mu - X %*% Beta) %*%
    #   RecurrentInverse_1 %*% t(Q) %*% Q %*% RecurrentInverse_1 %*% (mu - X %*% Beta)
    # print(MSPEtest)
    
    # if the if m/gamma is too large, converges too fast, not as accurate
    if(i >= max_rep || (m / t) < tol) {
      sigmaV <- gamma * sigmaE # gamma = varcomp of rand effect / varcomp of rand error
      returnList <- list(num_iter = i, signal_to_noise = gamma, sigmaE = sigmaE,
                         sigmaV = sigmaV, beta = Beta, prevGammas = GammaVector,
                         DerMSPE_1_values = DerMSPE_1_values,
                         DerMSPE_2_values = DerMSPE_2_values
      )
      return(returnList)
    }
    
    t <- as.numeric(t) * as.numeric(w) # update gamma
  }
}

# UUU <- Calc_SignaltoNoise(w = 2, tol = 1e-7, max_rep = 3, m = 1, 
#                     t = 12, mu = y, gamma = 1, K = K, Q = Q, X = X)

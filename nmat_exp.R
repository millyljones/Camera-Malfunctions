library(nimble)
library(expm)

factrial <- nimbleFunction(run = function(x = integer(0)){
  # nimble has protected 'factorial' somewhere
  
  returnType(double(0))
  return(exp(lfactorial(x)))
}
)

pow.mat <- nimbleFunction(run = function(C = double(2), 
                                          n = integer(0),
                                          n.states = integer(0)){
   # take matrix to the nth power
   returnType(double(2))
   
   Cn <- matrix(type = 'double', ncol = n.states, nrow = n.states)
   Cn[,] <- C[,]
   if (n == 0){
     return(diag(rep(1, n.states)))
   } else if (n == 1){
     return(C[,])
   } else {
     for (i in 2:n){
       Cn[,] <- Cn[,] %*% C[,]
     }
     return(Cn[,])
   }
 }
)

pow2.mat <- nimbleFunction(run = function(C = double(2), 
                                         n = integer(0),
                                         n.states = integer(0)){
  # take matrix to 2^n power
  returnType(double(2))
  
  Cn <- array(type = 'double', dim=c(n,n.states,n.states))
  Cn[1,,] <- C[,] %*% C[,]
  if (n == 0){
    return(diag(rep(1, n.states)))
  } else if (n == 1){
    return(Cn[1,,])
  } else {
    for (i in 2:n){
      Cn[i,,] <- Cn[i-1,,] %*% Cn[i-1,,]
    }
    return(Cn[n,,])
  }
}
)

r <- nimbleFunction(run = function(m = integer(0), 
                                   C = double(2),
                                   n.states = integer(0)){
  
  returnType(double(2))
  
  # initialise
  k <- m
  p <- matrix(type = 'double', ncol = n.states, nrow = n.states)
  q <- matrix(type = 'double', ncol = n.states, nrow = n.states)
  
  for (i in 0:m){
    mult <- factrial(k+m-i)*factrial(k)/(factrial(k+m)*factrial(k-i)*factrial(i))
    p[,] <- p[,] + mult*(pow.mat(C,i,n.states))
    q[,] <- q[,] + mult*(pow.mat(-C,i,n.states))
  }
  
  return(p[,] %*% inverse(q[,]))
  
}
)

norm1 <- nimbleFunction(run = function(C = double(2),
                                       n.states = integer(0)){
  
  returnType(double(0))
  
  norm <- sum(abs(C[,1]))
  
  if (n.states > 1){
    for (i in 2:n.states){
      norm.prop <- sum(abs(C[,i]))
      if (norm.prop > norm){
        norm <- norm.prop
      }
    }
  }
  
  return(norm)
  
})

balance_scale <- nimbleFunction(run = function(A = double(2), 
                                               n.states = integer(0)) {
  
  # balances a matrix through scaling to make the 1 norm of cols and rows similar
  # follows as best as I could understand the dgebal from LAPACK
  # returns a 3x6 matrix with first 3 columns the scaled matrix,
  # and last three columns the diagonal matrix of 
  # scaling used in each dimension
  
  returnType(double(2))
  
  sclfac <- 2.0
  factor <- 0.95
  converged <- FALSE
  
  scale <- numeric(length=n.states)
  scale[1:n.states] <- rep(1, n.states)
  
  if (sum(A[n.states,1:n.states])==0){
    ilo <- 1
    ihi <- 2
  } else {
    ilo <- 1
    ihi <- n.states
  }
  
  while (!converged) {
    
    converged <- TRUE
    
    for (i in ilo:ihi) {
      
      cval <- sum(abs(A[ilo:ihi, i])) 
      rval <- sum(abs(A[i, ilo:ihi])) 
      
      if (cval != 0 & rval != 0){
        g <- rval / sclfac
        f <- 1.0
        s <- cval + rval
        
        while (cval < g) {
          cval <- cval * sclfac
          rval <- rval / sclfac
          f <- f * sclfac
        }
        
        g <- cval / sclfac
        
        while (g >= rval) {
          cval <- cval / sclfac
          rval <- rval * sclfac
          f <- f / sclfac
          g <- cval / sclfac
        }
        
        if ((cval + rval) < factor * s) {
          scale[i] <- scale[i] * f
          A[i, ilo:ihi] <- A[i, ilo:ihi] / f
          A[ilo:ihi, i] <- A[ilo:ihi, i] * f
          converged <- FALSE
        }
      } 
    }
  }
  
  S <- matrix(nrow= n.states, ncol = n.states)
  S[1:n.states,1:n.states] <- diag(scale)  
  
  R <- matrix(nrow=n.states, ncol = 2*n.states)
  R[1:n.states, 1:n.states] <- A[1:n.states,1:n.states]
  R[1:n.states, n.states + (1:n.states)] <- S[1:n.states,1:n.states]
  
  return(R)
  
}
)

balance_perm <- nimbleFunction(run = function(A = double(2), 
                                              n.states = integer(0)){
  # function to create as close as possible a block diagonal matrix
  # by permuting rows and columns. Since A has a very pre-determined
  # form (transition matrix Q for the CTMC), this can be simplified
  # to the reduced function here.
  # Returns a 3x6 matrix with first three columns the permuted matrix
  # and last three columns the permutation matrix P that encodes the permutations
  
  returnType(double(2))
  
  P <- matrix(nrow=n.states, ncol = n.states)
  P[1:n.states, 1:n.states] <- diag(rep(1, n.states))
  
  if (A[2,1]==0 & A[1,3]!=0){
    A[c(2,3), 1:n.states] <- A[c(3,2), 1:n.states]
    A[1:n.states, c(2,3)] <- A[1:n.states, c(3,2)]
    P[1:n.states, 1:n.states] <- diag(rep(1, n.states))
    P[2,2] <- 0; P[2,3] <- 1
    P[3,2] <- 1; P[3,3] <- 0
  }
  
  R <- matrix(nrow=n.states, ncol = 2*n.states)
  R[1:n.states, 1:n.states] <- A[1:n.states,1:n.states]
  R[1:n.states, n.states + (1:n.states)] <- P[1:n.states,1:n.states]
  
  return(R)
  
}
)

nmat.exp.pade <- nimbleFunction(run = function(A = double(2),
                                               n.states = integer(0),
                                               balancing = integer(0, default = 0)){
  # code follows almost exactly expm.Higham08
  # in places modified as our matrix A is fixed square 3x3 and has predetermined
  # non-zero entries
  
  returnType(double(2))
  
  b <- c(64764752532480000,
         32382376266240000,
         7771770303897600,
         1187353796428800,
         129060195264000,
         10559470521600,
         670442572800,
         33522128640,
         1323241920,
         40840800,
         960960,
         16380,
         182,
         1)
  
  theta <- c(1.495585217958292e-2,
             2.539398330063230e-1,
             9.504178996162932e-1,
             2.097847961257068e0,
             5.371920851148152e0)
  
  m.opt <- c(3,5,7,9,13)
  
  if (balancing){
    
    AP <- matrix (type = 'double', ncol = n.states, nrow = n.states)
    AS <- matrix (type = 'double', ncol = n.states, nrow = n.states)
    P <- matrix (type = 'double', ncol = n.states, nrow = n.states)
    S <- matrix (type = 'double', ncol = n.states, nrow = n.states)
    R1 <- matrix (type = 'double', ncol = 2*n.states, nrow = n.states)
    R2 <- matrix (type = 'double', ncol = 2*n.states, nrow = n.states)
    
    R1[,] <- balance_perm(A[,], n.states)
    AP[,] <- R1[1:n.states, 1:n.states]
    P[,] <- R1[1:n.states, n.states + (1:n.states)]
    
    R2[,] <- balance_scale(AP[,], n.states)
    AS[,] <- R2[1:n.states, 1:n.states]
    S[,] <- R2[1:n.states, n.states + (1:n.states)]
    
    C <- matrix(type = 'double', ncol = n.states, nrow = n.states)
    C[,] <- AS[,] #- diag(rep(nu, n.states))
    
  } else {
    
    C <- matrix(type = 'double', ncol = n.states, nrow = n.states)
    C[,] <- A[,] #- diag(rep(nu, n.states))
    
  }
  
  C.norm <- norm1(C[,], n.states)
  
  C1 <- matrix(type = 'double', ncol = n.states, nrow = n.states)
  C2 <- matrix(type = 'double', ncol = n.states, nrow = n.states)
  C4 <- matrix(type = 'double', ncol = n.states, nrow = n.states)
  C6 <- matrix(type = 'double', ncol = n.states, nrow = n.states)
  U <- matrix(type = 'double', ncol = n.states, nrow = n.states)
  V <- matrix(type = 'double', ncol = n.states, nrow = n.states)
  r13 <- matrix(type = 'double', ncol = n.states, nrow = n.states)
  X <- matrix(type = 'double', ncol = n.states, nrow = n.states)
  
  if (C.norm <= theta[length(theta)]){ # if can use the short cut
    
    for (m in 1:5){
      if (C.norm <= theta[m]){
        X[,] = r(m.opt[m], C[,], n.states)
        X[,] = X[,] #exp(nu)*X[,]
        return(X[,])
      }
    }
    
  } else { # if need to use the scaling and squaring
    
    log.2 <- log(C.norm/theta[length(theta)])/log(2) # change of base formula
    # to get log base 2
    s <- ceiling(log.2)
    if (s > 0){
      C1[,] <- C[,]/(2^s)
    } else {
      C1[,] <- C[,]
    }
    C2[,] <- pow.mat(C1[,], 2, n.states)
    C4[,] <- pow.mat(C2[,], 2, n.states)
    C6[,] <- C2[,] %*% C4[,]
    U[,] <- C1[,] %*% ((C6[,]%*%(b[14]*C6[,] + b[12]*C4[,] + b[10]*C2[,])) +
                         b[8]*C6[,] + b[6]*C4[,] + b[4]*C2[,] + b[2]*diag(rep(1, n.states)))
    V[,] <- (C6[,] %*% (b[13]*C6[,] + b[11]*C4[,] + b[9]*C2[,])) + 
      b[7]*C6[,] + b[5]*C4[,] + b[3]*C2[,] + b[1]*diag(rep(1, n.states))
    
    r13[,] <- solve(V[,]-U[,], U[,]+V[,])
    # needs repeated squaring to solve
    # need to square it |s| times
    X[,] <- r13[,]
    if (s > 0){
      X[,] <- pow2.mat(X[,], s, n.states)
    }
    
    ## undo scaling
    if (balancing){
      
      d <- numeric(length = n.states)
      d <- diag(S)
      X[,] <- X[,] * d * rep(1/d, each = n.states)
    
      X[,] <- t(P[,]) %*% X[,] %*% P[,]
    }
    
    return(X[,])
  }
}
)

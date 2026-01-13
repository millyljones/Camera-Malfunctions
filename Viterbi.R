
# this file contains the functions to compute the probabilities
# that the system is in a given state at any point in time

# https://rdrr.io/cran/HiddenMarkov/man/mmpp-2nd-level-functions.html?utm_source=chatgpt.com
# check that these guys are doing approx the same thing we are doing with the MMPP

prob.det <- nimbleFunction(run = function(C = double(2),
                                          x = double(1),
                                          n.states = integer(0),
                                          Ltp = double(1), 
                                          Lfp = double(1), 
                                          type = integer(1),
                                          start = integer(1),
                                          collection = double(1)){
  # C matrix = Q - L
  # x is the times of all detections (including 0 at front and E at end)
  # for camera history (note not all the histories at the camera itself)
  # Ltp[i] = ltp is detection rate in state i for true positives
  # Lfp[i] = lfpi is detection rate in state i for false positives
  # type indicates true positive (1) or false neg (0) (first entry is NA, last is usually NA)
  # type has same length as x
  # start is a vector with a 1 indicating what state the camera started in
  # collection is a vector with a 1 indicating what state the camera ended in -
  # this is a vector of all 1's if we are not tracking state.at.collection
  # I know collection should be an integer but the function refused to compile so
  # a double it shall be

  # Uses the Viterbi algorithm to determine at every point in X (usually a detection point)
  # what probability the state is in. This includes the probabilities at 
  # t=0 and t=E. t=0 is always in state 1.
  # Returns a vector where entry j is the most probable state that point j is in
  
  returnType(integer(1))
  
  K <- length(x)
  
  E <- matrix(type = 'double', nrow = K, ncol = n.states) # initialise matrix
  E[1,] <- start[1:n.states]
  
  for (i in 2:K){

    x.diff <- x[i]-x[i-1]
    e_Cx <- matrix(type='double', nrow = n.states, ncol = n.states)
    e_Cx[,] <- nmat.exp.pade(C*x.diff, n.states)

    if (i != K){ # if not dealing with trailing observation
     
       E.mat <- matrix(type='double', nrow = n.states, ncol = n.states)
       E.mat[,] <- diag(E[i-1,]) %*% e_Cx[,] %*% ((type[i]*diag(Ltp[1:n.states])) + ((1-type[i])*diag(Lfp[1:n.states])))
     
       for (j in 1:n.states){ # nimble does not have an apply function
         E[i,j] <- max(E.mat[,j])  # take maximum over columns
       }
     
    } else { # if dealing with last observation
     
       E.mat <- matrix(type='double', nrow = n.states, ncol = n.states)
       E.mat[,] <- diag(E[i-1,]) %*% e_Cx[,] %*%  diag(collection[1:n.states])
     
       # if collection (1,1,1) then this doesn't have an impact and model doesn't
       # learn end state. If collection wlog (1,0,0) then only first column of E.mat
       # is non-zero.
     
       for (j in 1:n.states){ # nimble does not have an apply function
         E[i,j] <- max(E.mat[,j])  # take maximum over columns
       }

    }

     E[i,] <- E[i,]/sum(E[i,]) # numerical stability
  } # created matrix of E's
  
  states <- integer(length = K)
  states[K] <- which(E[K,] == max(E[K,]))[1] # complains during compile
  # # when I don't put [1]
  # 
  for (i.s in 1:(K-1)){ # C++ compilation can't index backwards, 
     # so sort that out with i.s to i transformation
     i <- K - i.s
     
     x.diff <- x[i+1]-x[i]
     e_Cx <- matrix(type='double', nrow = n.states, ncol = n.states)
     e_Cx[,] <- nmat.exp.pade(C*x.diff, n.states)
     
     states[i] <- which(E[i,]*e_Cx[,states[i+1]] == max(E[i,]*e_Cx[, states[i+1]]))[1]
  }
  
  return(states[])
  
})

prob.det.notype <- nimbleFunction(run = function(C = double(2),
                                          x = double(1),
                                          n.states = integer(0),
                                          L = double(1), 
                                          start = integer(1),
                                          collection = double(1)){
  # C matrix = Q - L
  # x is the times of all detections (including 0 at front and E at end)
  # for camera history (note not all the histories at the camera itself)
  # L[i] = is detection rate in state i
  # start is a vector with a 1 indicating what state the camera started in
  # collection is a vector with a 1 indicating what state the camera ended in -
  # this is a vector of all 1's if we are not tracking state.at.collection

  # Uses the Viterbi algorithm to determine at every point in X (usually a detection point)
  # what probability the state is in. This includes the probabilities at 
  # t=0 and t=E. t=0 is always in state 1.
  # Returns a vector where entry j is the most probable state that point j is in
  
  returnType(integer(1))
  
  K <- length(x)
  
  E <- matrix(type = 'double', nrow = K, ncol = n.states) # initialise matrix
  E[1,] <- start[1:n.states]
  
  for (i in 2:K){
    
    x.diff <- x[i]-x[i-1]
    e_Cx <- matrix(type='double', nrow = n.states, ncol = n.states)
    e_Cx[,] <- nmat.exp.pade(C*x.diff, n.states)
    
    if (i != K){ # if not dealing with trailing observation
      
      E.mat <- matrix(type='double', nrow = n.states, ncol = n.states)
      E.mat[,] <- diag(E[i-1,]) %*% e_Cx[,] %*% diag(L)
      
      for (j in 1:n.states){ # nimble does not have an apply function
        E[i,j] <- max(E.mat[,j])  # take maximum over columns
      }
      
    } else { # if dealing with last observation
      
      E.mat <- matrix(type='double', nrow = n.states, ncol = n.states)
      E.mat[,] <- diag(E[i-1,]) %*% e_Cx[,] %*% diag(collection[1:n.states])
      
      # if collection (1,1,1) then this doesn't have an impact and model doesn't
      # learn end state. If collection wlog (1,0,0) then only first column of E.mat
      # is non-zero.
      
      for (j in 1:n.states){ # nimble does not have an apply function
        E[i,j] <- max(E.mat[,j])  # take maximum over columns
      }
      
    }
    
    E[i,] <- E[i,]/sum(E[i,]) # numerical stability
  } # created matrix of E's
  
  states <- integer(length = K)
  states[K] <- which(E[K,] == max(E[K,]))[1] # complains during compile
  # when I don't put [1]
  
  for (i.s in 1:(K-1)){ # C++ compilation can't index backwards, 
    # so sort that out with i.s to i transformation
    i <- K - i.s
    
    x.diff <- x[i+1]-x[i]
    e_Cx <- matrix(type='double', nrow = n.states, ncol = n.states)
    e_Cx[,] <- nmat.exp.pade(C*x.diff, n.states)
    
    states[i] <- which(E[i,]*e_Cx[,states[i+1]] == max(E[i,]*e_Cx[, states[i+1]]))[1]
  }
  
  return(states)
  
})

prob_in_state = nimbleFunction(run = function(C = double(2), 
                                              y1 = double(0), 
                                              y2 = double(0), 
                                              state = integer(0),
                                              n.states = integer(0)){
  # C is the matrix Q-L for Q transition matrix and L the diagonal 
  # matrix of total detection rates in each state
  # y1 is distance from previous arrival/detection to t
  # y2 is distance from t to next arrival/detection
  # state is the state we want to find the probability of being in 
  
  # Returns a matrix where entry i,j is the probability of being in 
  # state `state` after distance y1 from a detection in state i, 
  # and where the next detection is in distance y2 and in state j
  
  returnType(double(2))
  
  y <- y1+y2
  
  e_Cy <- matrix(type = 'double', nrow = n.states, ncol = n.states)
  e_Cy[,] <- nmat.exp.pade(C*y, 3)
  
  e_Cy1 <- matrix(type = 'double', nrow = n.states, ncol = n.states)
  e_Cy1[,] <- nmat.exp.pade(C*y1, 3)
  
  e_Cy2 <- matrix(type = 'double', nrow = n.states, ncol = n.states)
  e_Cy2[,] <- nmat.exp.pade(C*y2, 3)
  
  Pt <- matrix(type = 'double', nrow = n.states, ncol = n.states)
  for (i in 1:n.states){
    for (j in 1:n.states){
      Pt[i,j] <- e_Cy1[i,state]*e_Cy2[state,j]/e_Cy[i,j]
    }
  }
  
  return(Pt[,])
  
})

prob_t <- nimbleFunction(run=function(dets = double(1),
                                      type = integer(1),
                                      res = integer(0),
                                      ltp = double(1),
                                      lfp1 = double(1),
                                      lfp3 = double(1),
                                      mu12 = double(1),
                                      mu13 = double(1),
                                      mu21 = double(1),
                                      mu32 = double(1),
                                      n.states = integer(0),
                                      state = integer(0),
                                      E = double(0),
                                      start = integer(1),
                                      collection = double(1)
){
  
  # dets is a vector of detection times (including 0 and E at start and end
  # respectively)
  # type is a vector of true positive and false positive labels (starts NA and
  # ends NA). type is same length as dets.
  # res indicates the number of points to compute probabilities at between 0 and E
  # essentially the resolution to plot against.
  # ltp - true pos rate; lfp1, lfp3 - false pos rates for states 1 and 3
  # mu12;mu13;mu21;mu32 - transition rates for underlying CTMT.
  # ltp through mu32 are vectors of same length (usually posterior samples)
  # n.states - number of states in the system
  # state - the state we are interested in computing the probability for
  # E - total length of the survey
  # start - vector with 1 indicating state of camera at start of survey
  # collection - vector with 1 indicating state of camera at end of survey
  # collection = c(1,1,1) if state at end of survey not recorded
  
  # Returns a matrix of probabilities
  # nrow is number of posterior samples
  # ncol is the number of points across time to compute
  
  returnType(double(2))
  
  prob_times <- seq(0.01, E-0.1, length.out = res) # for numerical stability
  
  N <- length(mu12) # number of posterior samples
  probs1 <- matrix(type = 'double', nrow = N, ncol = res)
  
  for (k in 1:N){ # loop over posterior samples
    
    Ltp <- c(ltp[k], 0, ltp[k])
    Lfp <- c(lfp1[k], 0, lfp3[k])
    L <- matrix(type = 'double', nrow = n.states, ncol = n.states)
    L[,] <- diag(Ltp + Lfp)
    mu12.i <- mu12[k]
    mu13.i <- mu13[k]
    mu21.i <- mu21[k]
    mu32.i <- mu32[k]
    Q <- matrix(type = 'double', nrow = n.states, ncol = n.states)
    Q[1,1] <- -(mu12.i + mu13.i)
    Q[1,2] <- mu12.i
    Q[1,3] <- mu13.i
    Q[2,1] <- mu21.i
    Q[2,2] <- -(mu21.i)
    Q[2,3] <- 0
    Q[3,1] <- 0
    Q[3,2] <- mu32.i
    Q[3,3] <- -(mu32.i)
    
    states <- numeric(length(dets))
    
    # Viterbi algorithm gives only most likely sequence of states
    states[] <- prob.det(Q-L, dets, n.states, 
                         Ltp, Lfp, 
                         type, start, collection)

    for (j in 1:res){ # loop through time points
      
      # find time since and to nearest detections
      time <- prob_times[j]
      y1 <- time - dets[dets <= time][sum(dets <= time)] # time since previous detection
      y2 <- dets[dets >= time][1] - time # time till next detection
      
      s1 <- states[sum(dets <= time)] # get previous state
      s2 <- states[1+sum(dets <= time)] # get next state
      
      probs1[k,j] <- prob_in_state(Q-L, y1, y2, state, n.states)[s1, s2]
    }
  }
  
  return(probs1[,])
  
})

prob_t.notype <- nimbleFunction(run=function(dets = double(1),
                                      res = integer(0),
                                      l1 = double(1),
                                      l3 = double(1),
                                      mu12 = double(1),
                                      mu13 = double(1),
                                      mu21 = double(1),
                                      mu32 = double(1),
                                      n.states = integer(0),
                                      state = integer(0),
                                      E = double(0), 
                                      start = integer(1),
                                      collection = double(1)
){
  
  # dets is a vector of detection times (including 0 and E at start and end
  # respectively)
  # res indicates the number of points to compute probabilities at between 0 and E
  # essentially the resolution to plot against.
  # l1, l3- detection rates for states 1 and 3
  # mu12;mu13;mu21;mu32 - transition rates for underlying CTMT.
  # l1 through mu32 are vectors of same length (usually posterior samples)
  # n.states - number of states in the system
  # state - the state we are interested in computing the probability for
  # E - total length of the survey
  # start - vector with 1 indicating state of camera at start of survey
  # collection - vector with 1 indicating state of camera at end of survey
  # collection = c(1,1,1) if state at end of survey not recorded
  
  # Returns a matrix of probabilities
  # nrow is number of posterior samples
  # ncol is the number of points across time to compute
  
  returnType(double(2))
  
  prob_times <- seq(0.1, E-0.1, length.out = res) # for numerical stability
  
  N <- length(mu12) # number of posterior samples
  probs1 <- matrix(type = 'double', nrow = N, ncol = res)
  
  for (k in 1:N){ # loop over posterior samples
    
    L <- matrix(type = 'double', nrow = n.states, ncol = n.states)
    L[,] <- diag(c(l1[k], 0, l3[k]))
    mu12.i <- mu12[k]
    mu13.i <- mu13[k]
    mu21.i <- mu21[k]
    mu32.i <- mu32[k]
    Q <- matrix(type = 'double', nrow = n.states, ncol = n.states)
    Q[1,1] <- -(mu12.i + mu13.i)
    Q[1,2] <- mu12.i
    Q[1,3] <- mu13.i
    Q[2,1] <- mu21.i
    Q[2,2] <- -(mu21.i)
    Q[2,3] <- 0
    Q[3,1] <- 0
    Q[3,2] <- mu32.i
    Q[3,3] <- -(mu32.i)
    
    states <- numeric(length(dets))
    # Viterbi algorithm gives only most likely sequence of states
    states[] <- prob.det.notype(Q-L, dets, n.states, 
                       c(l1[k], 0, l3[k]), 
                       start, collection)
    
    for (j in 1:res){ # loop through time points
      
      # find time since and to nearest detections
      time <- prob_times[j]
      y1 <- time - dets[dets <= time][sum(dets <= time)] # time since previous detection
      y2 <- dets[dets >= time][1] - time # time till next detection
      
      s1 <- states[sum(dets <= time)] # get previous state
      s2 <- states[1+sum(dets <= time)] # get next state
      
      probs1[k,j] <- prob_in_state(Q-L, y1, y2, state, n.states)[s1, s2]
    }
  }
  
  return(probs1[,])
  
})

cprob_t <- compileNimble(prob_t)
cprob_t.notype <- compileNimble(prob_t.notype)


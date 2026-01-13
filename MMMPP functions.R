library(nimble)
library(ggplot2)
library(stringr)
library(gridExtra)

# MMMPP functions
# contains the functions to run the MMMPP analysis of the camera
# trap histories
# This is the continuous analogue of the HMM formulation

build_data <- function(data, meta,
                       fix_transitions = list(mu21.0 = F, mu13.0 = F),
                       state.at.collection = T, 
                       intermediate.checks = NULL,
                       marks = T){
  
  ############################################################
  ############ Want to remove any 'Work' entries from data
  data <- data[!(data$Camera_Installation %in% 'Work'), ]
  
  ############################################################
  ############ Order camera numbers by those in meta
  meta <- meta[order(meta$Site_N),]
  
  ###########################################################
  ###### If state at collection need to have column camera_fail
  if (state.at.collection){
    if (!('camera_fail' %in% colnames(meta))){
      stop('If state.at.collection is true, expect meta to have column called `camera_fail`')
    }
    # camera_checks_end <- data.frame(Site_ID = meta$Site_ID,
    #                             Date.Time = meta$camera_off,
    #                             state = ifelse(is.na(meta$camera_fail), 1, 2),
    #                             repaired = ifelse(is.na(meta$camera_fail), NA, 0))
  }
  
  ############################################################
  ############ If intermediate checks are given:
  # Will create a new camera history for each period between
  # known states.
  meta_t <- meta
  meta_t$Site_ID_t <- meta_t$Site_ID
  data_t <- data
  data_t$Site_ID_t <- data$Site_ID
  # if camera_fail doesn't exist as a column, create one:
  if (!('camera_fail' %in% colnames(meta_t))){
    meta_t$camera_fail <- NA # create one of NA's for now
  }
  
  if (!is.null(intermediate.checks)){ 
    # get checks in camera and date order
    intermediate.checks <- intermediate.checks[order(intermediate.checks$Site_ID, intermediate.checks$Date.Time),]
    affected_cameras <- unique(intermediate.checks$Site_ID)

    for (cam in affected_cameras){
      inter.checks.cam  <- intermediate.checks[intermediate.checks$Site_ID==cam,]
      n.ac <- nrow(inter.checks.cam)+ 1 # number of entries for this camera
      ID_t <- paste0(cam, '.', as.character(1:(n.ac)))
      start.times <- c(meta$camera_on[meta$Site_ID==cam], inter.checks.cam$Date.Time)
      end.times <- c(inter.checks.cam$Date.Time, meta$camera_off[meta$Site_ID==cam])
      L_t <- difftime(end.times, start.times, units = 'days')
      
      ###### check this camera meets conditions
      if (fix_transitions$mu13.0 == T & any(inter.checks.cam$state==3) ){
        stop(paste0('Error: for camera ',
                    cam,
                    ', fix_transitions$mu13.0 is FALSE but intermediate.checks shows camera being in state 3'))
      }
      
      if (nrow(inter.checks.cam) > 1){
        if (state.at.collection == F){
          transitions <- paste0(as.character(c(1, inter.checks.cam$state[1:(n.ac-2)])), # camera starts working
                                as.character(c(inter.checks.cam$state[1:(n.ac-1)])) )
        } else {
          transitions <- paste0(as.character(c(1, inter.checks.cam$state[1:(n.ac-1)])), # camera starts working
                                as.character(c(inter.checks.cam$state[1:(n.ac-1)], 
                                               ifelse(!is.na(meta$camera_fail[meta$Site_ID==cam]),
                                                      2, 1)) )) # assuming only fail or fine
        }
      } else {
        if (state.at.collection == T){
          transitions <- paste0(as.character(c(1, inter.checks.cam$state[1])),
                                as.character(c(inter.checks.cam$state[1], 
                                               ifelse(!is.na(meta$camera_fail[meta$Site_ID==cam]), 2, 1))))
        } else {
          transitions <- paste0(as.character(c(1)), as.character(inter.checks.cam$state[1]))
        }
      }
      
      if (fix_transitions$mu21.0 == T & (any(transitions %in% c('21', '23')) | 
                                       any(inter.checks.cam$repaired == 1, na.rm = T))){
        stop(paste0('Error: for camera ',
                    cam,
                    ', fix_transitions$mu21.0 is FALSE but intermediate.checks shows cameras being repaired'))
      }
      
      ######## add new meta observations
      meta_t_add <- data.frame(matrix(nrow=n.ac, ncol = ncol(meta_t)))
      colnames(meta_t_add) <- colnames(meta_t)
      meta_t_add$Site_ID <- cam
      meta_t_add$camera_on <- start.times
      meta_t_add$camera_off <- end.times
      meta_t_add$camera_L <- as.numeric(L_t)
      meta_t_add$Site_N <- meta$Site_N[meta$Site_ID==cam]
      meta_t_add$Site_ID_t <- ID_t
      if (T){ # changed my mind, always want to do this if have checks
        camera_fail_t <- c(ifelse((inter.checks.cam$state==2), # vectorised
                                  end.times[1:(n.ac-1)], rep(NA, n.ac-1)), 
                           meta$camera_fail[meta$Site_ID==cam])
        meta_t_add$camera_fail <- as.POSIXct(camera_fail_t)
      }
      meta_t <- meta_t[meta_t$Site_ID!=cam,]
      meta_t <- rbind(meta_t, meta_t_add)
      
      # meta_t has been created, now need to ammend the data detections
      # to the new corresponding camera names
      data_t.cam <- data_t[data_t$Site_ID==cam,]
      for (i in 1:n.ac){ # loop through new camera histories
        data_t.cam$Site_ID_t[data_t.cam$Date.Time <= meta_t$camera_off[meta_t$Site_ID==cam][i] & 
                               data_t.cam$Date.Time >= meta_t$camera_on[meta_t$Site_ID==cam][i]] <-
          meta_t$Site_ID_t[meta_t$Site_ID==cam][i]
      }
      
      data_t <- data_t[data_t$Site_ID!=cam,]
      data_t <- rbind(data_t, data_t.cam)
    }
  }
  
  #############################################################
  ############ Get Survey constants:
  # re-order by camera ID
  meta_t <- meta_t[order(meta_t$Site_N, meta_t$camera_on),]
  meta_t$Site_N_t <- 1:nrow(meta_t)
  
  E <- as.numeric(meta_t$camera_L) # total effort
  ncam <- nrow(meta) # number of cameras
  nhist <- nrow(meta_t) # number of camera histories
  nobs <- rep(NA, nhist)
  # nobs is the number of detections + 1
  # to include time between last detection and next
  for (i in 1:nhist){
    nobs[i] <- nrow(data_t[data_t$Site_ID_t==meta_t$Site_ID_t[i],])+1
  }
  
  #############################################################
  ############## Create observation matrix
  X <- matrix(NA, nrow = nhist, ncol = max(nobs)+1)
  type <- matrix(NA, nrow = nhist, ncol = max(nobs)+1)

  X[,1] <- -1
  type[,1:2] <- -1
  
  for (i in 1:nhist){ # loop through cameras
    site <- meta_t$Site_ID_t[i]
    df.site <- data_t[data_t$Site_ID_t == site,]
    if (nrow(df.site) == 0){ # if we have no observations
      X[i,2] <- E[i]
    } else { # otherwise make matrix of inter-detection times
      X[i,2:(nobs[i]+1)] <- as.numeric(difftime(c(df.site$Date.Time[1:(nobs[i]-1)],
                                              meta_t$camera_off[meta_t$Site_ID_t == site]),
                                            c(meta_t$camera_on[meta_t$Site_ID_t == site],
                                              df.site$Date.Time[1:(nobs[i]-1)]),
                                            units = 'days'))
      if (marks){ # if tracking marks
        if (nobs[i]>1){
          type[i,3:(nobs[i]+1)] <- df.site$Species.ind
        }
      } else {
        type <- NULL
      }
    }
  }
  
  if (any(X[,-1]<0, na.rm=T)){
    stop('In X: negative inter-detection times. Check detections and meta timings are compatible.')
  }
  
  #############################################################
  ####### Now to do state at collections
  start <- matrix(0, nrow = nhist, ncol = 3)
  collection = matrix(0, nrow = nhist, ncol = 3)

  if (is.null(intermediate.checks)){ # if no checks during survey
    start[,1] <- 1 # all cameras start working
    if (!state.at.collection){# if no state at collection recorded
      collection[,] <- 1
    } else {
      collection[,ifelse(is.na(meta_t$camera_fail), 1, 2)] <- 1
    }
  }
  
  if (!is.null(intermediate.checks)){
    for (i in 1:nhist){
      cam <- meta_t$Site_ID_t[i] # for this camera history
      if (!(meta_t$Site_ID[meta_t$Site_ID_t%in%cam] %in% unique(intermediate.checks$Site_ID))){ # if does not appear in checks, assume camera starts in state 1
        start[i,1] <- 1
        if (!state.at.collection){
          collection[i,1:3] <- 1
        } else {
          collection[i, ifelse(is.na(meta_t$camera_fail[i]), 1, 2)] <- 1
        }
      } else { # camera does appear in checks
        intermediate.checks.cam <- intermediate.checks[intermediate.checks$Site_ID==meta_t$Site_ID[meta_t$Site_ID_t%in%cam],]
        end.time.check <- meta_t$camera_off[meta_t$Site_ID_t==cam]
        if (any(intermediate.checks.cam$Date.Time==end.time.check)){ # end time corresponds to inter survey check
          which.check <- which(intermediate.checks.cam$Date.Time==end.time.check)
          state.at.check <- intermediate.checks.cam$state[which.check]
          collection[i,state.at.check] <- 1
        } else {
          which.check <- nrow(intermediate.checks.cam)+1 # end.time.check corresponds to end of survey
          if (state.at.collection){
            collection[i, ifelse(is.na(meta_t$camera_fail[i]), 1, 2)] <- 1
          } else {
            collection[i, 1:3] <- 1
          }
        }
        if (which.check == 1){ # if first check up
          start[i,1] <- 1 # camera always start in state 1
        }  else {
          if (is.na(intermediate.checks.cam$repaired[which.check-1]) | 
              (intermediate.checks.cam$repaired[which.check-1]==1)){ # if was repaired at previous check or working:
            start[i,1] <- 1
          } else if (intermediate.checks.cam$repaired[which.check-1]==0) {
            start[i,intermediate.checks.cam$state[which.check-1]] <- 1
          }
        }
      }
    }
  }
  
  return(list(ncam = ncam,
              nhist = nhist,
              nobs = nobs,
              cam.n = meta_t$Site_N,
              data_t = data_t,
              meta_t = meta_t,
              X = X, 
              E = E,
              type = type,  
              state.at.collection = state.at.collection,
              collection = collection,
              start = start,
              marks = marks, 
              fix_transitions = fix_transitions))
  
}

set_model <- function(MCMCdata, camera_hetero = T, 
                      priors = list(detection = list(ltp.scale = 1, ltp.shape = 1,
                                                     l1fp.scale = 1, l1fp.shape = 1,
                                                     l3fp.scale = 1, l3fp.shape = 200),
                                    transition = list(mu12.scale = 0.01, mu12.shape = 1,
                                                      mu13.scale = 0.01, mu13.shape = 1,
                                                      mu21.scale = 0.01, mu21.shape = 1,
                                                      mu32.scale = 1, mu32.shape = 1)),
                      initial.vals = list(detection = list(lambdatp = rep(priors$detection$ltp.scale*priors$detection$ltp.shape, MCMCdata$ncam),
                                                            lambda1fp = rep(priors$detection$l1fp.scale*priors$detection$l1fp.shape, MCMCdata$ncam),
                                                            lambda3fp = priors$detection$l3fp.scale*priors$detection$l3fp.shape),
                                          transition = list(mu12 = priors$transition$mu12.scale*priors$transition$mu12.shape,
                                                            mu13 = priors$transition$mu13.scale*priors$transition$mu13.shape,
                                                            mu21 = priors$transition$mu21.scale*priors$transition$mu21.shape,
                                                            mu32 = priors$transition$mu32.scale*priors$transition$mu32.shape)),
                      pars_to_track = NULL){
  
 ############################ check that priors and initial values match the code restrictions                      
                      
  prior.names <- names(priors$detection)
  init.names <- names(initial.vals$detection)
  
  if (MCMCdata$marks){ # if marks = T
    if (!all(c('lambdatp', 'lambda1fp', 'lambda3fp') %in% init.names)){
      stop('In initial.vals: expecting named detection values lambdatp, lambda1fp, lambda3fp')
    }
    if (!all(c('ltp.scale', 'ltp.shape',
                                'l1fp.scale', 'l1fp.shape',
                                'l3fp.scale', 'l3fp.shape') %in% prior.names)){
      stop('In priors: expecting named detection values ltp.scale, ltp.shape, l1fp.scale,
           l1fp.shape, l3fp.scale', 'l3fp.shape')
    }
    if (camera_hetero){
      n <- MCMCdata$ncam
      init.length <- c(length(initial.vals$detection$lambdatp),
                        length(initial.vals$detection$lambda1fp),
                        length(initial.vals$detection$lambda3fp))
      if(!all(init.length[1:2] == n)){
        warning('In initial.vals: expecting lambdatp and lambda1fp to have length n. Copying first value.')
        initial.vals$detection$lambdatp <- rep(initial.vals$detection$lambdatp[1], n)
        initial.vals$detection$lambda1fp <- rep(initial.vals$detection$lambda1fp[1], n)
      } 
      if (init.length[3] != 1){
        warning('In initial.vals: expecting lambda3fp to have length 1. Taking only first value')
        initial.vals$detection$lambda3fp <- initial.vals$detection$lambda3fp[1]
      }
    } else {
      n <- MCMCdata$ncam
      init.length <- c(length(initial.vals$detection$lambdatp),
                        length(initial.vals$detection$lambda1fp),
                        length(initial.vals$detection$lambda3fp))
      if(!all(init.length == 1)){
        warning('In initial.vals: expecting lambdatp, lambda1fp, lambda3fp to have length 1. Taking only first value.')
        initial.vals$detection$lambdatp <- initial.vals$detection$lambdatp[1]
        initial.vals$detection$lambda1fp <- initial.vals$detection$lambda1fp[1]
        initial.vals$detection$lambda3fp <- initial.vals$detection$lambda3fp[1]
      }
    }
  } else { # if marks = F
    if (!all(c('lambda1', 'lambda3') %in% init.names)){
      stop('In initial.vals: expecting named detection values lambda1, lambda3')
    }
    if (!all(c('l1.scale', 'l1.shape',
                                'l3.scale', 'l3.shape') %in% prior.names)){
      stop('In priors: expecting named detection values l1.scale, l1.shape, l3.scale, l3.shape')
    }
    if (camera_hetero){
      n <- MCMCdata$ncam
      init.length <- c(length(initial.vals$detection$lambda1),
                        length(initial.vals$detection$lambda3))
      if(!(init.length[1] == n)){
        warning('In initial.vals: expecting lambda1 to have length n. Copying first value.')
        initial.vals$detection$lambda1 <- rep(initial.vals$detection$lambda1[1], n)
      } 
      if (init.length[2] != 1){
        warning('In initial.vals: expecting lambda3 to have length 1. Taking only first value')
        initial.vals$detection$lambda3 <- initial.vals$detection$lambda3[1]
      }
    } else {
      n <- MCMCdata$ncam
      init.length <- c(length(initial.vals$detection$lambda1),
                       length(initial.vals$detection$lambda3))
      if(!all(init.length == 1)){
        warning('In initial.vals: expecting lambda1 and lambda3 to have length 1. Taking only first value.')
        initial.vals$detection$lambda1 <- initial.vals$detection$lambda1[1]
        initial.vals$detection$lambda3 <- initial.vals$detection$lambda3[1]
      }
    }
  }
  
  # Transition parameters:
  # Always need mu12 to be non.zero:
  if (initial.vals$transition$mu12 == 0){
    stop('In initial.vals: Expecting mu12 to be non-zero.')
  }
  if (priors$transition$mu12.shape==0 | priors$transition$mu12.scale==0){
    stop('in priors: Expecting mu12.shape and mu12.scale to be non-zero.')
  }
  
  if (MCMCdata$fix_transitions$mu21.0){
    if (initial.vals$transition$mu21!=0){
      stop('In initial.vals: Fixing mu21 to zero. Expecting mu21 to be zero.')
    }
    warning('Fixing mu21 to zero. Ignoring mu21.scale and mu21.shape if provided.')
  } else {
    if (priors$transition$mu21.shape == 0 | priors$transition$mu21.scale == 0 | 
        initial.vals$transition$mu21 == 0){
      stop('In priors or initial.vals: Expecting non-zero elements for mu21, mu21.shape, mu21.scale.')
    }
  }
  
  if (MCMCdata$fix_transitions$mu13.0){
    if (initial.vals$transition$mu13!=0){
      stop('In initial.vals: Fixing theta13 to zero. Expecting mu13 to be zero.')
    }
    warning('Fixing theta13 to zero. Ignoring mu13.scale and mu13.shape if provided. Also ignoring mu32, mu32.shape, mu32.scale if provided.')
  } else {
    if (initial.vals$transition$mu13==0 | priors$transition$mu13.scale==0 | priors$transition$mu13.shape==0){
      stop('In initial.vals or priors: Expecting non-zero elements for mu13, mu13.shape, mu13.scale.')
    }
    if (priors$transition$mu32.shape == 0 | priors$transition$mu32.shape == 0 | 
        initial.vals$transition$mu32 == 0){
      stop('In priors or initial.vals: Expecting non-zero elements for mu32, mu32.shape, mu32.scale.')
    }
  }
  
  ########### Set constants, data and inits ############
  
  DEconstants <- c(list(ncam = MCMCdata$ncam, nhist = MCMCdata$nhist,
                        cam.n = MCMCdata$cam.n,
                        n_inter_detections = MCMCdata$nobs,
                        start = MCMCdata$start, 
                        n.states = 3,
                        collection = MCMCdata$collection,
                        type = MCMCdata$type,
                        marks = MCMCdata$marks, camera_hetero = camera_hetero,
                        mu21.0 = MCMCdata$fix_transitions$mu21.0,
                        mu13.0 = MCMCdata$fix_transitions$mu13.0
                        ),
                   priors$transition,
                   priors$detection)
  
  DEdata <- list(X = MCMCdata$X)
  
  DEinits = function(){
    c(initial.vals$transition, initial.vals$detection)
  }
  
  if (is.null(pars_to_track)){
    if (MCMCdata$marks){
      if (camera_hetero){
        pars_to_track = c('lambdatp', 'lambda3fp', 'lambda1fp', 'mu12', 'mu13', 'mu21', 'mu32')
      } else {
        pars_to_track = c('lambdatp', 'lambda1fp', 'lambda3fp', 'mu12', 'mu13', 'mu21', 'mu32')
      }
    } else {
      if (camera_hetero){
        pars_to_track = c('lambda1', 'lambda3', 'mu12', 'mu13', 'mu21', 'mu32')
      } else {
        pars_to_track = c('lambda1', 'lambda3', 'mu12', 'mu13', 'mu21', 'mu32')
      }
    }
  }
  
  return(list(DEconstants = DEconstants, 
              DEdata = DEdata, 
              DEinits = DEinits, 
              pars_to_track = pars_to_track, 
              MCMCdata = MCMCdata))
  
}


############################ MCMC functions

dchist4 <- nimbleFunction(
  #Define the inputs
  run = function(x = double(1),
                 mu12 = double(0),
                 mu21 = double(0),
                 mu13 = double(0),
                 mu32 = double(0),
                 ltp = double(0),
                 l1fp = double(0),
                 l3fp = double(0),
                 n.states = integer(0),
                 type. = double(1),
                 collection = double(1),
                 start = double(1),
                 log = integer(0, default = 0)){
    
    returnType(double(0)) # specify scalar return
    
    d_i <- length(x) - 2 # number of detections
    is.di.0 <- 1*(d_i == 0)
    
    # dynamic scaling so all rates are below 1:
    all.rates <- numeric(length = 7)
    all.rates[1:7] <- c(mu12, mu21, mu13, mu32, ltp, l1fp, l3fp)
    max.rate <- max(all.rates[1:7])
    fctr <- ceiling(log(max.rate)/log(10)) + 1 # just to be safe
    # re-scale the rates in the model:
    mu12 <- mu12/(1*10**fctr)
    mu13 <- mu13/(1*10**fctr)
    mu21 <- mu21/(1*10**fctr)
    mu32 <- mu32/(1*10**fctr)
    ltp <- ltp/(1*10**fctr)
    l1fp <- l1fp/(1*10**fctr)
    l3fp <- l3fp/(1*10**fctr)
    
    x. <- numeric(length = d_i+1)
    x.[1:(d_i+1)] <- (x[2:(d_i+2)])*(1*10**fctr) # truncate leading -1
    
    s_factors <- 0 # sum of factors
    s_factors <- s_factors + (d_i*fctr)
    
    # redefine type to remove the starting 2
    type <- numeric(length = 1*(is.di.0) + (d_i)*(1-is.di.0))
    type[1:length((1*(is.di.0) + 2*(1-is.di.0)):(d_i+1))] <- type.[(1 + 1*(is.di.0) + 2*(1-is.di.0)):(d_i+2)]
    
    Q <- matrix(type = 'double', nrow = n.states, ncol = n.states)
    Q[1,1] <- -(mu12+mu13)
    Q[1,2] <- mu12
    Q[1,3] <- mu13
    Q[2,1] <- mu21
    Q[2,2] <- -mu21
    Q[2,3] <- 0
    Q[3,1] <- 0
    Q[3,2] <- mu32
    Q[3,3] <- -mu32
    
    l1 <- ltp + l1fp
    l2 <- 0
    l3 <- ltp + l3fp
    L <- matrix(type = 'double', nrow = n.states, ncol = n.states)
    L[1:n.states, 1:n.states] <- diag(c(l1, l2, l3))
    
    d_chist <- matrix(type = 'double', nrow = n.states, ncol = n.states)
    d_chist <- diag(3)
    e_CxL <- matrix(type = 'double', nrow = n.states, ncol = n.states)
    dt1 <- matrix(type = 'double', nrow = n.states, ncol = n.states) # proposal
    dt <- matrix(type = 'double', nrow = n.states, ncol = n.states) # product
    dt <- diag(3)
    
    # Go through detections
    if (d_i != 0){
      for (j in 1:d_i){ # loop through detections
        e_CxL[,] <- nmat.exp.pade((Q[,]-L[,])*x.[j], n.states = n.states) %*% 
          diag(c(type[j]*ltp + (1-type[j])*l1fp, 
                 0, 
                 type[j]*ltp + (1-type[j])*l3fp))
        dt1[,] <- dt[,] %*% e_CxL[,]
        if (mu21 != 0){ 
          if (mu13 != 0){
            if (any(dt1[,1] == 0) | any(dt1[,3] == 0)){ # if mu12, mu12 != 0
              fctr <- max(ceiling(log(dt[,])/log(10))) 
              s_factors <- s_factors + fctr
              d_chist[,] <- d_chist[,] %*% (dt[,]/(1*10**fctr))
              dt[,] <- e_CxL[,]
            } else {
              dt[,] <- dt1[,]
            }
          } else {
            if (any(dt1[1:2,1] == 0)){ # if mu12 !=0 but mu13 =0
              fctr <- max(ceiling(log(dt[,])/log(10))) 
              s_factors <- s_factors + fctr
              d_chist[,] <- d_chist[,] %*% (dt[,]/(1*10**fctr))
              dt[,] <- e_CxL[,]
            } else {
              dt[,] <- dt1[,]
            }
          }
        } else {
          if (mu13 != 0){
            if (dt1[1,1] == 0 | dt1[1,3] == 0){ # mu12 = 0 but mu13 !=0
              fctr <- max(ceiling(log(dt[,])/log(10))) 
              s_factors <- s_factors + fctr
              d_chist[,] <- d_chist[,] %*% (dt[,]/(1*10**fctr))
              dt[,] <- e_CxL[,]
            } else {
              dt[,] <- dt1[,]
            }
          } else {
            if (dt1[1,1] == 0){ # mu12, mu13 =0
              fctr <- max(ceiling(log(dt[,])/log(10))) 
              s_factors <- s_factors + fctr
              d_chist[,] <- d_chist[,] %*% (dt[,]/(1*10**fctr))
              dt[,] <- e_CxL[,]
            } else {
              dt[,] <- dt1[,]
            }
          }
        }
      }
    }
    
    # Handle last detection set (last set won't trigger the failure)
    fctr <- max(ceiling(log(dt[,])/log(10))) 
    s_factors <- s_factors + fctr
    d_chist[,] <- d_chist[,] %*% (dt[,]/(1*10**fctr))
    
    # Trailing no observation
    e_CxL[,] <- nmat.exp.pade((Q[,]-L[,])*x.[d_i+1], n.states = n.states)
    fctr <- max(ceiling(log(e_CxL[,])/log(10))) 
    s_factors <- s_factors + fctr
    d_chist[,] <- d_chist[,] %*% (e_CxL[,]/(1*10**fctr))
    
    # multiply all together
    d.dot <- start[1:3] %*% d_chist[,] %*% collection[1:3]
    d <- log(d.dot)/log(10) + s_factors  # in base 10
    d <- log(10) * d # in base e
    
    if (log==1) {return(d[1,1])} else {return(exp(d[1,1]))}
    
  }
)

dchist4.notype <- nimbleFunction(
  #Define the inputs
  run = function(x = double(1),
                 mu12 = double(0),
                 mu21 = double(0),
                 mu13 = double(0),
                 mu32 = double(0),
                 l1 = double(0),
                 l3 = double(0),
                 n.states = integer(0),
                 collection = double(1),
                 start = double(1),
                 n_inter_detections = double(0),
                 log = integer(0, default = 0)){
    
    returnType(double(0)) # specify scalar return
    
    #n_inter_detections is ignored in the density function
    
    d_i <- length(x) - 2 # number of detections
    is.di.0 <- 1*(d_i == 0)
    
    # dynamic scaling so all rates are below 1:
    all.rates <- numeric(length = 6)
    all.rates[1:6] <- c(mu12, mu21, mu13, mu32, l1, l3)
    max.rate <- max(all.rates[1:6])
    fctr <- ceiling(log(max.rate)/log(10)) + 1 # just to be safe
    # re-scale the rates in the model:
    mu12 <- mu12/(1*10**fctr)
    mu13 <- mu13/(1*10**fctr)
    mu21 <- mu21/(1*10**fctr)
    mu32 <- mu32/(1*10**fctr)
    l1 <- l1/(1*10**fctr)
    l3 <- l3/(1*10**fctr)
    
    x. <- numeric(length = d_i+1)
    x.[1:(d_i+1)] <- (x[2:(d_i+2)])*(1*10**fctr) # truncate leading -1
    
    s_factors <- 0 # sum of factors
    s_factors <- s_factors + (d_i*fctr)
    
    Q <- matrix(type = 'double', nrow = n.states, ncol = n.states)
    Q[1,1] <- -(mu12+mu13)
    Q[1,2] <- mu12
    Q[1,3] <- mu13
    Q[2,1] <- mu21
    Q[2,2] <- -mu21
    Q[2,3] <- 0
    Q[3,1] <- 0
    Q[3,2] <- mu32
    Q[3,3] <- -mu32
    
    l2 <- 0
    L <- matrix(type = 'double', nrow = n.states, ncol = n.states)
    L[1:n.states, 1:n.states] <- diag(c(l1, l2, l3))
    
    d_chist <- matrix(type = 'double', nrow = n.states, ncol = n.states)
    d_chist <- diag(3)
    e_CxL <- matrix(type = 'double', nrow = n.states, ncol = n.states)
    dt1 <- matrix(type = 'double', nrow = n.states, ncol = n.states) # proposal
    dt <- matrix(type = 'double', nrow = n.states, ncol = n.states) # product
    dt <- diag(3)
    
    # Go through detections
    if (d_i != 0){
      for (j in 1:d_i){ # loop through detections
        e_CxL[,] <- nmat.exp.pade((Q[,]-L[,])*x.[j], n.states = n.states) %*% L
        dt1[,] <- dt[,] %*% e_CxL[,]
        if (mu21 != 0){ 
          if (mu13 != 0){
            if (any(dt1[,1] == 0) | any(dt1[,3] == 0)){ # if mu12, mu12 != 0
              fctr <- max(ceiling(log(dt[,])/log(10))) 
              s_factors <- s_factors + fctr
              d_chist[,] <- d_chist[,] %*% (dt[,]/(1*10**fctr))
              dt[,] <- e_CxL[,]
            } else {
              dt[,] <- dt1[,]
            }
          } else {
            if (any(dt1[1:2,1] == 0)){ # if mu12 !=0 but mu13 =0
              fctr <- max(ceiling(log(dt[,])/log(10))) 
              s_factors <- s_factors + fctr
              d_chist[,] <- d_chist[,] %*% (dt[,]/(1*10**fctr))
              dt[,] <- e_CxL[,]
            } else {
              dt[,] <- dt1[,]
            }
          }
        } else {
          if (mu13 != 0){
            if (dt1[1,1] == 0 | dt1[1,3] == 0){ # mu12 = 0 but mu13 !=0
              fctr <- max(ceiling(log(dt[,])/log(10))) 
              s_factors <- s_factors + fctr
              d_chist[,] <- d_chist[,] %*% (dt[,]/(1*10**fctr))
              dt[,] <- e_CxL[,]
            } else {
              dt[,] <- dt1[,]
            }
          } else {
            if (dt1[1,1] == 0){ # mu12, mu13 =0
              fctr <- max(ceiling(log(dt[,])/log(10))) 
              s_factors <- s_factors + fctr
              d_chist[,] <- d_chist[,] %*% (dt[,]/(1*10**fctr))
              dt[,] <- e_CxL[,]
            } else {
              dt[,] <- dt1[,]
            }
          }
        }
      }
    }
    
    # Handle last detection set (last set won't trigger the failure)
    fctr <- max(ceiling(log(dt[,])/log(10))) 
    s_factors <- s_factors + fctr
    d_chist[,] <- d_chist[,] %*% (dt[,]/(1*10**fctr))
    
    # Trailing no observation
    e_CxL[,] <- nmat.exp.pade((Q[,]-L[,])*x.[d_i+1], n.states = n.states)
    fctr <- max(ceiling(log(e_CxL[,])/log(10))) 
    s_factors <- s_factors + fctr
    d_chist[,] <- d_chist[,] %*% (e_CxL[,]/(1*10**fctr))
    
    # multiply all together
    d.dot <- start[1:3] %*% d_chist[,] %*% collection[1:3]
    d <- log(d.dot)/log(10) + s_factors  # in base 10
    d <- log(10) * d # in base e
    
    if (log==1) {return(d[1,1])} else {return(exp(d[1,1]))}
    
  }
)

rchist4 <- nimbleFunction(
  #Define the inputs
  run = function(n = integer(0),
                 mu12 = double(0),
                 mu21 = double(0),
                 mu13 = double(0),
                 mu32 = double(0),
                 ltp = double(0),
                 l1fp = double(0),
                 l3fp = double(0),
                 n.states = integer(0),
                 type. = double(1),
                 collection = double(1),
                 start = double(1)){
    # collection is how we handle the state at collection of camera
    
    returnType(double(1)) # specify scalar return
    
    d <- length(type.) + 2
    
    hist <- rexp(d, rate = ltp)
    
    return(hist)
    
  }
)

rchist4.notype <- nimbleFunction(
  #Define the inputs
  run = function(n = integer(0),
                 mu12 = double(0),
                 mu21 = double(0),
                 mu13 = double(0),
                 mu32 = double(0),
                 l1 = double(0),
                 l3 = double(0),
                 n.states = integer(0),
                 collection = double(1),
                 start = double(1),
                 n_inter_detections = double(0)){
    # collection is how we handle the state at collection of camera

    returnType(double(1)) # specify scalar return

    d <- n_inter_detections + 2

    hist <- rexp(d, rate = l1)

    return(hist)

  }
)

DEcode = nimbleCode({
  
  ############# Priors Transition Rates
  mu12 ~ dgamma(shape = mu12.shape, scale = mu12.scale) 
  
  if (mu21.0 == TRUE){
    mu21 <- 0
  } else {
    mu21 ~ dgamma(shape = mu21.shape, scale = mu21.scale)
  }
  if (mu13.0 == TRUE){
    mu13 <- 0
    mu32 <- 0
  } else {
    mu13 ~ dgamma(shape = mu13.shape, scale = mu13.scale)
    mu32 ~ dgamma(shape = mu32.shape, scale = mu32.scale)
  }
  
  ############### Priors detection rates
  if (marks){
    
    if (camera_hetero){
      for (i in 1:ncam){
        lambdatp[i] ~ dgamma(shape = ltp.shape, scale = ltp.scale) 
        lambda1fp[i] ~ dgamma(shape = l1fp.shape, scale = l1fp.scale)
        lambdatp.[i] <- lambdatp[i]
        lambda1fp.[i] <- lambda1fp[i]
      }
    } else {
      lambdatp ~ dgamma(shape = ltp.shape, scale = ltp.scale)
      lambda1fp ~ dgamma(shape = l1fp.shape, scale = l1fp.scale)
      for (i in 1:ncam){
        lambdatp.[i] <- lambdatp
        lambda1fp.[i] <- lambda1fp
      }
    }
    lambda3fp ~ dgamma(shape = l3fp.shape, scale = l3fp.scale)
    
  } else {
    
    if (camera_hetero){
      for (i in 1:ncam){
        lambda1[i] ~ dgamma(shape = l1.shape, scale = l1.scale) 
        lambda1.[i] <- lambda1[i]
      }
    } else {
      lambda1 ~ dgamma(shape = l1.shape, scale = l1.scale)
      for (i in 1:ncam){
        lambda1.[i] <- lambda1
      }
    }
    lambda3 ~ dgamma(shape = l3.shape, scale = l3.scale)
    
  }

  # density of observations
  for (i in 1:nhist){ # loop through the cameras
    if (marks){
         X[i,1:(n_inter_detections[i]+1)] ~ dchist4(mu12 = mu12, mu21 = mu21, mu13 = mu13, mu32 = mu32,
                                               ltp = lambdatp.[cam.n[i]],
                                               l1fp = lambda1fp.[cam.n[i]],
                                               l3fp = lambda3fp,
                                               n.states = n.states,
                                               type. = type[i,1:(n_inter_detections[i]+1)],
                                               collection = collection[i, 1:3],
                                               start = start[i, 1:3])
    } else {
      X[i,1:(n_inter_detections[i]+1)] ~ dchist4.notype(mu12 = mu12, mu21 = mu21, mu13 = mu13, mu32 = mu32,
                                                 l1 = lambda1.[cam.n[i]],
                                                 l3 = lambda3,
                                                 n.states = n.states,
                                                 collection = collection[i, 1:3],
                                                 start = start[i, 1:3],
                                                 n_inter_detections = n_inter_detections[i])
    }
  }
  
})

# MCMCsetup is output from set_model
run_MCMC <- function(MCMCsetup, niter = 110000, nburnin=10000, thin = 100, nchains = 3){
  
  list2env(MCMCsetup, environment())
  
  DEmodel = nimbleModel(DEcode, constants = DEconstants, data = DEdata, inits = DEinits())
  
  # Check the initial conditions for problems:
  logProbs.inits <- rep(NA, MCMCsetup$DEconstants$nhist)
  for (i in 1:MCMCsetup$DEconstants$nhist){
    logProbs.inits[i] <- DEmodel$calculate(paste0('X[', as.character(i), ', 1:', as.character(MCMCsetup$DEconstants$n_inter_detections[i]+1) ,']'))
  }
  if (any(logProbs.inits== -Inf)){
    affected_camera_hist <- MCMCsetup$MCMCdata$meta_t$Site_ID_t[which(logProbs.inits== -Inf)]
    affected_cameras <- MCMCsetup$MCMCdata$meta_t$Site_ID[which(logProbs.inits== -Inf)]
    message <- cat(paste0('logProbs for following cameras: ', paste0(affected_cameras, collapse = ', '), ' are -Inf.\n',
    'This indicates either problems with initial values or that detections at cameras violate modelling assumptions.\n',
    'Recommend checking assumptions, initial values, or removing affected cameras.\n',
    'MCMC may fail to converge if model not correctly initialised.\n', 
    "Would you like to continue? (Y/N)"))
    question1 <- readline(message)
    if(regexpr(question1, 'n', ignore.case = TRUE) == 1){
      stop('In setting up model: Checking initialisation.')
    }  
  }
  
  mcmcConf = configureMCMC(DEmodel, monitors = pars_to_track, print = TRUE)
  DEmcmc <- buildMCMC(mcmcConf, print = FALSE)
  
  cDEmodel <- compileNimble(DEmodel)
  cDEmcmc <- compileNimble(DEmcmc, project = cDEmodel)
  
  DEresults <- runMCMC(cDEmcmc, niter = niter, nburnin = nburnin, thin = thin, nchains = nchains)
  
  return(DEresults)
  
}

############### Plot results ################
# these are the equivalents to the plotting functions for the HMM

# View results through time
plot_time <- function(results, state = 1, MCMCdata, cameras = NULL, save.plots = F, res=100,
                      ...){
  
  # cameras is named cameras (if NULL all cameras are plotted)
  # results needs to be output from run_MCMC
  meta_t <- MCMCdata$meta_t
  df_t <- MCMCdata$data_t
  
  meta_t <- meta_t[order(meta_t$Site_N, meta_t$camera_on),]
  
  if (is.null(cameras)){
    # what cameras do we want to see (by name in meta_t$Site_ID)
    cameras <- unique(meta_t$Site_ID)
    to_plot <- unique(meta_t$Site_ID)
  } else {
    to_plot <- cameras
  }
  
  # concatenate the chains
  n.chains <- ifelse(is.list(results), length(results), 1)
  
  if (n.chains > 1){
    samples <- results[[1]]
    for (i in 2:n.chains){
      samples <- rbind(samples, results[[i]]) # merge the entries
    }
  } else {
    samples <- results
  }
  
  # save posterior state results
  post_samples <- nrow(samples)
  states_post <- array(dim = c(length(to_plot), post_samples, res))
  
  plot.list <- list()
  
  for (i in 1:length(to_plot)){# loop through cameras
    
    cam <- to_plot[i]
    
    df.i <- df_t[df_t$Site_ID %in% cam, ]
    df.i <- df.i[!(df.i$Camera_Installation %in% 'Work'),]
    
    # plot the detections
    if (MCMCdata$marks){
      df.i$Species.ind <- as.factor(df.i$Species.ind)
      p1 <- ggplot(df.i) + geom_point(aes(x = Date.Time, y = 0, col = Species.ind)) + 
        ggtitle(cam) + xlim(min(meta_t$camera_on[meta_t$Site_ID %in% cam]), 
                            max(meta_t$camera_off[meta_t$Site_ID %in% cam])) +
        theme(legend.position = "top")
    } else {
      p1 <- ggplot(df.i) + geom_point(aes(x = Date.Time, y = 0)) + 
        ggtitle(cam) + xlim(min(meta_t$camera_on[meta_t$Site_ID %in% cam]), 
                            max(meta_t$camera_off[meta_t$Site_ID %in% cam]))
    }
    
    # multiple failures are possible for each camera
    if (!all(is.na(meta_t$camera_fail[meta_t$Site_ID %in% cam]))){
      dt = data.frame(max = meta_t$camera_fail[(meta_t$Site_ID %in% cam) & !is.na(meta_t$camera_fail)])
      p1 <- p1 + geom_vline(data = dt, mapping = aes(xintercept = max), colour = 'red')
    }
    
    p1 <- p1 + theme(axis.text.y=element_blank(),
                     axis.ticks.y=element_blank()) + ylab('')
    
    # posteriors for this camera, which may be broken into distinct histories
    n.hist <- sum(meta_t$Site_ID %in% cam) # how many histories for this camera
    for (hist in 1:n.hist){
      
      hist_on <- as.numeric(meta_t$camera_on[meta_t$Site_ID_t %in% paste0(
        cam, ifelse(n.hist==1, '', paste0('.', as.character(hist)))
      )])
      hist_off <- as.numeric(meta_t$camera_off[meta_t$Site_ID_t %in% paste0(
        cam, ifelse(n.hist==1, '', paste0('.', as.character(hist)))
      )])
      camera_on <- as.numeric(min(meta_t$camera_on[meta_t$Site_ID %in% cam]))
      camera_off <- as.numeric(max(meta_t$camera_off[meta_t$Site_ID %in% cam]))
      res_hist <- sum(hist_on <= seq(camera_on, camera_off, length.out=res) &
                        seq(camera_on, camera_off, length.out=res) <= hist_off)
      
      # what camera_history numerber are we looking at
      site <- meta_t$Site_N_t[meta_t$Site_ID %in% cam][hist]
      
      # compute the state posteriors for this camera history
      x <- c(0, cumsum(MCMCdata$X[site, 2:(1+MCMCdata$nobs[site])]))
      
      if (MCMCdata$marks){
        if (MCMCdata$nobs[site]==1){
          type <- c(NA, NULL, NA)
        } else {
          type <- c(NA,  MCMCdata$type[site, 3:(1+MCMCdata$nobs[site])], NA)
        }
        # are we dealing with camera_hetero?
        camera_hetero <- sum(unlist(str_split(colnames(results), pattern='\\[')) %in% 'lambdatp') > 1
        if (camera_hetero){
          states_post_hist <- cprob_t(x, type, res_hist, 
                                      results[,paste0('lambdatp[', as.character(site), ']')],
                                      results[,paste0('lambda1fp[', as.character(site), ']')],
                                      results[,'lambda3fp'],
                                      results[, 'mu12'],
                                      results[, 'mu13'],
                                      results[, 'mu21'],
                                      results[, 'mu32'],
                                      3, state, MCMCdata$E[site], 
                                      MCMCdata$start[site,], MCMCdata$collection[site,])
        } else {
          states_post_hist <- cprob_t(x, type, res_hist, 
                                      results[,paste0('lambdatp')],
                                      results[,paste0('lambda1fp')],
                                      results[,'lambda3fp'],
                                      results[, 'mu12'],
                                      results[, 'mu13'],
                                      results[, 'mu21'],
                                      results[, 'mu32'],
                                      3, state, MCMCdata$E[site], 
                                      MCMCdata$start[site,], MCMCdata$collection[site,])
        }
      } else {
        # are we dealing with camera_hetero?
        camera_hetero <- sum(unlist(str_split(colnames(results), pattern='\\[')) %in% 'lambdatp') > 1
        if (camera_hetero){
          states_post_hist <- cprob_t.notype(x, res_hist, 
                                             results[,paste0('lambda1[', as.character(site), ']')],
                                             results[,paste0('lambda3')],
                                             results[, 'mu12'],
                                             results[, 'mu13'],
                                             results[, 'mu21'],
                                             results[, 'mu32'],
                                             3, state, MCMCdata$E[site],
                                             MCMCdata$start[site,], MCMCdata$collection[site,])
        } else {
          states_post_hist <- cprob_t.notype(x, res_hist, 
                                             results[,paste0('lambda1')],
                                             results[,paste0('lambda3')],
                                             results[, 'mu12'],
                                             results[, 'mu13'],
                                             results[, 'mu21'],
                                             results[, 'mu32'],
                                             3, state, MCMCdata$E[site],
                                             MCMCdata$start[site,], MCMCdata$collection[site,])
        }
      }
      
      hist.start <- sum(!is.na(states_post[i,1,]))+1
      hist.end <- hist.start+res_hist-1
      states_post[i, , hist.start:hist.end] <- states_post_hist
      
    }
    
    m.states <- colMeans(states_post[i,,])
    ucl.states <- apply(states_post[i,,], 2, quantile, probs = 0.975)
    lcl.states <- apply(states_post[i,,], 2, quantile, probs = 0.025)
    
    df.i <- data.frame(times = seq(camera_on, camera_off, length.out=res), m.states = m.states,
                       ucl.states = ucl.states, lcl.states = lcl.states)
    df.i$times <- as.POSIXct(df.i$times, format = '%Y-%m-%d %H:%M:%OS', tz ='UTZ')
    
    p2 <- ggplot(data = df.i) + geom_line(aes(x = times, y = m.states)) + 
      geom_line(aes(x=times, y=ucl.states), linetype = 'dashed') + 
      geom_line(aes(x=times, y=lcl.states), linetype = 'dashed') + 
      ylim(0,1.01) + 
      theme(axis.text.y=element_blank(),
            axis.ticks.y=element_blank()) + ylab('')
    
    p <- grid.arrange(p1,p2, ncol=1)
    
    plot.list[[i]] <- p
    
    if (save.plots){
      ggsave(p, filename = paste0(meta$Site_ID[i], '.png'), ...)
    }
  }
  
  dimnames(states_post)[[1]] <- cameras
  
  #################### overall use plot:
  
  post_state_total <- apply(states_post, c(1,2), mean)
  m.state_total <- apply(post_state_total, 1, mean)
  ucl.state_total <- apply(post_state_total, 1, quantile, probs = 0.975)
  lcl.state_total <- apply(post_state_total, 1, quantile, probs = 0.025)

  df <- data.frame(matrix(nrow=length(cameras), ncol=4))
  colnames(df) <- c('Mean', 'UCL', 'LCL', 'Camera')
  
  df$Camera <- cameras
  df$Mean <- m.state_total
  df$LCL <- lcl.state_total
  df$UCL <- ucl.state_total
  
  df$Camera <- as.factor(df$Camera)
  
  p <- ggplot(df, aes(x = Camera, y = Mean)) + 
    geom_point() + geom_errorbar(aes(ymin=LCL, ymax=UCL)) + 
    ylab(paste0('Time in state ', as.character(state))) +
    theme(axis.text.x = element_text(angle=90, vjust=.5, hjust=1))+
    ylim(0,1)
  
  use.plot <- p
  
  if (save.plots){
    ggsave(p, filename = paste0(meta$Site_ID[i], '.png'), ...)
  }
  
  return(list(Posteriors = states_post, 
              Plots_time = plot.list, 
              Plot_proportion = use.plot))
}

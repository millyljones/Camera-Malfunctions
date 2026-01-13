library(nimble)
library(ggplot2)
library(mcmcplots)

build_data <- function(data, meta, 
                       T.unit = 1, 
                       fix_transitions = list(theta21.0 = F, theta13.0 = F),
                       state.at.collection = T,
                       intermediate.checks = NULL,
                       marks = T
){
  
  ############################################################
  ############ Want to remove any 'Work' entries from data
  data <- data[!(data$Camera_Installation %in% 'Work'), ]
  
  ############################################################
  ############ Order camera numbers by those in meta
  meta <- meta[order(meta$Site_N),]
  
  #############################################################
  ############ Get number of time units for survey:
  E <- meta$camera_L # total effort
  n.units <- E %/% T.unit # whole time units
  left.units <- E %% T.unit 
  fractional.left <- left.units / T.unit # fractional time units
  n.units <- n.units + 1*(left.units > 0) # total time units
  
  #############################################################
  ############## Create observation matrix
  n = nrow(meta) # number of cameras
  if (marks){
    X <- array(dim = c(n, max(n.units), 2))
  } else {
    X <- matrix(nrow = n, ncol = max(n.units))
  }
  
  for (i in 1:n){ # loop through the cameras
    
    ID = meta$Site_ID[i]
    df.i <- data[data$Site_ID==ID,] # get this camera's data
    # take time difference between each detections and start of survey
    df.i$Date.Time <- difftime(df.i$Date.Time, meta$camera_on[i], units = 'days') 
    
    for (j in 1:n.units[i]){ # loop through time units
      t.from <- (j-1)*T.unit
      t.to <- j*T.unit
      df.i.t <- df.i[df.i$Date.Time <= t.to & df.i$Date.Time >= t.from,] # all detections within time frame
      if (marks){
        X[i,j,1] <- sum(df.i.t$Species.ind == 1) # number of true positive
        X[i,j,2] <- sum(df.i.t$Species.ind == 0) # number of false positive
      } else {
        X[i,j] <- nrow(df.i.t) # total number of detections
      }
    }
  }
  
  #############################################################
  ##### Create inputs for the MCMC
  
  # For each time-unit, is this a whole unit or fractional? (to account for 
  # trailing half-days at the end of survey periods)
  frac <- matrix(nrow=n, ncol = max(n.units))
  for (i in 1:n){
    for (j in 1:(n.units[i] - 1)){
      frac[i,j] <- 1
    }
    frac[i,n.units[i]] <- ifelse(fractional.left[i] == 0, 1, fractional.left[i])
  }
  
  # Fill in any known states during the survey period
  # we always assume that camera starts by functioning normally so state[i,1] = 1
  state <- matrix(nrow = nrow(X), ncol = ncol(X))
  state[,1] <- 1 # known
  
  # if we know state at collection (info contained in meta)
  if (state.at.collection){
    for (i in 1:nrow(state)){
      state[i,n.units[i]] <- 1*(is.na(meta$camera_fail[i])) + 
        2*(!is.na(meta$camera_fail[i]))
      if (state[i, n.units[i]] ==2 & (ifelse(marks, sum(X[i,n.units[i],]), X[i,n.units[i]]) != 0)){
        state[i, n.units[i]] <- NA
        warning(paste0('Camera', as.character(meta$Site_ID[i]), ': Camera listed as broken at final time point but
                       detections made on this unit. Ignoring camera state and assuming camera state is unknown. 
                       Consider reducing size of time units.'))
      }
    }
  }
  
  # if we have any other information about camera functionality:
  # intermediate.checks is a matrix with cols Site_ID, Date.Time, state, and repaired
  if (!is.null(intermediate.checks) & fix_transitions$theta21.0 == T){
    if (any(intermediate.checks$repaired == 1)){
      stop('In intermediate.checks: fixing theta21 to 0, but showing cameras being repaired.')
    }
  }
  
  if (!is.null(intermediate.checks)){
    n.checks <- nrow(intermediate.checks) # how many checks
    for (i in 1:n.checks){
      n.i <- meta$Site_N[meta$Site_ID == intermediate.checks$Site_ID[i]] # what camera?
      day.i <- as.numeric(difftime(intermediate.checks$Date.Time[i], meta$camera_on[meta$Site_N == n.i], 
                                   units = 'days')) # what day?
      j <- (day.i %/% T.unit) + 1 # what day in analysis?
      state.i <- intermediate.checks$state[i] # state for this day
      repaired.i <- intermediate.checks$repaired[i]
      if (state.i == 1){ # if camera working at check up
        state[n.i, j] <- 1
      } else if (state.i == 2){ # if camera was broken at check up
        if (repaired.i == 0){ # and was not repaired
          state[n.i, j] <- 2
          if (ifelse(marks, sum(X[n.i, j,]), X[n.i, j]) != 0){
            warning(paste0('In intermediate checks: camera ', as.character(intermediate.checks$Site_ID[i]), ' listed as 
                           broken at time point ', as.character(j), ' but detections made on this
                           unit. Ignoring this entry in intermediate checks, camera state for this unit is assumed unknown. 
                           Consider reducing size of time units.'))
            state[n.i, j] <- NA
          }
        } else { # if camera was repaired
          prev.dets <- ifelse(marks, sum(X[n.i, j-1,]), X[n.i, j-1])
          state[n.i, j] <- 1 
          frac[n.i, j] <- 1 - (day.i%%T.unit)/(T.unit) # fractional working time
          if (prev.dets == 0 & is.na(state[n.i, j-1])){ # and no detections made the previous day
            state[n.i, j-1] <- 2 # previous day broken
          } 
        }
      } else if (state.i == 3){ # if camera fritzing at check up
        state[n.i, j] <- 3
        warning('Observed state 3 in intermediate checks: 
                may  have warnings about initialisation of MCMC chains. Currently working on
                solutions for this.')
      }
    }
  }
  
  if (fix_transitions$theta13.0 == T){
    # check no states are 3
    if (any(state == 3, na.rm=T)){
      stop('In state: theta13 is fixed to zero, but state matrix has 3`s.')
    }
  }
  
  if (fix_transitions$theta21.0 == T){
    # check no states go from 2 to 1
    for (i in 1:nrow(state)){
      known.states <- state[i, 1:n.units[i]][!is.na(state[i, 1:n.units[i]])]
      if (any(known.states==2)){
        if (any(known.states[-(1:(which(known.states==2)[1]))] == 1)){
          stop('In state: Fixed theta21 to zero, but state has transitions 2 to 1.')
        }
      }
    }
  }
  
  # Every day we don't have info on assume camera is working as initial condition
  state.init <- matrix(nrow = n, ncol = max(n.units))
  for (i in 1:nrow(state)){
    state.init[i,1:n.units[i]] <- 1
    if (fix_transitions$theta21.0 == T){
      # after any 2's must have all the rest be 2's
      if (any(state[i,1:n.units[i]] == 2, na.rm=T)){
        unit2 <- which(state[i,1:n.units[i]]==2)[1]
        state.init[i, min((unit2+1), n.units[i]):n.units[i]] <- 2
      }
    }
    state.init[!is.na(state.init) & !(is.na(state))] <- NA # NA's where we do have observations
  }
  
  return(list(n = n, n.units = n.units, X = X, state = state, state.init = state.init, frac = frac,
              marks = marks, fix_transitions = fix_transitions))
  
}

# MCMCdata is output from build_data
set_model <- function(MCMCdata, camera_hetero = T, 
                      priors = list(detection = list(ltp.scale = 1, ltp.shape = 1,
                                                     l1fp.scale = 1, l1fp.shape = 1,
                                                     l3fp.scale = 1, l3fp.shape = 200),
                                    transition = list(Alpha = matrix(c(0.99, 0.005, 0.005, 
                                                                       0.1, 0.9, 0,
                                                                       0, 0.9, 0.1), byrow = T, ncol = 3))),
                      initial.vals = list(detection = list(lambdatp = rep(1, MCMCdata$n),
                                                           lambda1fp = rep(1, MCMCdata$n),
                                                           lambda3fp = 800),
                                          transition = list(Theta = matrix(c(0.99, 0.005, 0.005,
                                                                              0.1, 0.9, 0,
                                                                              0, 0.9, 0.1), byrow=T, ncol=3))),
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
      n <- MCMCdata$n
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
      n <- MCMCdata$n
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
      n <- MCMCdata$n
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
      n <- MCMCdata$n
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
  if (!all(rowSums(initial.vals$transition$Theta) == 1)){
    stop('In initial.vals: rows of Theta do not sum to 1')
  }
  if (initial.vals$transition$Theta[2,3]!=0 | initial.vals$transition$Theta[3,1]!=0 | 
      priors$transition$Alpha[2,3]!=0 | priors$transition$Alpha[3,1]!=0){
    stop('In initial.vals or priors: Expect Alpha[2,3], Alpha[3,1], Theta[2,3], and Theta[3,1] to be zero.')
  }
  
  if (MCMCdata$fix_transitions$theta21.0){
    if (initial.vals$transition$Theta[2,1]!=0){
      stop('In initial.vals: Expecting Theta[2,1] to be zero.')
    }
    warning('Fixing theta21 to zero. Ignoring Alpha[2,] and Theta[2,].')
  } else {
    if (priors$transition$Alpha[2,2] == 0 | priors$transition$Alpha[2,1] == 0 | 
        initial.vals$transition$Theta[2,2] == 0 | initial.vals$transition$Theta[2,1] == 0){
      stop('In priors or initial.vals: Expecting non-zero elements for Alpha[2,2] and Alpha[2,1], and Theta[2,2] and Theta[2,1].')
    }
  }
  
  if (MCMCdata$fix_transitions$theta13.0){
    if (initial.vals$transition$Theta[1,3]!=0 | priors$transition$Alpha[1,3]!=0){
      stop('In initial.vals or priors: Fixing theta31 to zero. Expecting Theta[1,3] and Alpha[1,3] to be zero.')
    }
    warning('Fixing theta13 to zero. Alpha[3,] and Theta[3,] will have no effect on model.')
  } else {
    if (initial.vals$transition$Theta[1,3]==0 | priors$transition$Alpha[1,3]==0){
      stop('In initial.vals or priors: Expecting Theta[1,3] and Alpha[1,3] to be non zero.')
    }
    if (priors$transition$Alpha[3,2] == 0 | priors$transition$Alpha[3,3] == 0 | 
        initial.vals$transition$Theta[3,2] == 0 | initial.vals$transition$Theta[3,3] == 0){
      stop('In priors or initial.vals: Expecting non-zero elements for Alpha[3,2] and Alpha[3,3], and Theta[3,2] and Theta[3,3].')
    }
  }
  
  # transform the transition parameters to what is needed for the code
  Alpha = priors$transition$Alpha
  Theta = initial.vals$transition$Theta
  priors$transition = list(alpha1 = Alpha[1,][Alpha[1,]!=0],
                           alpha2 = Alpha[2,][Alpha[2,]!=0],
                           alpha3 = Alpha[3,][Alpha[3,]!=0])
  initial.vals$transition = list(theta1 = Theta[1,][Theta[1,]!=0],
                           theta2 = Theta[2,][Theta[2,]!=0],
                           theta3 = Theta[3,][Theta[3,]!=0])
  
  DEconstants <- c(list(n = MCMCdata$n, n.units = MCMCdata$n.units, 
                        marks = MCMCdata$marks, camera_hetero = camera_hetero,
                        theta21.0 = MCMCdata$fix_transitions$theta21.0,
                        theta13.0 = MCMCdata$fix_transitions$theta13.0,
                        frac = MCMCdata$frac),
                   priors$transition,
                   priors$detection)
  
  DEdata <- list(X = MCMCdata$X, state = MCMCdata$state)
  
  DEinits = function(){
    c(list(state = MCMCdata$state.init), initial.vals$transition, initial.vals$detection)
  }
  
  if (is.null(pars_to_track)){
    if (MCMCdata$marks){
      if (camera_hetero){
        pars_to_track = c('state', 'Theta', 'lambdatp', 'lambda3fp', 'lambda1fp', 'U1', 'U2', 'U3')
      } else {
        pars_to_track = c('state', 'Theta', 'lambdatp.', 'lambda1fp.', 'lambda3fp', 'U1', 'U2', 'U3')
      }
    } else {
      if (camera_hetero){
        pars_to_track = c('state', 'Theta', 'lambda1', 'lambda3', 'U1', 'U2', 'U3')
      } else {
        pars_to_track = c('state', 'Theta', 'lambda1.', 'lambda3', 'U1', 'U2', 'U3')
      }
    }
  }
  
  return(list(DEconstants = DEconstants, 
              DEdata = DEdata, 
              DEinits = DEinits, 
              pars_to_track = pars_to_track))
  
}

DEcode = nimbleCode({
  
  ############ Priors: Detection rates
  
  if (marks){
    
    if (camera_hetero){
      for (i in 1:n){
        lambdatp[i] ~ dgamma(shape = ltp.shape, scale = ltp.scale) 
        lambda1fp[i] ~ dgamma(shape = l1fp.shape, scale = l1fp.scale)
        lambdatp.[i] <- lambdatp[i]
        lambda1fp.[i] <- lambda1fp[i]
      }
    } else {
      lambdatp ~ dgamma(shape = ltp.shape, scale = ltp.scale)
      lambda1fp ~ dgamma(shape = l1fp.shape, scale = l1fp.scale)
      lambdatp.[1:n] <- rep(lambdatp, n)
      lambda1fp.[1:n] <- rep(lambdafp, n)
    }
    lambda3fp ~ dgamma(shape = l3fp.shape, scale = l3fp.scale)
    
  } else {
    
    if (camera_hetero){
      for (i in 1:n){
        lambda1[i] ~ dgamma(shape = l1.shape, scale = l1.scale) 
        lambda1.[i] <- lambda1[i]
      }
    } else {
      lambda1 ~ dgamma(shape = l1.shape, scale = l1.scale)
      lambda1.[1:n] <- rep(lambda1, n)
    }
    lambda3 ~ dgamma(shape = l3.shape, scale = l3.scale)
    
  }
  
  ############ Priors: Transition rates
  
  if (theta21.0 == T){
    theta2[1:2] <- c(0,1)
  } else {
    theta2[1:2] ~ ddirch(alpha = alpha2[1:2])
  }
  
  if (theta13.0 == T){
    theta1[1:2] ~ ddirch(alpha = alpha1[1:2])
    theta3[1:2] <- c(0,1)
    Theta[1,1:3] <- c(theta1[1], theta1[2], 0)
  } else {
    theta1[1:3] ~ ddirch(alpha = alpha1[1:3]) 
    theta3[1:2] ~ ddirch(alpha = alpha3[1:2])
    Theta[1,1:3] <- c(theta1[1], theta1[2], theta1[3])
  }
  
  Theta[2,1:3] <- c(theta2[1], theta2[2], 0)
  Theta[3,1:3] <- c(0, theta3[1], theta3[2])
  
  ############### Latent camera states
  
  for (i in 1:n){
    state[i,1] <- 1 
    for (j in 2:n.units[i]){
      state[i,j] ~ dcat(Theta[state[i,j-1],1:3])
    }
  }
  
  ############### Observation Probabilities
  
  for (i in 1:n){ # loop through cameras
    for (j in 1:n.units[i]){ # loop through time units
      if (marks){
        
        X[i,j,1] ~ dpois(frac[i,j] * (0 + lambdatp.[i] * (state[i,j]!=2)))
        X[i,j,2] ~ dpois(frac[i,j] * (0 + (lambda1fp.[i] * (state[i,j]==1)) + 
                                        (lambda3fp * (state[i,j]==3))))
      } else {
        
        X[i,j] ~ dpois(frac[i,j] * (0 + (lambda1.[i] * (state[i,j] == 1)) + 
                                      lambda3 * (state[i,j] == 3)))
      }
    }
  }
  
  ############### Calculate use per camera
  for (i in 1:n){
    U1[i] <- sum(state[i, 1:n.units[i]] == 1)/n.units[i]
    U2[i] <- sum(state[i, 1:n.units[i]] == 2)/n.units[i]
    U3[i] <- sum(state[i, 1:n.units[i]] == 3)/n.units[i]
  }
  
})

# MCMCsetup is output from set_model
run_MCMC <- function(MCMCsetup, niter = 110000, nburnin=10000, thin = 100, nchains = 3){

  list2env(MCMCsetup, environment())
  
  DEmodel = nimbleModel(DEcode, constants = DEconstants, data = DEdata, inits = DEinits())
  
  mcmcConf = configureMCMC(DEmodel, monitors = pars_to_track, print = TRUE)
  DEmcmc <- buildMCMC(mcmcConf, print = FALSE)
  
  cDEmodel <- compileNimble(DEmodel)
  cDEmcmc <- compileNimble(DEmcmc, project = cDEmodel)
  
  DEresults <- runMCMC(cDEmcmc, niter = niter, nburnin = nburnin, thin = thin, nchains = nchains)
  
  return(DEresults)

}

################################## Plots to view outputs

# View results through time
plot_time <- function(results, state = 2, df, meta, MCMCdata, cameras = NULL, save.plots = F, ...){
  
  # cameras is named cameras
  # results needs to be a list
  
  meta <- meta[order(meta$Site_N),]
  
  if (!is.null(cameras)){
    # what cameras do we want to see (by number)
    to_plot <- meta$Site_N[meta$Site_ID %in% cameras]
  } else {
    to_plot <- c(1:nrow(meta))
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
  
  n.units <- MCMCdata$n.units
      
  for (i in to_plot){
    
    df.i <- df[df$Site_ID == meta$Site_ID[i], ]
    df.i <- df.i[!(df.i$Camera_Installation %in% 'Work'),]
    
    # plot the detections
    if (MCMCdata$marks){
      df.i$Species.ind <- as.factor(df.i$Species.ind)
      p1 <- ggplot(df.i) + geom_point(aes(x = Date.Time, y = 0, col = Species.ind)) + 
        ggtitle(meta$Site_ID[i]) + xlim(meta$camera_on[i], meta$camera_off[i]) +
        theme(legend.position = "top")
    } else {
      p1 <- ggplot(df.i) + geom_point(aes(x = Date.Time, y = 0)) + 
        ggtitle(meta$Site_ID[i]) + xlim(meta$camera_on[i], meta$camera_off[i]) 
    }
    if (!is.na(meta$camera_fail[i])){
      dt = data.frame(max = max(meta$camera_fail[i],
                                df.i$Date.Time))
      p1 <- p1 + geom_vline(data = dt, mapping = aes(xintercept = max), colour = 'red')
    }
    p1 <- p1 + theme(axis.text.y=element_blank(),
                     axis.ticks.y=element_blank()) + ylab('')
      
    # posteriors for this camera
    states <- samples[,paste0('state[', as.character(i), ', ', as.character(1:n.units[i]), ']')]
    m.states <- colMeans(1*(states == state))
    ucl.states <- apply(1*(states==state), 2, quantile, probs = 0.975)
    lcl.states <- apply(1*(states==state), 2, quantile, probs = 0.025)
    
    df.i <- data.frame(times = 1:n.units[i], m.states = m.states,
                       ucl.states = ucl.states, lcl.states = lcl.states)
    
    p2 <- ggplot(data = df.i) + geom_line(aes(x = times, y = m.states)) + 
      geom_line(aes(x=times, y=ucl.states), linetype = 'dashed') + 
      geom_line(aes(x=times, y=lcl.states), linetype = 'dashed') + 
      ylim(0,1) + 
      theme(axis.text.y=element_blank(),
            axis.ticks.y=element_blank()) + ylab('')
      
    p <- grid.arrange(p1,p2, ncol=1)

    if (save.plots){
      ggsave(p, filename = paste0(meta$Site_ID[i], '.png'), ...)
    }
  }
}

# View the effort in each state of the camera
plot_use <- function(results, meta, MCMCdata, cameras = NULL, save.plots = F, ...){
  
  # cameras is named cameras
  # results needs to be a list
  
  meta <- meta[order(meta$Site_N),]
  
  if (!is.null(cameras)){
    # what cameras do we want to see (by number)
    cameras <- cameras[order(match(cameras, meta$Site_ID))] # put in order of meta
    to_plot <- meta$Site_N[meta$Site_ID %in% cameras]
  } else {
    cameras <- meta$Site_ID
    to_plot <- c(1:nrow(meta))
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
  
  n.units <- MCMCdata$n.units

  U1 <- samples[,paste0('U1[', as.character(to_plot), ']')]
  U2 <- samples[,paste0('U2[', as.character(to_plot), ']')]
  U3 <- samples[,paste0('U3[', as.character(to_plot), ']')]
  
  m.U1 <- colMeans(U1)
  m.U2 <- colMeans(U2)
  m.U3 <- colMeans(U3)
  
  q1.U1 <- apply(U1, 2, quantile, probs = 0.025)
  q1.U2 <- apply(U2, 2, quantile, probs = 0.025)
  q1.U3 <- apply(U3, 2, quantile, probs = 0.025)
  q2.U1 <- apply(U1, 2, quantile, probs = 0.975)
  q2.U2 <- apply(U2, 2, quantile, probs = 0.975)
  q2.U3 <- apply(U3, 2, quantile, probs = 0.975)
  
  df <- data.frame(matrix(nrow=3*length(to_plot), ncol=5))
  colnames(df) <- c('State', 'Mean', 'UCL', 'LCL', 'Camera')
  
  df$Camera <- rep(cameras, 3)
  df$State <- rep(1:3, each = length(to_plot))
  df$Mean <- c(m.U1, m.U2, m.U3)
  df$LCL <- c(q1.U1, q1.U2, q1.U3)
  df$UCL <- c(q2.U1, q2.U2, q2.U3)
  
  df$Camera <- as.factor(df$Camera)
  df$State <- as.factor(df$State)
  
  p <- ggplot(df, aes(x = Camera, y = Mean, color = State)) + 
    geom_point() + geom_errorbar(aes(ymin=LCL, ymax=UCL, color = State)) + 
    ylab('Time in state') +
    theme(axis.text.x = element_text(angle=90, vjust=.5, hjust=1))
  
  plot(p)
  
  if (save.plots){
    ggsave(p, filename = paste0(meta$Site_ID[i], '.png'), ...)
  }
  
} 


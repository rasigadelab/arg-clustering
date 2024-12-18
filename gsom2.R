###########################################
# GENERALIZED SOM
#'
#' With DNA diffusion
#'
#' Multi-species variant
#' Concatenate State and Pars vectors
#' Use fact that names in vector can be repeated, then separate vectors for various pops

### Original functions JP ###
LVcompet <- function(Time, State, Pars) { # State & Pars are named vectors with population i suffix _i
  # Number of species (pop)
  nSpecies <- Pars["nSpecies"]
  # Total population (all species)
  Ntotal <- sum(State[c(
    paste0("x00_", 1:nSpecies),
    paste0("x01_", 1:nSpecies),
    paste0("x10_", 1:nSpecies),
    paste0("x11_", 1:nSpecies))])
  
  stepState <- State # to store the output of the simulation
  stepState <- sapply(stepState, function(x) 0)
  
  for(i in 1:nSpecies) {
    # Argument preparation
    State_pop_i <- State[names(stepState) %like% paste0("_",i)] # variable of species i
    names(State_pop_i) <-  gsub(paste0("_",i), "", names(State_pop_i)) # remove pop suffix for populationStep()
    Pars_pop_i <- Pars[names(Pars) %like% paste0("_",i)] # parameters of species i
    names(Pars_pop_i) <-  gsub(paste0("_",i), "", names(Pars_pop_i)) # remove pop suffix for populationStep()
    # Simulation step
    stepState[names(stepState) %like% paste0("_",i)] <- populationStep(
      Time,
      State_pop_i,
      c(Pars_pop_i, Ntotal = Ntotal)
    )[[1]]
  }
  
  return(list(
    stepState
  ))
} # Use parameter & variable names

populationStep <- function(Time, State, Pars) { # population step for one species
  with(as.list(c(State, Pars)), { # with(data, expression)
    
    # Moments 1 and 2 (about zero) and parameters
    if(x11 > 0) {
      Ef1 <- NEf1 / x11
      Ef2 <- NEf2 / x11
    } else { # uniform arc distance
      Ef1 <- 1/4
      Ef2 <- 1/12
    }
    
    # Find arc distance distribution parameters
    ab <- as.double(absolve(Ef1, Ef2 - Ef1^2))
    # print(ab)
    alpha <- ab[1]
    beta  <- ab[2]
    
    # Probability of comobilization
    pab <- comobg(gamma, alpha, beta)
    
    # Conditional moments: double-gene event
    Md <- as.double(muvarg(gamma, alpha, beta))
    Ed1 <- Md[1]
    Ed2 <- Md[2] + Ed1^2
    
    # # Conditional moments: single-gene deletion
    # Es1 <- (Ef1 - pab * Ed1) / (1 - pab)
    # Es2 <- (Ef2 - pab * Ed2) / (1 - pab)
    
    # Conditional moments: single-gene insertion (uniform arc distance)
    Eu1 <- 1/4
    Eu2 <- 1/12
    
    # Intermediate values
    eta <- pab / Es # /Es is here to compensate Es in rRelease term
    eta1 <- 1 - eta
    
    # CORRECTION 21-07-30, P(!A,B)
    # Conditional moments: single-gene deletion
    Es1 <- (Es * Ef1 - pab * Ed1) / (Es * eta1)
    Es2 <- (Es * Ef2 - pab * Ed2) / (Es * eta1)
    
    # Introduce population dynamics
    N <- x00 + x01 + x10 + x11
    rRelease <- Es * (rMob + rDup)
    rLoss <- rUpt * N + rClear
    
    dDna01 <- rRelease * (x01 + eta1 * x11) - rLoss * Dna01
    dDna10 <- rRelease * (x10 + eta1 * x11) - rLoss * Dna10
    dDna11 <- rRelease * (eta * x11)        - rLoss * Dna11
    
    # Gains
    tau00_01 <- rUpt * x00 * Dna01
    tau00_10 <- rUpt * x00 * Dna10
    tau00_11 <- rUpt * x00 * Dna11
    tau01_11 <- rUpt * x01 * Dna10
    tau10_11 <- rUpt * x10 * Dna01
    
    # Losses
    l <- Es * rMob
    tau01_00 <- l * x01
    tau10_00 <- l * x10
    tau11_00 <- l * eta  * x11
    tau11_01 <- l * eta1 * x11
    tau11_10 <- l * eta1 * x11
    
    # Variable death rates
    #dR <- deathRateFunc(Time)
    #dR <- deathRateFunc_tris(Time,db=0.2,di=0.1,dp=100,trt=3,lp=200) # ok
    dR <- deathRateFunc_tris(Time,db=deathBasal,di=deathIncrease,dp=deathPeriod,trt=treatmentType,lp=lifePeriod) # ok
    
    damp <- 1 - Ntotal / K # logistic growth rate
    r00 <- (b00 * damp - dR[1])
    r01 <- (b01 * damp - dR[2])
    r10 <- (b10 * damp - dR[3])
    r11 <- (b11 * damp - dR[4])
    
    # Populations
    dx00 <- r00 * x00 + (tau01_00 + tau10_00 + tau11_00) - (tau00_01 + tau00_10 + tau00_11)
    dx01 <- r01 * x01 + (tau00_01 + tau11_01) - (tau01_11 + tau01_00)
    dx10 <- r10 * x10 + (tau00_10 + tau11_10) - (tau10_11 + tau10_00)
    dx11 <- r11 * x11 + (tau00_11 + tau01_11 + tau10_11) - (tau11_00 + tau11_01 + tau11_10)
    
    # Component derivatives (intermediate forms have underscore suffix)
    dNf_ <- r11 * x11
    dNd_ <- tau00_11 - tau11_00
    dNs_ <- -(tau11_01 + tau11_10)
    dNu_ <- +(tau01_11 + tau10_11)
    
    dNEf1 <- dNf_ * Ef1 + dNd_ * Ed1 + dNs_ * Es1 + dNu_ * Eu1
    dNEf2 <- dNf_ * Ef2 + dNd_ * Ed2 + dNs_ * Es2 + dNu_ * Eu2
    
    return(list(
      c(dx00, dx01, dx10, dx11, dNEf1, dNEf2, dDna01, dDna10, dDna11),
      # a = a, b = b,
      eta = eta#,
      # Ef1 = Ef1, Ed1 = Ed1, Es1 = Es1#,
      # x11gain = (tau00_11 + tau01_11 + tau10_11),
      # x11loss = (tau11_00 + tau11_01 + tau11_10),
      # x11conta = tau01_11 + tau10_11,
      # r01 = r01,
      # dx01 = dx01
    )) # end return
  }) # end with()
}

### Iteration method ###
LVcompet_dt <- function(Time, State, Pars, dt=0.001) { # State & Pars are named vectors with population i suffix _i
  # Number of species (pop)
  nSpecies <- Pars["nSpecies"]
  # Total population (all species)
  Ntotal <- sum(State[c(
    paste0("x00_", 1:nSpecies),
    paste0("x01_", 1:nSpecies),
    paste0("x10_", 1:nSpecies),
    paste0("x11_", 1:nSpecies))])
  
  stepState <- State # to store the output of the simulation
  stepState <- sapply(stepState, function(x) 0)
  
  for(i in 1:nSpecies) {
    # Argument preparation
    State_pop_i <- State[names(stepState) %like% paste0("_",i)] # variable of species i
    names(State_pop_i) <-  gsub(paste0("_",i), "", names(State_pop_i)) # remove pop suffix for populationStep()
    Pars_pop_i <- Pars[names(Pars) %like% paste0("_",i)] # parameters of species i
    names(Pars_pop_i) <-  gsub(paste0("_",i), "", names(Pars_pop_i)) # remove pop suffix for populationStep()
    # Simulation step
    stepState[names(stepState) %like% paste0("_",i)] <- populationStep_dt(
      Time,
      State_pop_i,
      c(Pars_pop_i, Ntotal = Ntotal),
      dt
    )[[1]]
  }
  
  return(list(
    stepState
  ))
}

populationStep_dt <- function(Time, State, Pars, dt) { # population step for one species
  with(as.list(c(State, Pars)), { # with(data, expression)
    
    # Moments 1 and 2 (about zero) and parameters
    if(x11 > 0) {
      Ef1 <- NEf1 / x11
      Ef2 <- NEf2 / x11
    } else { # uniform arc distance
      Ef1 <- 1/4
      Ef2 <- 1/12
    }
    
    # Find arc distance distribution parameters
    ab <- as.double(absolve(Ef1, Ef2 - Ef1^2))
    # print(ab)
    alpha <- ab[1]
    beta  <- ab[2]
    
    # Probability of comobilization
    pab <- comobg(gamma, alpha, beta)
    
    # Conditional moments: double-gene event
    Md <- as.double(muvarg(gamma, alpha, beta))
    Ed1 <- Md[1]
    Ed2 <- Md[2] + Ed1^2
    
    # Conditional moments: single-gene deletion
    Es1 <- (Ef1 - pab * Ed1) / (1 - pab)
    Es2 <- (Ef2 - pab * Ed2) / (1 - pab)
    
    # Conditional moments: single-gene insertion (uniform arc distance)
    Eu1 <- 1/4
    Eu2 <- 1/12
    
    # Intermediate values
    eta <- pab / Es # /Es is here to compensate Es in rRelease term
    eta1 <- 1 - eta
    
    # Introduce population dynamics
    N <- x00 + x01 + x10 + x11
    rRelease <- Es * (rMob + rDup)
    rLoss <- rUpt * N + rClear
    
    dDna01 <- (rRelease * (x01 + eta1 * x11) - rLoss * Dna01) * dt # no cell uptake ?
    dDna10 <- (rRelease * (x10 + eta1 * x11) - rLoss * Dna10) * dt
    dDna11 <- (rRelease * (eta * x11)        - rLoss * Dna11) * dt
    
    # Gains
    tau00_01 <- rUpt * x00 * Dna01 
    tau00_10 <- rUpt * x00 * Dna10 
    tau00_11 <- rUpt * x00 * Dna11 
    tau01_11 <- rUpt * x01 * Dna10 
    tau10_11 <- rUpt * x10 * Dna01 
    
    # Losses
    l <- (Es * rMob) 
    tau01_00 <- l * x01 
    tau10_00 <- l * x10 
    tau11_00 <- l * eta  * x11 
    tau11_01 <- l * eta1 * x11
    tau11_10 <- l * eta1 * x11 
    
    # Variable death rates
    dR <- deathRateFunc(Time)
    
    damp <- 1 - Ntotal / K
    r00 <- b00 * damp - dR[1] 
    r01 <- b01 * damp - dR[2] 
    r10 <- b10 * damp - dR[3] 
    r11 <- b11 * damp - dR[4] 
    
    # Populations
    dx00 <- (r00 * x00 + (tau01_00 + tau10_00 + tau11_00) - (tau00_01 + tau00_10 + tau00_11)) * dt
    dx01 <- (r01 * x01 + (tau00_01 + tau11_01) - (tau01_11 + tau01_00)) * dt
    dx10 <- (r10 * x10 + (tau00_10 + tau11_10) - (tau10_11 + tau10_00)) * dt
    dx11 <- (r11 * x11 + (tau00_11 + tau01_11 + tau10_11) - (tau11_00 + tau11_01 + tau11_10)) * dt
    
    # Component derivatives (intermediate forms have underscore suffix)
    dNf_ <- (r11 * x11) * dt
    dNd_ <- (tau00_11 - tau11_00) * dt
    dNs_ <- (-(tau11_01 + tau11_10)) * dt
    dNu_ <- (+(tau01_11 + tau10_11)) * dt
    
    dNEf1 <- dNf_ * Ef1 + dNd_ * Ed1 + dNs_ * Es1 + dNu_ * Eu1
    dNEf2 <- dNf_ * Ef2 + dNd_ * Ed2 + dNs_ * Es2 + dNu_ * Eu2
    
    x00 <- x00 + dx00
    x01 <- x01 + dx01
    x10 <- x10 + dx10
    x11 <- x11 + dx11
    NEf1 <- NEf1 + dNEf1
    NEf2 <- NEf2 + dNEf2
    Dna01 <- Dna01 + dDna01
    Dna10 <- Dna10 + dDna10
    Dna11 <- Dna11 + dDna11
    if(x00 < 0){x00 = 0}
    if(x01 < 0){x01 = 0}
    if(x10 < 0){x10 = 0}
    if(x11 < 0){x11 = 0}
    if(Dna01 < 0){Dna01 = 0}
    if(Dna10 < 0){Dna10 = 0}
    if(Dna11 < 0){Dna11 = 0}
    
    return(list(
      #c(dx00, dx01, dx10, dx11, dNEf1, dNEf2, dDna01, dDna10, dDna11),
      c(x00, x01, x10, x11, NEf1, NEf2, Dna01, Dna10, Dna11),
      # a = a, b = b,
      eta = eta#,
      # Ef1 = Ef1, Ed1 = Ed1, Es1 = Es1#,
      # x11gain = (tau00_11 + tau01_11 + tau10_11),
      # x11loss = (tau11_00 + tau11_01 + tau11_10),
      # x11conta = tau01_11 + tau10_11,
      # r01 = r01,
      # dx01 = dx01
    )) # end return
  }) # end with()
}

### 2 patch model ###
LVcompet_2patch <- function(Time, State, Pars, DEBUG = FALSE) { # State & Pars are named vectors with population i suffix _i
  # Number of species (pop)
  nSpecies <- Pars["nSpecies"]
  # Total population (all species)
  Ntotal_env1 <- sum(State[c(
    paste0("x00_env1_sp", 1:nSpecies),
    paste0("x01_env1_sp", 1:nSpecies),
    paste0("x10_env1_sp", 1:nSpecies),
    paste0("x11_env1_sp", 1:nSpecies))])
  Ntotal_env2 <- sum(State[c(
    paste0("x00_env2_sp", 1:nSpecies),
    paste0("x01_env2_sp", 1:nSpecies),
    paste0("x10_env2_sp", 1:nSpecies),
    paste0("x11_env2_sp", 1:nSpecies))])
  
  stepState <- State # to store the output of the simulation
  stepState <- sapply(stepState, function(x) 0)
  
  for(i in 1:nSpecies) {
    # Argument preparation
    State_pop_i <- State[names(stepState) %like% paste0("_sp",i)] # variable of species i
    names(State_pop_i) <-  gsub(paste0("_sp",i), "", names(State_pop_i)) # remove pop suffix for populationStep()
    Pars_pop_i <- Pars[names(Pars) %like% paste0("_sp",i)] # parameters of species i
    names(Pars_pop_i) <-  gsub(paste0("_sp",i), "", names(Pars_pop_i)) # remove pop suffix for populationStep()
    # Simulation step
    stepState[names(stepState) %like% paste0("_sp",i)] <- populationStep_2patch(
      Time,
      State_pop_i,
      c(Pars_pop_i, Ntotal_env1 = Ntotal_env1, Ntotal_env2 = Ntotal_env2),
      DEBUG = DEBUG
    )[[1]]
  }
  
  return(list(
    stepState
  ))
}

populationStep_2patch <- function(Time, State, Pars, DEBUG = FALSE) { # population step for one species
  with(as.list(c(State, Pars)), { # with(data, expression)
    
    # Extract parameters and variables (!for test only!)
    # for(i in 1:length(State)){
    #   assign(gsub("_sp1","",names(State[i])), State[i])
    # }
    # for(i in 1:length(Pars)){
    #   assign(gsub("_sp1","",names(Pars[i])), Pars[i])
    # }
    
    for(envId in 1:2){ # iterate through environments
      # Get state of environment envId
      x00 <- get(paste0("x00_env", envId))
      x01 <- get(paste0("x01_env", envId))
      x10 <- get(paste0("x10_env", envId))
      x11 <- get(paste0("x11_env", envId))
      NEf1 <- get(paste0("NEf1_env", envId))
      NEf2 <- get(paste0("NEf2_env", envId))
      Dna01 <- get(paste0("Dna01_env", envId))
      Dna10 <- get(paste0("Dna10_env", envId))
      Dna11 <- get(paste0("Dna11_env", envId))
      Ntotal <- get(paste0("Ntotal_env", envId))
      
      if(DEBUG){
        #if(x00<0 | x01<0 | x10<0 | x11<0 | Dna01<0 | Dna10<0 | Dna11<0){ # prevent negative populations
        if(x00<(-0.1) | x01<(-0.1) | x10<(-0.1) | x11<(-0.1) | Dna01<(-0.01) | Dna10<(-0.01) | Dna11<(-0.01)){ # prevent negative population (with small error marge)
          print(paste0(
            " time:",Time," env:",envId," x00:",x00," x01:",x01," x10:",x10," x11:",x11," Dna01:",Dna01," Dna10:",Dna10," Dna11:",Dna11
          ))
          stop("ERROR: cell or Dna < 0")
        }
      }
      
      # Moments 1 and 2 (about zero) and parameters
      if(x11 > 1E-6) { # >0 may lead to simulation failure when x11 tends to 0
        Ef1 <- NEf1 / x11
        Ef2 <- NEf2 / x11
      } else { # uniform arc distance
        Ef1 <- 1/4
        Ef2 <- 1/12
      }
      
      # Find arc distance distribution parameters
      ab <- as.double(absolve(Ef1, Ef2 - Ef1^2))
      # print(ab)
      alpha <- ab[1]
      beta  <- ab[2]
      
      # Probability of comobilization
      pab <- comobg(gamma, alpha, beta)
      
      # Conditional moments: double-gene event
      Md <- as.double(muvarg(gamma, alpha, beta))
      Ed1 <- Md[1]
      Ed2 <- Md[2] + Ed1^2
      
      # Conditional moments: single-gene deletion (OBSOLETE)
      # Es1 <- (Ef1 - pab * Ed1) / (1 - pab)
      # Es2 <- (Ef2 - pab * Ed2) / (1 - pab)
      
      # Conditional moments: single-gene insertion (uniform arc distance)
      Eu1 <- 1/4
      Eu2 <- 1/12
      
      # Intermediate values
      eta <- pab / Es # /Es is here to compensate Es in rRelease term
      eta1 <- 1 - eta
      
      # TENTATIVE CORRECTION 21-07-30, P(!A,B)
      # Conditional moments: single-gene deletion
      Es1 <- (Es * Ef1 - pab * Ed1) / (Es * eta1)
      Es2 <- (Es * Ef2 - pab * Ed2) / (Es * eta1)
      
      if(DEBUG){
        if(Ef1<0 | Ef1>0.2501 | Ef2<0 | Ef2>0.8334 | # >1/4 & >1/12 + error marge allowed
           Ed1<0 | Ed1>0.2501 | Ed2<0 | Ed2>0.8334 |
           #Es1<0 | Es1>1/4 | Es2<0 | Es2>1/12 |
           #Eu1<0 | Eu1>1/4 | Eu2<0 | Eu2>1/12 |
           pab<0 | pab>1 | eta<0 | eta>1){
          print(paste0(
            " time:",Time," env:",envId," x00:",x00," x01:",x01," x10:",x10," x11:",x11," Dna01:",Dna01," Dna10:",Dna10," Dna11:",Dna11
          ))
          print(paste0(
            " alpha:",alpha," beta:",beta," pab:",pab," Ef1:",Ef1," Ef2:",Ef2,
            " Ed1:",Ed1," Ed2:",Ed2," Es1:",Es1, " Es2:",Es2," Eu1:",Eu1," Eu2:",Eu2," eta:",eta
          ))
          stop("ERROR: probability or moments out of bounds")
        }
      }
      
      # Introduce population dynamics
      N <- x00 + x01 + x10 + x11
      rRelease <- Es * (rMob + rDup)
      rLoss <- rUpt * N + rClear
      
      dDna01 <- rRelease * (x01 + eta1 * x11) - rLoss * Dna01
      dDna10 <- rRelease * (x10 + eta1 * x11) - rLoss * Dna10
      dDna11 <- rRelease * (eta * x11)        - rLoss * Dna11
      
      # Gains
      tau00_01 <- rUpt * x00 * Dna01
      tau00_10 <- rUpt * x00 * Dna10
      tau00_11 <- rUpt * x00 * Dna11
      tau01_11 <- rUpt * x01 * Dna10
      tau10_11 <- rUpt * x10 * Dna01
      
      # Losses
      l <- Es * rMob
      tau01_00 <- l * x01
      tau10_00 <- l * x10
      tau11_00 <- l * eta  * x11
      tau11_01 <- l * eta1 * x11
      tau11_10 <- l * eta1 * x11
      
      # Variable death rates
      if(envId==1){ # constant antibiotic
        #dR <- deathRateFunc(Time)
        dR <- c(deathBasal+deathIncrease,deathBasal+deathIncrease,deathBasal+deathIncrease,deathBasal)
        damp <- 1 - Ntotal / K1 # logistic growth rate
      } else { # no antibiotic
        dR <- c(deathBasal,deathBasal,deathBasal,deathBasal)
        damp <- 1 - Ntotal / K2 # logistic growth rate
      }
      
      r00 <- (b00 * damp - dR[1])
      r01 <- (b01 * damp - dR[2])
      r10 <- (b10 * damp - dR[3])
      r11 <- (b11 * damp - dR[4])
      
      # Populations
      dx00 <- r00 * x00 + (tau01_00 + tau10_00 + tau11_00) - (tau00_01 + tau00_10 + tau00_11)
      dx01 <- r01 * x01 + (tau00_01 + tau11_01) - (tau01_11 + tau01_00)
      dx10 <- r10 * x10 + (tau00_10 + tau11_10) - (tau10_11 + tau10_00)
      dx11 <- r11 * x11 + (tau00_11 + tau01_11 + tau10_11) - (tau11_00 + tau11_01 + tau11_10)
      
      # Component derivatives (intermediate forms have underscore suffix)
      dNf_ <- r11 * x11 # unchanged arc distance
      dNd_ <- tau00_11 - tau11_00 # double gene events (gain & loss)
      dNs_ <- -(tau11_01 + tau11_10) # double gene to single gene
      dNu_ <- +(tau01_11 + tau10_11) # single gene to double gene
      
      dNEf1 <- dNf_ * Ef1 + dNd_ * Ed1 + dNs_ * Es1 + dNu_ * Eu1
      dNEf2 <- dNf_ * Ef2 + dNd_ * Ed2 + dNs_ * Es2 + dNu_ * Eu2
      
      if(envId==1){
        # Cell migration
        dx00 <- dx00 + G * (x00_env2 - x00) 
        dx01 <- dx01 + G * (x01_env2 - x01)
        dx10 <- dx10 + G * (x10_env2 - x10)
        dx11 <- dx11 + G * (x11_env2 - x11)
        # Arc distance modifications
        if(x11_env2 > 0) {
          Ef1_env2 <- NEf1_env2 / x11_env2 # Moments 1 and 2 (about zero) and parameters
          Ef2_env2 <- NEf2_env2 / x11_env2
        } else { # uniform arc distance
          Ef1_env2 <- 1/4
          Ef2_env2 <- 1/12
        }
        # dNin_ <- G * x11_env2 # x11 cells comming from env2
        # dNout_ <- G * x11 # x11 cells leaving env1
        # dNEf1 <- dNEf1 + dNin_ * Ef1_env2 - dNout_ * Ef1 # arc distance variation due to x11 cells from environment2
        # dNEf2 <- dNEf2 + dNin_ * Ef2_env2 - dNout_ * Ef2
        dNEf1 <- dNEf1 + G * (x11_env2 * Ef1_env2 - x11 * Ef1)
        dNEf2 <- dNEf2 + G * (x11_env2 * Ef2_env2 - x11 * Ef2)
      } # end env1 specs
      if(envId==2){
        # Cell migration
        dx00 <- dx00 + G * (x00_env1 - x00) 
        dx01 <- dx01 + G * (x01_env1 - x01)
        dx10 <- dx10 + G * (x10_env1 - x10)
        dx11 <- dx11 + G * (x11_env1 - x11)
        # Arc distance modifications
        if(x11_env1 > 0) {
          Ef1_env1 <- NEf1_env1 / x11_env1 # Moments 1 and 2 (about zero) and parameters
          Ef2_env1 <- NEf2_env1 / x11_env1
        } else { # uniform arc distance
          Ef1_env1 <- 1/4
          Ef2_env1 <- 1/12
        }
        # dNin_ <- G * x11_env1 # x11 cells comming from env1
        # dNout_ <- G * x11 # cells leaving env2
        # dNEf1 <- dNEf1 + dNin_ * Ef1_env1 - dNout_ * Ef1 # arc distance variation due to x11 cells from environment2
        # dNEf2 <- dNEf2 + dNin_ * Ef2_env1 - dNout_ * Ef2
        dNEf1 <- dNEf1 + G * (x11_env1 * Ef1_env1 - x11 * Ef1)
        dNEf2 <- dNEf2 + G * (x11_env1 * Ef2_env1 - x11 * Ef2)
      } # end env2 specs
      
      state_env <- c(dx00, dx01, dx10, dx11, dNEf1, dNEf2, dDna01, dDna10, dDna11)
      names(state_env) <- paste0(c("dx00", "dx01", "dx10", "dx11", "dNEf1", "dNEf2", "dDna01", "dDna10", "dDna11"), "_env", envId)
      assign(paste0("state_env",envId), state_env)
    } # end environment loop
    to_return <- c(state_env1, state_env2)
    
    return(list(
      to_return,
      #c(dx00, dx01, dx10, dx11, dNEf1, dNEf2, dDna01, dDna10, dDna11),
      # a = a, b = b,
      eta = eta#,
      # Ef1 = Ef1, Ed1 = Ed1, Es1 = Es1#,
      # x11gain = (tau00_11 + tau01_11 + tau10_11),
      # x11loss = (tau11_00 + tau11_01 + tau11_10),
      # x11conta = tau01_11 + tau10_11,
      # r01 = r01,
      # dx01 = dx01
    )) # end return
  }) # end with()
}

### 2 patch model - iteration method ###
LVcompet_2patch_dt <- function(Time, State, Pars, dt=0.001) { # State & Pars are named vectors with population i suffix _i
  # Number of species (pop)
  nSpecies <- Pars["nSpecies"]
  # Total population (all species)
  Ntotal_env1 <- sum(State[c(
    paste0("x00_env1_sp", 1:nSpecies),
    paste0("x01_env1_sp", 1:nSpecies),
    paste0("x10_env1_sp", 1:nSpecies),
    paste0("x11_env1_sp", 1:nSpecies))])
  Ntotal_env2 <- sum(State[c(
    paste0("x00_env2_sp", 1:nSpecies),
    paste0("x01_env2_sp", 1:nSpecies),
    paste0("x10_env2_sp", 1:nSpecies),
    paste0("x11_env2_sp", 1:nSpecies))])
  
  stepState <- State # to store the output of the simulation
  stepState <- sapply(stepState, function(x) 0)
  
  for(i in 1:nSpecies) {
    # Argument preparation
    State_pop_i <- State[names(stepState) %like% paste0("_sp",i)] # variable of species i
    names(State_pop_i) <-  gsub(paste0("_sp",i), "", names(State_pop_i)) # remove pop suffix for populationStep()
    Pars_pop_i <- Pars[names(Pars) %like% paste0("_sp",i)] # parameters of species i
    names(Pars_pop_i) <-  gsub(paste0("_sp",i), "", names(Pars_pop_i)) # remove pop suffix for populationStep()
    # Simulation step
    stepState[names(stepState) %like% paste0("_sp",i)] <- populationStep_2patch_dt(
      Time,
      State_pop_i,
      c(Pars_pop_i, Ntotal_env1 = Ntotal_env1, Ntotal_env2 = Ntotal_env2),
      dt
    )[[1]]
  }
  
  return(list(
    stepState
  ))
}

populationStep_2patch_dt <- function(Time, State, Pars, dt) { # population step for one species
  with(as.list(c(State, Pars)), { # with(data, expression)
    
    # Extract parameters and variables (!for test only!)
    # for(i in 1:length(State)){
    #   assign(gsub("_sp1","",names(State[i])), State[i])
    # }
    # for(i in 1:length(Pars)){
    #   assign(gsub("_sp1","",names(Pars[i])), Pars[i])
    # }
    
    for(envId in 1:2){ # iterate through environments
      # Get state of environment envId
      x00 <- get(paste0("x00_env", envId))
      x01 <- get(paste0("x01_env", envId))
      x10 <- get(paste0("x10_env", envId))
      x11 <- get(paste0("x11_env", envId))
      NEf1 <- get(paste0("NEf1_env", envId))
      NEf2 <- get(paste0("NEf2_env", envId))
      Dna01 <- get(paste0("Dna01_env", envId))
      Dna10 <- get(paste0("Dna10_env", envId))
      Dna11 <- get(paste0("Dna11_env", envId))
      Ntotal <- get(paste0("Ntotal_env", envId))
      
      # Moments 1 and 2 (about zero) and parameters
      if(x11 > 0) {
        Ef1 <- NEf1 / x11
        Ef2 <- NEf2 / x11
      } else { # uniform arc distance
        Ef1 <- 1/4
        Ef2 <- 1/12
      }
      
      # Find arc distance distribution parameters
      ab <- as.double(absolve(Ef1, Ef2 - Ef1^2))
      # print(ab)
      alpha <- ab[1]
      beta  <- ab[2]
      
      # Probability of comobilization
      pab <- comobg(gamma, alpha, beta)
      
      # Conditional moments: double-gene event
      Md <- as.double(muvarg(gamma, alpha, beta))
      Ed1 <- Md[1]
      Ed2 <- Md[2] + Ed1^2
      
      # Conditional moments: single-gene deletion
      Es1 <- (Ef1 - pab * Ed1) / (1 - pab)
      Es2 <- (Ef2 - pab * Ed2) / (1 - pab)
      
      # Conditional moments: single-gene insertion (uniform arc distance)
      Eu1 <- 1/4
      Eu2 <- 1/12
      
      # Intermediate values
      eta <- pab / Es # /Es is here to compensate Es in rRelease term
      eta1 <- 1 - eta
      
      # Introduce population dynamics
      N <- x00 + x01 + x10 + x11
      rRelease <- Es * (rMob + rDup)
      rLoss <- rUpt * N + rClear
      
      dDna01 <- rRelease * (x01 + eta1 * x11) - rLoss * Dna01
      dDna10 <- rRelease * (x10 + eta1 * x11) - rLoss * Dna10
      dDna11 <- rRelease * (eta * x11)        - rLoss * Dna11
      
      # Gains
      tau00_01 <- rUpt * x00 * Dna01
      tau00_10 <- rUpt * x00 * Dna10
      tau00_11 <- rUpt * x00 * Dna11
      tau01_11 <- rUpt * x01 * Dna10
      tau10_11 <- rUpt * x10 * Dna01
      
      # Losses
      l <- Es * rMob
      tau01_00 <- l * x01
      tau10_00 <- l * x10
      tau11_00 <- l * eta  * x11
      tau11_01 <- l * eta1 * x11
      tau11_10 <- l * eta1 * x11
      
      # Variable death rates
      if(envId==1){
        dR <- deathRateFunc(Time)
      } else {
        dR <- c(deathBasal,deathBasal,deathBasal,deathBasal)
      }
      
      damp <- 1 - Ntotal / K # logistic growth rate
      r00 <- (b00 * damp - dR[1])
      r01 <- (b01 * damp - dR[2])
      r10 <- (b10 * damp - dR[3])
      r11 <- (b11 * damp - dR[4])
      
      # Populations
      dx00 <- r00 * x00 + (tau01_00 + tau10_00 + tau11_00) - (tau00_01 + tau00_10 + tau00_11)
      dx01 <- r01 * x01 + (tau00_01 + tau11_01) - (tau01_11 + tau01_00)
      dx10 <- r10 * x10 + (tau00_10 + tau11_10) - (tau10_11 + tau10_00)
      dx11 <- r11 * x11 + (tau00_11 + tau01_11 + tau10_11) - (tau11_00 + tau11_01 + tau11_10)
      
      # Component derivatives (intermediate forms have underscore suffix)
      dNf_ <- r11 * x11
      dNd_ <- tau00_11 - tau11_00
      dNs_ <- -(tau11_01 + tau11_10)
      dNu_ <- +(tau01_11 + tau10_11)
      
      dNEf1 <- dNf_ * Ef1 + dNd_ * Ed1 + dNs_ * Es1 + dNu_ * Eu1
      dNEf2 <- dNf_ * Ef2 + dNd_ * Ed2 + dNs_ * Es2 + dNu_ * Eu2
      
      if(envId==1){
        # Migration
        dx00 <- dx00 - G * x00 + G * x00_env2
        dx01 <- dx01 - G * x01 + G * x01_env2
        dx10 <- dx10 - G * x10 + G * x10_env2
        dx11 <- dx11 - G * x11 + G * x11_env2
        # Arc distance modifications
        if(x11_env2 > 0) {
          Ef1_env2 <- NEf1_env2 / x11_env2 # Moments 1 and 2 (about zero) and parameters
          Ef2_env2 <- NEf2_env2 / x11_env2
        } else { # uniform arc distance
          Ef1_env2 <- 1/4
          Ef2_env2 <- 1/12
        }
        dNin_ <- G * x11_env2 # x11 cells comming from env2
        dNout_ <- G * x11 # x11 cells leaving env1
        dNEf1 <- dNEf1 + dNin_ * Ef1_env2 - dNout_ * Ef1_env2 # arc distance variation due to x11 cells from environment2
        dNEf2 <- dNEf2 + dNin_ * Ef2_env2 - dNout_ * Ef2_env2
      } # end env1 specs
      if(envId==2){
        # Migration
        dx00 <- dx00 - G * x00 + G * x00_env1
        dx01 <- dx01 - G * x01 + G * x01_env1
        dx10 <- dx10 - G * x10 + G * x10_env1
        dx11 <- dx11 - G * x11 + G * x11_env1
        # Arc distance modifications
        if(x11_env1 > 0) {
          Ef1_env1 <- NEf1_env1 / x11_env1 # Moments 1 and 2 (about zero) and parameters
          Ef2_env1 <- NEf2_env1 / x11_env1
        } else { # uniform arc distance
          Ef1_env1 <- 1/4
          Ef2_env1 <- 1/12
        }
        dNin_ <- G * x11_env1 # x11 cells comming from env1
        dNout_ <- G * x11 # cells leaving env2
        dNEf1 <- dNEf1 + dNin_ * Ef1_env1 - dNout_ * Ef1_env2 # arc distance variation due to x11 cells from environment2
        dNEf2 <- dNEf2 + dNin_ * Ef2_env1 - dNout_ * Ef2_env2
      } # end env2 specs
      
      # Update state
      x00 <- x00 + dx00 * dt
      x01 <- x01 + dx01 * dt
      x10 <- x10 + dx10 * dt
      x11 <- x11 + dx11 * dt
      NEf1 <- NEf1 + dNEf1 * dt
      NEf2 <- NEf2 + dNEf2 * dt
      Dna01 <- Dna01 + dDna01 * dt
      Dna10 <- Dna10 + dDna10 * dt
      Dna11 <- Dna11 + dDna11 * dt
      # Prevent negative values
      if(x00 < 0){x00 = 0}
      if(x01 < 0){x01 = 0}
      if(x10 < 0){x10 = 0}
      if(x11 < 0){x11 = 0}
      if(Dna01 < 0){Dna01 = 0}
      if(Dna10 < 0){Dna10 = 0}
      if(Dna11 < 0){Dna11 = 0}
      
      state_env <- c(x00, x01, x10, x11, NEf1, NEf2, Dna01, Dna10, Dna11)
      names(state_env) <- paste0(c("x00", "x01", "x10", "x11", "NEf1", "NEf2", "Dna01", "Dna10", "Dna11"), "_env", envId)
      assign(paste0("state_env",envId), state_env)
    } # end environment loop
    to_return <- c(state_env1, state_env2)
    
    return(list(
      to_return,
      #c(dx00, dx01, dx10, dx11, dNEf1, dNEf2, dDna01, dDna10, dDna11),
      # a = a, b = b,
      eta = eta#,
      # Ef1 = Ef1, Ed1 = Ed1, Es1 = Es1#,
      # x11gain = (tau00_11 + tau01_11 + tau10_11),
      # x11loss = (tau11_00 + tau11_01 + tau11_10),
      # x11conta = tau01_11 + tau10_11,
      # r01 = r01,
      # dx01 = dx01
    )) # end return
  }) # end with()
}

### Wild Type input / x00_niche ###
LVcompet_wtinput <- function(Time, State, Pars, DEBUG=FALSE) { # State & Pars are named vectors with population i suffix _i
  # Number of species (pop)
  nSpecies <- Pars["nSpecies"]
  # Total population (all species)
  Ntotal <- sum(State[c(
    paste0("x00_", 1:nSpecies),
    paste0("x01_", 1:nSpecies),
    paste0("x10_", 1:nSpecies),
    paste0("x11_", 1:nSpecies))])
  
  stepState <- State # to store the output of the simulation
  stepState <- sapply(stepState, function(x) 0)
  
  for(i in 1:nSpecies) {
    # Argument preparation
    State_pop_i <- State[names(stepState) %like% paste0("_",i)] # variable of species i
    names(State_pop_i) <-  gsub(paste0("_",i), "", names(State_pop_i)) # remove pop suffix for populationStep()
    Pars_pop_i <- Pars[names(Pars) %like% paste0("_",i)] # parameters of species i
    names(Pars_pop_i) <-  gsub(paste0("_",i), "", names(Pars_pop_i)) # remove pop suffix for populationStep()
    # Simulation step
    stepState[names(stepState) %like% paste0("_",i)] <- populationStep_wtinput(
      Time,
      State_pop_i,
      c(Pars_pop_i, Ntotal = Ntotal),
      DEBUG=DEBUG
    )[[1]]
  }
  
  return(list(
    stepState
  ))
}

populationStep_wtinput <- function(Time, State, Pars, DEBUG=FALSE) { # population step for one species
  with(as.list(c(State, Pars)), { # with(data, expression)
    
    if(DEBUG){
      #if(x00<0 | x01<0 | x10<0 | x11<0 | Dna01<0 | Dna10<0 | Dna11<0){ # prevent negative populations
      if(x00<(-0.1) | x01<(-0.1) | x10<(-0.1) | x11<(-0.1) | Dna01<(-0.01) | Dna10<(-0.01) | Dna11<(-0.01)){ # prevnt negative population (with small error marge)
        print(paste0(
          "x00:",x00," x01:",x01," x10:",x10," x11:",x11,
          " Dna01:",Dna01," Dna10:",Dna10," Dna11:",Dna11,
          " Ntotal:",Ntotal
        ))
        stop("ERROR: cell or Dna < 0")
      }
    }
    
    # Moments 1 and 2 (about zero) and parameters
    if(x11 > 1E-6) { # using x11>0 may lead to simulation failure when x11 tends to 0
      # ! arc distance results are not valid for small x11 population
      Ef1 <- NEf1 / x11
      Ef2 <- NEf2 / x11
    } else { # uniform arc distance
      Ef1 <- 1/4
      Ef2 <- 1/12
    }
    
    # Find arc distance distribution parameters
    ab <- as.double(absolve(Ef1, Ef2 - Ef1^2))
    alpha <- ab[1]
    beta  <- ab[2]
    
    # Probability of comobilization
    # comobg: unconditional probability of cotransfer
    pab <- comobg(gamma, alpha, beta)
    
    # Conditional moments: double-gene event
    # i.e. arc distance mean and variance conditional to a HGT event involving both genes
    Md <- as.double(muvarg(gamma, alpha, beta))
    Ed1 <- Md[1]
    Ed2 <- Md[2] + Ed1^2
    
    # Conditional moments: single-gene deletion, P(!(AB)) (!OBSOLETE!)
    #Es1 <- (Ef1 - pab * Ed1) / (1 - pab)
    #Es2 <- (Ef2 - pab * Ed2) / (1 - pab)
    
    # Conditional moments: single-gene insertion (uniform arc distance)
    Eu1 <- 1/4
    Eu2 <- 1/12
    
    # Intermediate values
    # Es: probability P(A)=P(B) that an HGT event involves A and/or B
    # eta: probability P(A|B) comobilization conditional to HGT event involves A and/or B
    # pab: probability comobilization when HGT event
    eta <- pab / Es # conditional probability comobilization
    eta1 <- 1 - eta # conditional probability single mobilization
    
    # CORRECTION 21-07-30, P(!A,B)
    # Conditional moments: single-gene deletion
    Es1 <- (Es * Ef1 - pab * Ed1) / (Es * eta1)
    Es2 <- (Es * Ef2 - pab * Ed2) / (Es * eta1)
    
    if(DEBUG){
      if(Ef1<0 | Ef1>0.2501 | Ef2<0 | Ef2>0.8334 | # >1/4 & >1/12 + error marge allowed
         Ed1<0 | Ed1>0.2501 | Ed2<0 | Ed2>0.8334 |
         #Es1<0 | Es1>1/4 | Es2<0 | Es2>1/12 |
         #Eu1<0 | Eu1>1/4 | Eu2<0 | Eu2>1/12 |
         pab<0 | pab>1 | eta<0 | eta>1){
        print(paste0(
          " time:",Time," x00:",x00," x01:",x01," x10:",x10," x11:",x11," Dna01:",Dna01," Dna10:",Dna10," Dna11:",Dna11
        ))
        print(paste0(
          " alpha:",alpha," beta:",beta," pab:",pab," Ef1:",Ef1," Ef2:",Ef2,
          " Ed1:",Ed1," Ed2:",Ed2," Es1:",Es1, " Es2:",Es2," Eu1:",Eu1," Eu2:",Eu2," eta:",eta
        ))
        stop("ERROR: probability or moments out of bounds")
      }
    }
    
    # Add x00 niche cells (constant x00 population)
    x00_tot <- x00 + x00_niche
    
    # Introduce population dynamics
    N <- x00_tot + x01 + x10 + x11 # total cells of this species
    rRelease <- Es * (rMob + rDup) # DNA release rate
    rLoss <- rUpt * N + rClear # DNA removal rate
    
    # eDNA dynamics
    dDna01 <- rRelease * (x01 + eta1 * x11) - rLoss * Dna01
    dDna10 <- rRelease * (x10 + eta1 * x11) - rLoss * Dna10
    dDna11 <- rRelease * (eta * x11)        - rLoss * Dna11
    
    # Gene gains
    tau00_01 <- rUpt * x00_tot * Dna01
    tau00_10 <- rUpt * x00_tot * Dna10
    tau00_11 <- rUpt * x00_tot * Dna11
    tau01_11 <- rUpt * x01 * Dna10
    tau10_11 <- rUpt * x10 * Dna01
    
    # Gene losses
    l <- Es * rMob
    tau01_00 <- l * x01
    tau10_00 <- l * x10
    tau11_00 <- l * eta  * x11 # double-gene loss
    tau11_01 <- l * eta1 * x11 # single gene loss
    tau11_10 <- l * eta1 * x11 # single gene loss
    
    # Variable death rates
    dR <- deathRateFunc_tris(Time, db=deathBasal, di=deathIncrease, dp=deathPeriod, trt=treatmentType, lp=lifePeriod)
    
    if(DEBUG){
      if(is.na(dR[1]) | is.na(dR[2]) | is.na(dR[3]) | is.na(dR[4])){
        stop("ERROR: death rate function output NA")
      }
    }
    
    # Bacterial growth
    damp <- 1 - Ntotal / K # logistic growth rate
    r00 <- (b00 * damp - dR[1]) # original formula
    r01 <- (b01 * damp - dR[2]) # original formula
    r10 <- (b10 * damp - dR[3]) # original formula
    r11 <- (b11 * damp - dR[4]) # original formula
    
    # Cell dynamics
    dx00 <- r00 * x00 + (tau01_00 + tau10_00 + tau11_00) - (tau00_01 + tau00_10 + tau00_11) + WTin * (K - x00)
    dx01 <- r01 * x01 + (tau00_01 + tau11_01) - (tau01_11 + tau01_00) - WTin * x01
    dx10 <- r10 * x10 + (tau00_10 + tau11_10) - (tau10_11 + tau10_00) - WTin * x10
    dx11 <- r11 * x11 + (tau00_11 + tau01_11 + tau10_11) - (tau11_00 + tau11_01 + tau11_10) - WTin * x11
    
    # Component derivatives (intermediate forms have underscore suffix)
    dNf_ <- (r11 - WTin) * x11 # unchanged arc distance
    dNd_ <- tau00_11 - tau11_00 # double gene gain and loss
    dNs_ <- -(tau11_01 + tau11_10) # double gene to single gene
    dNu_ <- +(tau01_11 + tau10_11) # single gene to double gene
    
    # Arc distance dynamics
    dNEf1 <- dNf_ * Ef1 + dNd_ * Ed1 + dNs_ * Es1 + dNu_ * Eu1
    dNEf2 <- dNf_ * Ef2 + dNd_ * Ed2 + dNs_ * Es2 + dNu_ * Eu2
    
    return(list(
      c(dx00, dx01, dx10, dx11, dNEf1, dNEf2, dDna01, dDna10, dDna11),
      # a = a, b = b,
      eta = eta#,
      # Ef1 = Ef1, Ed1 = Ed1, Es1 = Es1#,
      # x11gain = (tau00_11 + tau01_11 + tau10_11),
      # x11loss = (tau11_00 + tau11_01 + tau11_10),
      # x11conta = tau01_11 + tau10_11,
      # r01 = r01,
      # dx01 = dx01
    )) # end return
  }) # end with()
}

### Wild Type input / x00_niche + additional outputs for debugging ###
LVcompet_wtinput_bis <- function(Time, State, Pars, DEBUG=FALSE) { # State & Pars are named vectors with population i suffix _i
  # Number of species (pop)
  nSpecies <- Pars["nSpecies"]
  # Total population (all species)
  Ntotal <- sum(State[c(
    paste0("x00_", 1:nSpecies),
    paste0("x01_", 1:nSpecies),
    paste0("x10_", 1:nSpecies),
    paste0("x11_", 1:nSpecies))])
  
  stepState <- State # to store the output of the simulation
  stepState <- sapply(stepState, function(x) 0)
  
  for(i in 1:nSpecies) {
    # Argument preparation
    State_pop_i <- State[names(stepState) %like% paste0("_",i)] # variable of species i
    names(State_pop_i) <-  gsub(paste0("_",i), "", names(State_pop_i)) # remove pop suffix for populationStep()
    Pars_pop_i <- Pars[names(Pars) %like% paste0("_",i)] # parameters of species i
    names(Pars_pop_i) <-  gsub(paste0("_",i), "", names(Pars_pop_i)) # remove pop suffix for populationStep()
    # Simulation step
    pop_step <- populationStep_wtinput_bis(
      Time,
      State_pop_i,
      c(Pars_pop_i, Ntotal = Ntotal),
      DEBUG=DEBUG
    )
    stepState[names(stepState) %like% paste0("_",i)] <- pop_step[[1]] # state variables
    if(i == 1){
      intermediate_variables_1 <- pop_step[[2]]
    }
  }
  
  return(list(
    stepState,
    intermediate_variables_1
  ))
}

populationStep_wtinput_bis <- function(Time, State, Pars, DEBUG=FALSE) { # population step for one species
  with(as.list(c(State, Pars)), { # with(data, expression)
    
    if(DEBUG){
      #if(x00<0 | x01<0 | x10<0 | x11<0 | Dna01<0 | Dna10<0 | Dna11<0){
      if(x00<(-0.1) | x01<(-0.1) | x10<(-0.1) | x11<(-0.1) | Dna01<(-0.01) | Dna10<(-0.01) | Dna11<(-0.01)){
        print(paste0(
          "x00:",x00," x01:",x01," x10:",x10," x11:",x11,
          " Dna01:",Dna01," Dna10:",Dna10," Dna11:",Dna11,
          " Ntotal:",Ntotal
        ))
        stop("ERROR: cell or Dna < 0")
      }
    }
    
    # Moments 1 and 2 (about zero) and parameters
    if(x11 > 0) {
      Ef1 <- NEf1 / x11
      Ef2 <- NEf2 / x11
    } else { # uniform arc distance
      Ef1 <- 1/4
      Ef2 <- 1/12
    }
    
    # Find arc distance distribution parameters
    ab <- as.double(absolve(Ef1, Ef2 - Ef1^2))
    alpha <- ab[1]
    beta  <- ab[2]
    
    # Probability of comobilization
    # comobg: unconditional probability of cotransfer
    pab <- comobg(gamma, alpha, beta)
    
    # Conditional moments: double-gene event
    # i.e. arc distance mean and variance conditional to a HGT event involving both genes
    Md <- as.double(muvarg(gamma, alpha, beta))
    Ed1 <- Md[1]
    Ed2 <- Md[2] + Ed1^2
    
    # Conditional moments: single-gene deletion
    #Es1 <- (Ef1 - pab * Ed1) / (1 - pab)
    #Es2 <- (Ef2 - pab * Ed2) / (1 - pab)
    
    # Conditional moments: single-gene insertion (uniform arc distance)
    Eu1 <- 1/4
    Eu2 <- 1/12
    
    # Intermediate values
    # 1 / Es: probability P(A)=P(B) probability that an HGT event involves A and/or B
    # pab: probability P(A|B), probability comobilization conditional to HGT event involves A and/or B
    eta <- pab / Es # conditional probability comobilization
    eta1 <- 1 - eta # conditional probability single mobilization
    
    # TENTATIVE CORRECTION 21-07-30, P(!A,B)
    # Conditional moments: single-gene deletion
    Es1 <- (Es * Ef1 - pab * Ed1) / (Es * eta1)
    Es2 <- (Es * Ef2 - pab * Ed2) / (Es * eta1)
    
    if(DEBUG){
      if(Ef1<0 | Ef1>0.2501 | Ef2<0 | Ef2>0.0834 | # 1/4 & 1/12 + error tolerance
         Ed1<0 | Ed1>0.2501 | Ed2<0 | Ed2>0.0834 |
         #Es1<0 | Es1>1/4 | Es2<0 | Es2>1/12 |
         #Eu1<0 | Eu1>1/4 | Eu2<0 | Eu2>1/12 |
         pab<0 | pab>1 | eta<0 | eta>1){
        print(paste0(
          " time:",Time," x00:",x00," x01:",x01," x10:",x10," x11:",x11," Dna01:",Dna01," Dna10:",Dna10," Dna11:",Dna11
        ))
        print(paste0(
          " alpha:",alpha," beta:",beta," pab:",pab," Ef1:",Ef1," Ef2:",Ef2,
          " Ed1:",Ed1," Ed2:",Ed2," Es1:",Es1, " Es2:",Es2," Eu1:",Eu1," Eu2:",Eu2," eta:",eta
        ))
        stop("ERROR: probability or moments out of bounds")
      }
    }
    
    # Add x00 niche cells (constant x00 population)
    x00_tot <- x00 + x00_niche
    
    # Introduce population dynamics
    N <- x00_tot + x01 + x10 + x11 # total cells of these species
    rRelease <- Es * (rMob + rDup) # DNA release rate
    rLoss <- rUpt * N + rClear # original formula
    
    # dDna01 <- rRelease * (x01 + eta1 * x11) - rLoss * Dna01 # original formula
    # dDna10 <- rRelease * (x10 + eta1 * x11) - rLoss * Dna10 # original formula
    dDna01 <- rRelease * (x01 + eta1 * x11) - rLoss * Dna01 # (50/50 gene 01 or 10)
    dDna10 <- rRelease * (x10 + eta1 * x11) - rLoss * Dna10 # (50/50 gene 01 or 10)
    dDna11 <- rRelease * (eta * x11)        - rLoss * Dna11 # original formula
    
    # Gains
    # tau00_01 <- rUpt * x00_tot * Dna01
    # tau00_10 <- rUpt * x00_tot * Dna10
    # tau00_11 <- rUpt * x00_tot * Dna11
    tau00_01 <- rUpt * x00_tot * Dna01
    tau00_10 <- rUpt * x00_tot * Dna10
    tau00_11 <- rUpt * x00_tot * Dna11
    tau01_11 <- rUpt * x01 * Dna10
    tau10_11 <- rUpt * x10 * Dna01
    
    # Losses
    l <- Es * rMob
    tau01_00 <- l * x01
    tau10_00 <- l * x10
    tau11_00 <- l * eta  * x11 # double-gene loss
    tau11_01 <- l * eta1 * x11 # single gene loss
    tau11_10 <- l * eta1 * x11 # single gene loss
    
    # Variable death rates
    #dR <- deathRateFunc(Time)#, deathBasal, deathIncrease, deathPeriod)
    #dR <- deathRateFunc_bis(Time, deathBasal, deathIncrease, deathPeriod, treatment=treatmentType, lifePeriod)
    dR <- deathRateFunc_tris(Time, db=deathBasal, di=deathIncrease, dp=deathPeriod, trt=treatmentType, lp=lifePeriod)
    
    #print(paste0("time: ",Time," dR_x00: ",dR[1]))
    #print(paste0("db:",deathBasal, " di:",deathIncrease, " dp:",deathPeriod, " trt:",treatmentType, " lp:",lifePeriod))
    
    damp <- 1 - Ntotal / K # logistic growth rate
    r00 <- (b00 * damp - dR[1]) # original formula
    r01 <- (b01 * damp - dR[2]) # original formula
    r10 <- (b10 * damp - dR[3]) # original formula
    r11 <- (b11 * damp - dR[4]) # original formula
    
    # Populations
    # dx00 <- r00 * x00 + (tau01_00 + tau10_00 + tau11_00) - (tau00_01 + tau00_10 + tau00_11)
    # dx01 <- r01 * x01 + (tau00_01 + tau11_01) - (tau01_11 + tau01_00)
    # dx10 <- r10 * x10 + (tau00_10 + tau11_10) - (tau10_11 + tau10_00)
    # dx11 <- r11 * x11 + (tau00_11 + tau01_11 + tau10_11) - (tau11_00 + tau11_01 + tau11_10)
    dx00 <- r00 * x00 + (tau01_00 + tau10_00 + tau11_00) - (tau00_01 + tau00_10 + tau00_11) + WTin * (K - x00)
    dx01 <- r01 * x01 + (tau00_01 + tau11_01) - (tau01_11 + tau01_00) - WTin * x01
    dx10 <- r10 * x10 + (tau00_10 + tau11_10) - (tau10_11 + tau10_00) - WTin * x10
    dx11 <- r11 * x11 + (tau00_11 + tau01_11 + tau10_11) - (tau11_00 + tau11_01 + tau11_10) - WTin * x11
    
    # Component derivatives (intermediate forms have underscore suffix)
    #dNf_ <- r11 * x11 # unchanged arc distance
    #dNf_ <- (r11 - WTin) * x11 # unchanged arc distance
    dNf_ <- (r11 - WTin) * x11 #- (tau11_00 + tau11_01 + tau11_10) # unchanged arc distance
    dNd_ <- tau00_11 - tau11_00 # double gene gain and loss
    dNs_ <- -(tau11_01 + tau11_10) # double gene to single gene
    dNu_ <- +(tau01_11 + tau10_11) # single gene to double gene
    
    dNEf1 <- dNf_ * Ef1 + dNd_ * Ed1 + dNs_ * Es1 + dNu_ * Eu1
    dNEf2 <- dNf_ * Ef2 + dNd_ * Ed2 + dNs_ * Es2 + dNu_ * Eu2
    
    return(list(
      c(dx00, dx01, dx10, dx11, dNEf1, dNEf2, dDna01, dDna10, dDna11),
      c(pab, eta, dNf_, Ef1, Ef2, dNd_, Ed1, Ed2, dNs_, Es1, Es2, dNu_, Eu1, Eu2)
      # a = a, b = b,
      #eta = eta#,
      # Ef1 = Ef1, Ed1 = Ed1, Es1 = Es1#,
      # x11gain = (tau00_11 + tau01_11 + tau10_11),
      # x11loss = (tau11_00 + tau11_01 + tau11_10),
      # x11conta = tau01_11 + tau10_11,
      # r01 = r01,
      # dx01 = dx01
    )) # end return
  }) # end with()
}

### Shared eDNA pool ###
LVcompet_sharedDNA <- function(Time, State, Pars, DEBUG=FALSE) { # State & Pars are named vectors with population i suffix _i
  # Number of species (pop)
  nSpecies <- Pars["nSpecies"]
  # Total population (all species)
  Ntotal <- sum(State[c(
    paste0("x00_", 1:nSpecies),
    paste0("x01_", 1:nSpecies),
    paste0("x10_", 1:nSpecies),
    paste0("x11_", 1:nSpecies))])
  
  # Total eDNA pool (all species)
  Dna11_total <- sum(State[paste0("Dna11_", 1:nSpecies)])
  Dna01_total <- sum(State[paste0("Dna01_", 1:nSpecies)])
  Dna10_total <- sum(State[paste0("Dna10_", 1:nSpecies)])
  
  if(DEBUG){
    if(Dna01_total<(-0.01) | Dna10_total<(-0.01) | Dna11_total<(-0.01)){ # prevnt negative population (with small error marge)
      print(paste0(
        " Dna01_total:",Dna01_total," Dna10_total:",Dna10_total," Dna11_total:",Dna11_total
      ))
      stop("ERROR: cell or Dna < 0")
    }
  }
  
  stepState <- State # to store the output of the simulation
  stepState <- sapply(stepState, function(x) 0)
  
  for(i in 1:nSpecies) {
    # Argument preparation
    State_pop_i <- State[names(stepState) %like% paste0("_",i)] # variable of species i
    names(State_pop_i) <-  gsub(paste0("_",i), "", names(State_pop_i)) # remove pop suffix for populationStep()
    Pars_pop_i <- Pars[names(Pars) %like% paste0("_",i)] # parameters of species i
    names(Pars_pop_i) <-  gsub(paste0("_",i), "", names(Pars_pop_i)) # remove pop suffix for populationStep()
    # Simulation step
    stepState[names(stepState) %like% paste0("_",i)] <- populationStep_sharedDNA(
      Time,
      State_pop_i,
      c(Pars_pop_i, Ntotal=Ntotal, Dna11_total=Dna11_total, Dna01_total=Dna10_total, Dna10_total=Dna10_total),
      DEBUG=DEBUG
    )[[1]]
  }
  # stepState contains the derivatives of the current step
  # eDNA derivatives account for DNA release and consumption by each species to the total eDNA
  # DNA total derivatrive is the sum ofthe contribution of each species
  dDna11_total <- sum(stepState[names(stepState) %like% "Dna11_"]) - Pars["rClear_1"] * Dna11_total
  dDna01_total <- sum(stepState[names(stepState) %like% "Dna01_"]) - Pars["rClear_1"] * Dna01_total
  dDna10_total <- sum(stepState[names(stepState) %like% "Dna10_"]) - Pars["rClear_1"] * Dna10_total
  stepState[names(stepState) %like% "Dna11_"] <- dDna11_total / nSpecies
  stepState[names(stepState) %like% "Dna01_"] <- dDna01_total / nSpecies
  stepState[names(stepState) %like% "Dna10_"] <- dDna10_total / nSpecies
  
  return(list(
    stepState
  ))
}

populationStep_sharedDNA <- function(Time, State, Pars, DEBUG=FALSE) { # population step for one species
  with(as.list(c(State, Pars)), { # with(data, expression)
    
    if(DEBUG){
      #if(x00<0 | x01<0 | x10<0 | x11<0 | Dna01<0 | Dna10<0 | Dna11<0){ # prevent negative populations
      if(x00<(-0.1) | x01<(-0.1) | x10<(-0.1) | x11<(-0.1)){#} | Dna01<(-0.01) | Dna10<(-0.01) | Dna11<(-0.01)){ # prevnt negative population (with small error marge)
        print(paste0(
          "x00:",x00," x01:",x01," x10:",x10," x11:",x11,
          " Dna01:",Dna01," Dna10:",Dna10," Dna11:",Dna11,
          " Ntotal:",Ntotal
        ))
        stop("ERROR: cell or Dna < 0")
      }
    }
    
    # Moments 1 and 2 (about zero) and parameters
    if(x11 > 1E-6) { # using x11>0 may lead to simulation failure when x11 tends to 0
      # ! arc distance results are not valid for small x11 population
      Ef1 <- NEf1 / x11
      Ef2 <- NEf2 / x11
    } else { # uniform arc distance
      Ef1 <- 1/4
      Ef2 <- 1/12
    }
    
    # Find arc distance distribution parameters
    ab <- as.double(absolve(Ef1, Ef2 - Ef1^2))
    alpha <- ab[1]
    beta  <- ab[2]
    
    # Probability of comobilization
    # comobg: unconditional probability of cotransfer
    pab <- comobg(gamma, alpha, beta)
    
    # Conditional moments: double-gene event
    # i.e. arc distance mean and variance conditional to a HGT event involving both genes
    Md <- as.double(muvarg(gamma, alpha, beta))
    Ed1 <- Md[1]
    Ed2 <- Md[2] + Ed1^2
    
    # Conditional moments: single-gene deletion, P(!(AB)) (!OBSOLETE!)
    #Es1 <- (Ef1 - pab * Ed1) / (1 - pab)
    #Es2 <- (Ef2 - pab * Ed2) / (1 - pab)
    
    # Conditional moments: single-gene insertion (uniform arc distance)
    Eu1 <- 1/4
    Eu2 <- 1/12
    
    # Intermediate values
    # Es: probability P(A)=P(B) that an HGT event involves A and/or B
    # eta: probability P(A|B) comobilization conditional to HGT event involves A and/or B
    # pab: probability comobilization when HGT event
    eta <- pab / Es # conditional probability comobilization
    eta1 <- 1 - eta # conditional probability single mobilization
    
    # CORRECTION 21-07-30, P(!A,B)
    # Conditional moments: single-gene deletion
    Es1 <- (Es * Ef1 - pab * Ed1) / (Es * eta1)
    Es2 <- (Es * Ef2 - pab * Ed2) / (Es * eta1)
    
    if(DEBUG){
      if(Ef1<0 | Ef1>0.2501 | Ef2<0 | Ef2>0.8334 | # >1/4 & >1/12 + error marge allowed
         Ed1<0 | Ed1>0.2501 | Ed2<0 | Ed2>0.8334 |
         #Es1<0 | Es1>1/4 | Es2<0 | Es2>1/12 |
         #Eu1<0 | Eu1>1/4 | Eu2<0 | Eu2>1/12 |
         pab<0 | pab>1 | eta<0 | eta>1){
        print(paste0(
          " time:",Time," x00:",x00," x01:",x01," x10:",x10," x11:",x11," Dna01:",Dna01," Dna10:",Dna10," Dna11:",Dna11
        ))
        print(paste0(
          " alpha:",alpha," beta:",beta," pab:",pab," Ef1:",Ef1," Ef2:",Ef2,
          " Ed1:",Ed1," Ed2:",Ed2," Es1:",Es1, " Es2:",Es2," Eu1:",Eu1," Eu2:",Eu2," eta:",eta
        ))
        stop("ERROR: probability or moments out of bounds")
      }
    }
    
    # Add x00 niche cells (constant x00 population)
    x00_tot <- x00 + x00_niche
    
    # Introduce population dynamics
    N <- x00_tot + x01 + x10 + x11 # total cells of this species
    rRelease <- Es * (rMob + rDup) # DNA release rate
    #rLoss <- rUpt * N + rClear # DNA removal rate
    
    # eDNA dynamics
    #dDna01 <- rRelease * (x01 + eta1 * x11) - rLoss * Dna01
    #dDna10 <- rRelease * (x10 + eta1 * x11) - rLoss * Dna10
    #dDna11 <- rRelease * (eta * x11)        - rLoss * Dna11
    Dna_total <- Dna01_total + Dna10_total + Dna11_total
    rUpt_reg <- rUpt * Dna_total / (Dna_total + Dna_cst) # limit DNA uptake when eDNA is low
    if(Dna_total > 0){
      dDna01 <- rRelease * (x01 + eta1 * x11) - rUpt_reg * N * Dna01_total / Dna_total
      dDna10 <- rRelease * (x10 + eta1 * x11) - rUpt_reg * N * Dna10_total / Dna_total
      dDna11 <- rRelease * (eta * x11)        - rUpt_reg * N * Dna11_total / Dna_total
      
      tau00_01 <- rUpt_reg * x00_tot * Dna01_total / Dna_total
      tau00_10 <- rUpt_reg * x00_tot * Dna10_total / Dna_total
      tau00_11 <- rUpt_reg * x00_tot * Dna11_total / Dna_total
      tau01_11 <- rUpt_reg * x01 * Dna10_total / Dna_total
      tau10_11 <- rUpt_reg * x10 * Dna01_total / Dna_total
    } else {
      dDna01 <- rRelease * (x01 + eta1 * x11)
      dDna10 <- rRelease * (x10 + eta1 * x11)
      dDna11 <- rRelease * (eta * x11)
      
      tau00_01 <- 0
      tau00_10 <- 0
      tau00_11 <- 0
      tau01_11 <- 0
      tau10_11 <- 0
    }
    # Gene gains (old)
    # tau00_01 <- rUpt * x00_tot * Dna01_total
    # tau00_10 <- rUpt * x00_tot * Dna10_total
    # tau00_11 <- rUpt * x00_tot * Dna11_total
    # tau01_11 <- rUpt * x01 * Dna10_total
    # tau10_11 <- rUpt * x10 * Dna01_total
    
    # Gene losses
    l <- Es * rMob
    tau01_00 <- l * x01
    tau10_00 <- l * x10
    tau11_00 <- l * eta  * x11 # double-gene loss
    tau11_01 <- l * eta1 * x11 # single gene loss
    tau11_10 <- l * eta1 * x11 # single gene loss
    
    # Variable death rates
    dR <- deathRateFunc_tris(Time, db=deathBasal, di=deathIncrease, dp=deathPeriod, trt=treatmentType, lp=lifePeriod)
    
    if(DEBUG){
      if(is.na(dR[1]) | is.na(dR[2]) | is.na(dR[3]) | is.na(dR[4])){
        stop("ERROR: death rate function output NA")
      }
    }
    
    # Bacterial growth
    damp <- 1 - Ntotal / K # logistic growth rate
    r00 <- (b00 * damp - dR[1]) # original formula
    r01 <- (b01 * damp - dR[2]) # original formula
    r10 <- (b10 * damp - dR[3]) # original formula
    r11 <- (b11 * damp - dR[4]) # original formula
    
    # Cell dynamics
    dx00 <- r00 * x00 + (tau01_00 + tau10_00 + tau11_00) - (tau00_01 + tau00_10 + tau00_11) + WTin * (K - x00)
    dx01 <- r01 * x01 + (tau00_01 + tau11_01) - (tau01_11 + tau01_00) - WTin * x01
    dx10 <- r10 * x10 + (tau00_10 + tau11_10) - (tau10_11 + tau10_00) - WTin * x10
    dx11 <- r11 * x11 + (tau00_11 + tau01_11 + tau10_11) - (tau11_00 + tau11_01 + tau11_10) - WTin * x11
    
    # Component derivatives (intermediate forms have underscore suffix)
    dNf_ <- (r11 - WTin) * x11 # unchanged arc distance
    dNd_ <- tau00_11 - tau11_00 # double gene gain and loss
    dNs_ <- -(tau11_01 + tau11_10) # double gene to single gene
    dNu_ <- +(tau01_11 + tau10_11) # single gene to double gene
    
    # Arc distance dynamics
    dNEf1 <- dNf_ * Ef1 + dNd_ * Ed1 + dNs_ * Es1 + dNu_ * Eu1
    dNEf2 <- dNf_ * Ef2 + dNd_ * Ed2 + dNs_ * Es2 + dNu_ * Eu2
    
    return(list(
      c(dx00, dx01, dx10, dx11, dNEf1, dNEf2, dDna01, dDna10, dDna11),
      # a = a, b = b,
      eta = eta#,
      # Ef1 = Ef1, Ed1 = Ed1, Es1 = Es1#,
      # x11gain = (tau00_11 + tau01_11 + tau10_11),
      # x11loss = (tau11_00 + tau11_01 + tau11_10),
      # x11conta = tau01_11 + tau10_11,
      # r01 = r01,
      # dx01 = dx01
    )) # end return
  }) # end with()
}

### rUpt = DNA_tot / (DNA_tot + cst) ###
LVcompet_uptReg <- function(Time, State, Pars, DEBUG=FALSE) { # State & Pars are named vectors with population i suffix _i
  # Number of species (pop)
  nSpecies <- Pars["nSpecies"]
  # Total population (all species)
  Ntotal <- sum(State[c(
    paste0("x00_", 1:nSpecies),
    paste0("x01_", 1:nSpecies),
    paste0("x10_", 1:nSpecies),
    paste0("x11_", 1:nSpecies))])
  
  stepState <- State # to store the output of the simulation
  stepState <- sapply(stepState, function(x) 0)
  
  for(i in 1:nSpecies) {
    # Argument preparation
    State_pop_i <- State[names(stepState) %like% paste0("_",i)] # variable of species i
    names(State_pop_i) <-  gsub(paste0("_",i), "", names(State_pop_i)) # remove pop suffix for populationStep()
    Pars_pop_i <- Pars[names(Pars) %like% paste0("_",i)] # parameters of species i
    names(Pars_pop_i) <-  gsub(paste0("_",i), "", names(Pars_pop_i)) # remove pop suffix for populationStep()
    # Simulation step
    stepState[names(stepState) %like% paste0("_",i)] <- populationStep_uptReg(
      Time,
      State_pop_i,
      c(Pars_pop_i, Ntotal = Ntotal),
      DEBUG=DEBUG
    )[[1]]
  }
  
  return(list(
    stepState
  ))
}

populationStep_uptReg <- function(Time, State, Pars, DEBUG=FALSE) { # population step for one species
  with(as.list(c(State, Pars)), { # with(data, expression)
    
    if(DEBUG){
      #if(x00<0 | x01<0 | x10<0 | x11<0 | Dna01<0 | Dna10<0 | Dna11<0){ # prevent negative populations
      if(x00<(-0.1) | x01<(-0.1) | x10<(-0.1) | x11<(-0.1) | Dna01<(-0.01) | Dna10<(-0.01) | Dna11<(-0.01)){ # prevnt negative population (with small error marge)
        print(paste0(
          "x00:",x00," x01:",x01," x10:",x10," x11:",x11,
          " Dna01:",Dna01," Dna10:",Dna10," Dna11:",Dna11,
          " Ntotal:",Ntotal
        ))
        stop("ERROR: cell or Dna < 0")
      }
    }
    
    # Moments 1 and 2 (about zero) and parameters
    if(x11 > 1E-6) { # using x11>0 may lead to simulation failure when x11 tends to 0
      # ! arc distance results are not valid for small x11 population
      Ef1 <- NEf1 / x11
      Ef2 <- NEf2 / x11
    } else { # uniform arc distance
      Ef1 <- 1/4
      Ef2 <- 1/12
    }
    
    # Find arc distance distribution parameters
    ab <- as.double(absolve(Ef1, Ef2 - Ef1^2))
    alpha <- ab[1]
    beta  <- ab[2]
    
    # Probability of comobilization
    # comobg: unconditional probability of cotransfer
    pab <- comobg(gamma, alpha, beta)
    
    # Conditional moments: double-gene event
    # i.e. arc distance mean and variance conditional to a HGT event involving both genes
    Md <- as.double(muvarg(gamma, alpha, beta))
    Ed1 <- Md[1]
    Ed2 <- Md[2] + Ed1^2
    
    # Conditional moments: single-gene deletion, P(!(AB)) (!OBSOLETE!)
    #Es1 <- (Ef1 - pab * Ed1) / (1 - pab)
    #Es2 <- (Ef2 - pab * Ed2) / (1 - pab)
    
    # Conditional moments: single-gene insertion (uniform arc distance)
    Eu1 <- 1/4
    Eu2 <- 1/12
    
    # Intermediate values
    # Es: probability P(A)=P(B) that an HGT event involves A and/or B
    # eta: probability P(A|B) comobilization conditional to HGT event involves A and/or B
    # pab: probability comobilization when HGT event
    eta <- pab / Es # conditional probability comobilization
    eta1 <- 1 - eta # conditional probability single mobilization
    
    # CORRECTION 21-07-30, P(!A,B)
    # Conditional moments: single-gene deletion
    Es1 <- (Es * Ef1 - pab * Ed1) / (Es * eta1)
    Es2 <- (Es * Ef2 - pab * Ed2) / (Es * eta1)
    
    if(DEBUG){
      if(Ef1<0 | Ef1>0.2501 | Ef2<0 | Ef2>0.8334 | # >1/4 & >1/12 + error marge allowed
         Ed1<0 | Ed1>0.2501 | Ed2<0 | Ed2>0.8334 |
         #Es1<0 | Es1>1/4 | Es2<0 | Es2>1/12 |
         #Eu1<0 | Eu1>1/4 | Eu2<0 | Eu2>1/12 |
         pab<0 | pab>1 | eta<0 | eta>1){
        print(paste0(
          " time:",Time," x00:",x00," x01:",x01," x10:",x10," x11:",x11," Dna01:",Dna01," Dna10:",Dna10," Dna11:",Dna11
        ))
        print(paste0(
          " alpha:",alpha," beta:",beta," pab:",pab," Ef1:",Ef1," Ef2:",Ef2,
          " Ed1:",Ed1," Ed2:",Ed2," Es1:",Es1, " Es2:",Es2," Eu1:",Eu1," Eu2:",Eu2," eta:",eta
        ))
        stop("ERROR: probability or moments out of bounds")
      }
    }
    
    # Add x00 niche cells (constant x00 population)
    x00_tot <- x00 + x00_niche
    
    # Introduce population dynamics
    N <- x00_tot + x01 + x10 + x11 # total cells of this species
    rRelease <- Es * (rMob + rDup) # DNA release rate
    #rLoss <- rUpt * N + rClear # DNA removal rate
    Dna_total <- Dna01 + Dna10 + Dna11
    rUpt_reg <- rUpt * Dna_total / (Dna_total + Dna_cst) # limit DNA uptake when DNA low
    
    # eDNA dynamics
    if(Dna_total > 0){
      dDna01 <- rRelease * (x01 + eta1 * x11) - rUpt_reg * N * Dna01 / Dna_total - rClear * Dna01
      dDna10 <- rRelease * (x10 + eta1 * x11) - rUpt_reg * N * Dna10 / Dna_total - rClear * Dna10
      dDna11 <- rRelease * (eta * x11)        - rUpt_reg * N * Dna11 / Dna_total - rClear * Dna11
      # Gene gains
      tau00_01 <- rUpt_reg * x00_tot * Dna01 / Dna_total
      tau00_10 <- rUpt_reg * x00_tot * Dna10 / Dna_total
      tau00_11 <- rUpt_reg * x00_tot * Dna11 / Dna_total
      tau01_11 <- rUpt_reg * x01 * Dna10 / Dna_total
      tau10_11 <- rUpt_reg * x10 * Dna01 / Dna_total
    } else {
      dDna01 <- rRelease * (x01 + eta1 * x11)
      dDna10 <- rRelease * (x10 + eta1 * x11)
      dDna11 <- rRelease * (eta * x11)
      tau00_01 <- 0
      tau00_10 <- 0
      tau00_11 <- 0
      tau01_11 <- 0
      tau10_11 <- 0
    }
    
    # Gene losses
    l <- Es * rMob
    tau01_00 <- l * x01
    tau10_00 <- l * x10
    tau11_00 <- l * eta  * x11 # double-gene loss
    tau11_01 <- l * eta1 * x11 # single gene loss
    tau11_10 <- l * eta1 * x11 # single gene loss
    
    # Variable death rates
    dR <- deathRateFunc_tris(Time, db=deathBasal, di=deathIncrease, dp=deathPeriod, trt=treatmentType, lp=lifePeriod)
    
    if(DEBUG){
      if(is.na(dR[1]) | is.na(dR[2]) | is.na(dR[3]) | is.na(dR[4])){
        stop("ERROR: death rate function output NA")
      }
    }
    
    # Bacterial growth
    damp <- 1 - Ntotal / K # logistic growth rate
    r00 <- (b00 * damp - dR[1]) # original formula
    r01 <- (b01 * damp - dR[2]) # original formula
    r10 <- (b10 * damp - dR[3]) # original formula
    r11 <- (b11 * damp - dR[4]) # original formula
    
    # Cell dynamics
    dx00 <- r00 * x00 + (tau01_00 + tau10_00 + tau11_00) - (tau00_01 + tau00_10 + tau00_11) + WTin * (K - x00)
    dx01 <- r01 * x01 + (tau00_01 + tau11_01) - (tau01_11 + tau01_00) - WTin * x01
    dx10 <- r10 * x10 + (tau00_10 + tau11_10) - (tau10_11 + tau10_00) - WTin * x10
    dx11 <- r11 * x11 + (tau00_11 + tau01_11 + tau10_11) - (tau11_00 + tau11_01 + tau11_10) - WTin * x11
    
    # Component derivatives (intermediate forms have underscore suffix)
    dNf_ <- (r11 - WTin) * x11 # unchanged arc distance
    dNd_ <- tau00_11 - tau11_00 # double gene gain and loss
    dNs_ <- -(tau11_01 + tau11_10) # double gene to single gene
    dNu_ <- +(tau01_11 + tau10_11) # single gene to double gene
    
    # Arc distance dynamics
    dNEf1 <- dNf_ * Ef1 + dNd_ * Ed1 + dNs_ * Es1 + dNu_ * Eu1
    dNEf2 <- dNf_ * Ef2 + dNd_ * Ed2 + dNs_ * Es2 + dNu_ * Eu2
    
    return(list(
      c(dx00, dx01, dx10, dx11, dNEf1, dNEf2, dDna01, dDna10, dDna11),
      # a = a, b = b,
      eta = eta#,
      # Ef1 = Ef1, Ed1 = Ed1, Es1 = Es1#,
      # x11gain = (tau00_11 + tau01_11 + tau10_11),
      # x11loss = (tau11_00 + tau11_01 + tau11_10),
      # x11conta = tau01_11 + tau10_11,
      # r01 = r01,
      # dx01 = dx01
    )) # end return
  }) # end with()
}

### 2 patch model - rUpt reg ###
LVcompet_2patchReg <- function(Time, State, Pars, DEBUG = FALSE) { # State & Pars are named vectors with population i suffix _i
  # Number of species (pop)
  nSpecies <- Pars["nSpecies"]
  # Total population (all species)
  Ntotal_env1 <- sum(State[c(
    paste0("x00_env1_sp", 1:nSpecies),
    paste0("x01_env1_sp", 1:nSpecies),
    paste0("x10_env1_sp", 1:nSpecies),
    paste0("x11_env1_sp", 1:nSpecies))])
  Ntotal_env2 <- sum(State[c(
    paste0("x00_env2_sp", 1:nSpecies),
    paste0("x01_env2_sp", 1:nSpecies),
    paste0("x10_env2_sp", 1:nSpecies),
    paste0("x11_env2_sp", 1:nSpecies))])
  
  stepState <- State # to store the output of the simulation
  stepState <- sapply(stepState, function(x) 0)
  
  for(i in 1:nSpecies) {
    # Argument preparation
    State_pop_i <- State[names(stepState) %like% paste0("_sp",i)] # variable of species i
    names(State_pop_i) <-  gsub(paste0("_sp",i), "", names(State_pop_i)) # remove pop suffix for populationStep()
    Pars_pop_i <- Pars[names(Pars) %like% paste0("_sp",i)] # parameters of species i
    names(Pars_pop_i) <-  gsub(paste0("_sp",i), "", names(Pars_pop_i)) # remove pop suffix for populationStep()
    # Simulation step
    stepState[names(stepState) %like% paste0("_sp",i)] <- populationStep_2patchReg(
      Time,
      State_pop_i,
      c(Pars_pop_i, Ntotal_env1 = Ntotal_env1, Ntotal_env2 = Ntotal_env2),
      DEBUG = DEBUG
    )[[1]]
  }
  
  return(list(
    stepState
  ))
}

populationStep_2patchReg <- function(Time, State, Pars, DEBUG = FALSE) { # population step for one species
  with(as.list(c(State, Pars)), { # with(data, expression)
    
    # Extract parameters and variables (!for test only!)
    # for(i in 1:length(State)){
    #   assign(gsub("_sp1","",names(State[i])), State[i])
    # }
    # for(i in 1:length(Pars)){
    #   assign(gsub("_sp1","",names(Pars[i])), Pars[i])
    # }
    
    for(envId in 1:2){ # iterate through environments
      # Get state of environment envId
      x00 <- get(paste0("x00_env", envId))
      x01 <- get(paste0("x01_env", envId))
      x10 <- get(paste0("x10_env", envId))
      x11 <- get(paste0("x11_env", envId))
      NEf1 <- get(paste0("NEf1_env", envId))
      NEf2 <- get(paste0("NEf2_env", envId))
      Dna01 <- get(paste0("Dna01_env", envId))
      Dna10 <- get(paste0("Dna10_env", envId))
      Dna11 <- get(paste0("Dna11_env", envId))
      Ntotal <- get(paste0("Ntotal_env", envId))
      
      if(DEBUG){
        #if(x00<0 | x01<0 | x10<0 | x11<0 | Dna01<0 | Dna10<0 | Dna11<0){ # prevent negative populations
        if(x00<(-0.1) | x01<(-0.1) | x10<(-0.1) | x11<(-0.1) | Dna01<(-0.01) | Dna10<(-0.01) | Dna11<(-0.01)){ # prevent negative population (with small error marge)
          print(paste0(
            " time:",Time," env:",envId," x00:",x00," x01:",x01," x10:",x10," x11:",x11," Dna01:",Dna01," Dna10:",Dna10," Dna11:",Dna11
          ))
          stop("ERROR: cell or Dna < 0")
        }
      }
      
      # Moments 1 and 2 (about zero) and parameters
      if(x11 > 1E-6) { # >0 may lead to simulation failure when x11 tends to 0
        Ef1 <- NEf1 / x11
        Ef2 <- NEf2 / x11
      } else { # uniform arc distance
        Ef1 <- 1/4
        Ef2 <- 1/12
      }
      
      # Find arc distance distribution parameters
      ab <- as.double(absolve(Ef1, Ef2 - Ef1^2))
      # print(ab)
      alpha <- ab[1]
      beta  <- ab[2]
      
      # Probability of comobilization
      pab <- comobg(gamma, alpha, beta)
      
      # Conditional moments: double-gene event
      Md <- as.double(muvarg(gamma, alpha, beta))
      Ed1 <- Md[1]
      Ed2 <- Md[2] + Ed1^2
      
      # Conditional moments: single-gene deletion (OBSOLETE)
      # Es1 <- (Ef1 - pab * Ed1) / (1 - pab)
      # Es2 <- (Ef2 - pab * Ed2) / (1 - pab)
      
      # Conditional moments: single-gene insertion (uniform arc distance)
      Eu1 <- 1/4
      Eu2 <- 1/12
      
      # Intermediate values
      eta <- pab / Es # /Es is here to compensate Es in rRelease term
      eta1 <- 1 - eta
      
      # TENTATIVE CORRECTION 21-07-30, P(!A,B)
      # Conditional moments: single-gene deletion
      Es1 <- (Es * Ef1 - pab * Ed1) / (Es * eta1)
      Es2 <- (Es * Ef2 - pab * Ed2) / (Es * eta1)
      
      if(DEBUG){
        if(Ef1<0 | Ef1>0.2501 | Ef2<0 | Ef2>0.8334 | # >1/4 & >1/12 + error marge allowed
           Ed1<0 | Ed1>0.2501 | Ed2<0 | Ed2>0.8334 |
           #Es1<0 | Es1>1/4 | Es2<0 | Es2>1/12 |
           #Eu1<0 | Eu1>1/4 | Eu2<0 | Eu2>1/12 |
           pab<0 | pab>1 | eta<0 | eta>1){
          print(paste0(
            " time:",Time," env:",envId," x00:",x00," x01:",x01," x10:",x10," x11:",x11," Dna01:",Dna01," Dna10:",Dna10," Dna11:",Dna11
          ))
          print(paste0(
            " alpha:",alpha," beta:",beta," pab:",pab," Ef1:",Ef1," Ef2:",Ef2,
            " Ed1:",Ed1," Ed2:",Ed2," Es1:",Es1, " Es2:",Es2," Eu1:",Eu1," Eu2:",Eu2," eta:",eta
          ))
          stop("ERROR: probability or moments out of bounds")
        }
      }
      
      # Introduce population dynamics
      N <- x00 + x01 + x10 + x11
      rRelease <- Es * (rMob + rDup)
      #rLoss <- rUpt * N + rClear
      
      Dna_total <- Dna01 + Dna10 + Dna11
      rUpt_reg <- rUpt * Dna_total / (Dna_total + Dna_cst) # limit DNA uptake when DNA low
      
      # eDNA dynamics
      if(Dna_total > 0){
        dDna01 <- rRelease * (x01 + eta1 * x11) - rUpt_reg * N * Dna01 / Dna_total - rClear * Dna01
        dDna10 <- rRelease * (x10 + eta1 * x11) - rUpt_reg * N * Dna10 / Dna_total - rClear * Dna10
        dDna11 <- rRelease * (eta * x11)        - rUpt_reg * N * Dna11 / Dna_total - rClear * Dna11
        # Gene gains
        tau00_01 <- rUpt_reg * x00 * Dna01 / Dna_total
        tau00_10 <- rUpt_reg * x00 * Dna10 / Dna_total
        tau00_11 <- rUpt_reg * x00 * Dna11 / Dna_total
        tau01_11 <- rUpt_reg * x01 * Dna10 / Dna_total
        tau10_11 <- rUpt_reg * x10 * Dna01 / Dna_total
      } else {
        dDna01 <- rRelease * (x01 + eta1 * x11)
        dDna10 <- rRelease * (x10 + eta1 * x11)
        dDna11 <- rRelease * (eta * x11)
        tau00_01 <- 0
        tau00_10 <- 0
        tau00_11 <- 0
        tau01_11 <- 0
        tau10_11 <- 0
      }
      
      # Losses
      l <- Es * rMob
      tau01_00 <- l * x01
      tau10_00 <- l * x10
      tau11_00 <- l * eta  * x11
      tau11_01 <- l * eta1 * x11
      tau11_10 <- l * eta1 * x11
      
      # Variable death rates
      if(envId==1){ # constant antibiotic
        #dR <- deathRateFunc(Time)
        dR <- c(deathBasal+deathIncrease,deathBasal+deathIncrease,deathBasal+deathIncrease,deathBasal)
        damp <- 1 - Ntotal / K1 # logistic growth rate
      } else { # no antibiotic
        dR <- c(deathBasal,deathBasal,deathBasal,deathBasal)
        damp <- 1 - Ntotal / K2 # logistic growth rate
      }
      
      r00 <- (b00 * damp - dR[1])
      r01 <- (b01 * damp - dR[2])
      r10 <- (b10 * damp - dR[3])
      r11 <- (b11 * damp - dR[4])
      
      # Populations
      dx00 <- r00 * x00 + (tau01_00 + tau10_00 + tau11_00) - (tau00_01 + tau00_10 + tau00_11)
      dx01 <- r01 * x01 + (tau00_01 + tau11_01) - (tau01_11 + tau01_00)
      dx10 <- r10 * x10 + (tau00_10 + tau11_10) - (tau10_11 + tau10_00)
      dx11 <- r11 * x11 + (tau00_11 + tau01_11 + tau10_11) - (tau11_00 + tau11_01 + tau11_10)
      
      # Component derivatives (intermediate forms have underscore suffix)
      dNf_ <- r11 * x11 # unchanged arc distance
      dNd_ <- tau00_11 - tau11_00 # double gene events (gain & loss)
      dNs_ <- -(tau11_01 + tau11_10) # double gene to single gene
      dNu_ <- +(tau01_11 + tau10_11) # single gene to double gene
      
      dNEf1 <- dNf_ * Ef1 + dNd_ * Ed1 + dNs_ * Es1 + dNu_ * Eu1
      dNEf2 <- dNf_ * Ef2 + dNd_ * Ed2 + dNs_ * Es2 + dNu_ * Eu2
      
      if(envId==1){
        # Cell migration
        dx00 <- dx00 + G * (x00_env2 - x00) 
        dx01 <- dx01 + G * (x01_env2 - x01)
        dx10 <- dx10 + G * (x10_env2 - x10)
        dx11 <- dx11 + G * (x11_env2 - x11)
        # Arc distance modifications
        if(x11_env2 > 0) {
          Ef1_env2 <- NEf1_env2 / x11_env2 # Moments 1 and 2 (about zero) and parameters
          Ef2_env2 <- NEf2_env2 / x11_env2
        } else { # uniform arc distance
          Ef1_env2 <- 1/4
          Ef2_env2 <- 1/12
        }
        # dNin_ <- G * x11_env2 # x11 cells comming from env2
        # dNout_ <- G * x11 # x11 cells leaving env1
        # dNEf1 <- dNEf1 + dNin_ * Ef1_env2 - dNout_ * Ef1 # arc distance variation due to x11 cells from environment2
        # dNEf2 <- dNEf2 + dNin_ * Ef2_env2 - dNout_ * Ef2
        dNEf1 <- dNEf1 + G * (x11_env2 * Ef1_env2 - x11 * Ef1)
        dNEf2 <- dNEf2 + G * (x11_env2 * Ef2_env2 - x11 * Ef2)
      } # end env1 specs
      if(envId==2){
        # Cell migration
        dx00 <- dx00 + G * (x00_env1 - x00) 
        dx01 <- dx01 + G * (x01_env1 - x01)
        dx10 <- dx10 + G * (x10_env1 - x10)
        dx11 <- dx11 + G * (x11_env1 - x11)
        # Arc distance modifications
        if(x11_env1 > 0) {
          Ef1_env1 <- NEf1_env1 / x11_env1 # Moments 1 and 2 (about zero) and parameters
          Ef2_env1 <- NEf2_env1 / x11_env1
        } else { # uniform arc distance
          Ef1_env1 <- 1/4
          Ef2_env1 <- 1/12
        }
        # dNin_ <- G * x11_env1 # x11 cells comming from env1
        # dNout_ <- G * x11 # cells leaving env2
        # dNEf1 <- dNEf1 + dNin_ * Ef1_env1 - dNout_ * Ef1 # arc distance variation due to x11 cells from environment2
        # dNEf2 <- dNEf2 + dNin_ * Ef2_env1 - dNout_ * Ef2
        dNEf1 <- dNEf1 + G * (x11_env1 * Ef1_env1 - x11 * Ef1)
        dNEf2 <- dNEf2 + G * (x11_env1 * Ef2_env1 - x11 * Ef2)
      } # end env2 specs
      
      state_env <- c(dx00, dx01, dx10, dx11, dNEf1, dNEf2, dDna01, dDna10, dDna11)
      names(state_env) <- paste0(c("dx00", "dx01", "dx10", "dx11", "dNEf1", "dNEf2", "dDna01", "dDna10", "dDna11"), "_env", envId)
      assign(paste0("state_env",envId), state_env)
    } # end environment loop
    to_return <- c(state_env1, state_env2)
    
    return(list(
      to_return,
      #c(dx00, dx01, dx10, dx11, dNEf1, dNEf2, dDna01, dDna10, dDna11),
      # a = a, b = b,
      eta = eta#,
      # Ef1 = Ef1, Ed1 = Ed1, Es1 = Es1#,
      # x11gain = (tau00_11 + tau01_11 + tau10_11),
      # x11loss = (tau11_00 + tau11_01 + tau11_10),
      # x11conta = tau01_11 + tau10_11,
      # r01 = r01,
      # dx01 = dx01
    )) # end return
  }) # end with()
}

### 2 patch model - rUpt reg - migration equal ###
# migration env1->env2 and env2->env1 equal even if K1 and K2 different
LVcompet_2patchReg_bis <- function(Time, State, Pars, DEBUG = FALSE) { # State & Pars are named vectors with population i suffix _i
  # Number of species (pop)
  nSpecies <- Pars["nSpecies"]
  # Total population (all species)
  Ntotal_env1 <- sum(State[c(
    paste0("x00_env1_sp", 1:nSpecies),
    paste0("x01_env1_sp", 1:nSpecies),
    paste0("x10_env1_sp", 1:nSpecies),
    paste0("x11_env1_sp", 1:nSpecies))])
  Ntotal_env2 <- sum(State[c(
    paste0("x00_env2_sp", 1:nSpecies),
    paste0("x01_env2_sp", 1:nSpecies),
    paste0("x10_env2_sp", 1:nSpecies),
    paste0("x11_env2_sp", 1:nSpecies))])
  
  stepState <- State # to store the output of the simulation
  stepState <- sapply(stepState, function(x) 0)
  
  for(i in 1:nSpecies) {
    # Argument preparation
    State_pop_i <- State[names(stepState) %like% paste0("_sp",i)] # variable of species i
    names(State_pop_i) <-  gsub(paste0("_sp",i), "", names(State_pop_i)) # remove pop suffix for populationStep()
    Pars_pop_i <- Pars[names(Pars) %like% paste0("_sp",i)] # parameters of species i
    names(Pars_pop_i) <-  gsub(paste0("_sp",i), "", names(Pars_pop_i)) # remove pop suffix for populationStep()
    # Simulation step
    stepState[names(stepState) %like% paste0("_sp",i)] <- populationStep_2patchReg_bis(
      Time,
      State_pop_i,
      c(Pars_pop_i, Ntotal_env1 = Ntotal_env1, Ntotal_env2 = Ntotal_env2),
      DEBUG = DEBUG
    )[[1]]
  }
  
  return(list(
    stepState
  ))
}

populationStep_2patchReg_bis <- function(Time, State, Pars, DEBUG = FALSE) { # population step for one species
  with(as.list(c(State, Pars)), { # with(data, expression)
    
    # Extract parameters and variables (!for test only!)
    # for(i in 1:length(State)){
    #   assign(gsub("_sp1","",names(State[i])), State[i])
    # }
    # for(i in 1:length(Pars)){
    #   assign(gsub("_sp1","",names(Pars[i])), Pars[i])
    # }
    
    for(envId in 1:2){ # iterate through environments
      # Get state of environment envId
      x00 <- get(paste0("x00_env", envId))
      x01 <- get(paste0("x01_env", envId))
      x10 <- get(paste0("x10_env", envId))
      x11 <- get(paste0("x11_env", envId))
      NEf1 <- get(paste0("NEf1_env", envId))
      NEf2 <- get(paste0("NEf2_env", envId))
      Dna01 <- get(paste0("Dna01_env", envId))
      Dna10 <- get(paste0("Dna10_env", envId))
      Dna11 <- get(paste0("Dna11_env", envId))
      Ntotal <- get(paste0("Ntotal_env", envId))
      
      if(DEBUG){
        #if(x00<0 | x01<0 | x10<0 | x11<0 | Dna01<0 | Dna10<0 | Dna11<0){ # prevent negative populations
        if(x00<(-0.1) | x01<(-0.1) | x10<(-0.1) | x11<(-0.1) | Dna01<(-0.01) | Dna10<(-0.01) | Dna11<(-0.01)){ # prevent negative population (with small error marge)
          print(paste0(
            " time:",Time," env:",envId," x00:",x00," x01:",x01," x10:",x10," x11:",x11," Dna01:",Dna01," Dna10:",Dna10," Dna11:",Dna11
          ))
          stop("ERROR: cell or Dna < 0")
        }
      }
      
      # Moments 1 and 2 (about zero) and parameters
      if(x11 > 1E-6) { # >0 may lead to simulation failure when x11 tends to 0
        Ef1 <- NEf1 / x11
        Ef2 <- NEf2 / x11
      } else { # uniform arc distance
        Ef1 <- 1/4
        Ef2 <- 1/12
      }
      
      # Find arc distance distribution parameters
      ab <- as.double(absolve(Ef1, Ef2 - Ef1^2))
      # print(ab)
      alpha <- ab[1]
      beta  <- ab[2]
      
      # Probability of comobilization
      pab <- comobg(gamma, alpha, beta)
      
      # Conditional moments: double-gene event
      Md <- as.double(muvarg(gamma, alpha, beta))
      Ed1 <- Md[1]
      Ed2 <- Md[2] + Ed1^2
      
      # Conditional moments: single-gene deletion (OBSOLETE)
      # Es1 <- (Ef1 - pab * Ed1) / (1 - pab)
      # Es2 <- (Ef2 - pab * Ed2) / (1 - pab)
      
      # Conditional moments: single-gene insertion (uniform arc distance)
      Eu1 <- 1/4
      Eu2 <- 1/12
      
      # Intermediate values
      eta <- pab / Es # /Es is here to compensate Es in rRelease term
      eta1 <- 1 - eta
      
      # TENTATIVE CORRECTION 21-07-30, P(!A,B)
      # Conditional moments: single-gene deletion
      Es1 <- (Es * Ef1 - pab * Ed1) / (Es * eta1)
      Es2 <- (Es * Ef2 - pab * Ed2) / (Es * eta1)
      
      if(DEBUG){
        if(Ef1<0 | Ef1>0.2501 | Ef2<0 | Ef2>0.8334 | # >1/4 & >1/12 + error marge allowed
           Ed1<0 | Ed1>0.2501 | Ed2<0 | Ed2>0.8334 |
           #Es1<0 | Es1>1/4 | Es2<0 | Es2>1/12 |
           #Eu1<0 | Eu1>1/4 | Eu2<0 | Eu2>1/12 |
           pab<0 | pab>1 | eta<0 | eta>1){
          print(paste0(
            " time:",Time," env:",envId," x00:",x00," x01:",x01," x10:",x10," x11:",x11," Dna01:",Dna01," Dna10:",Dna10," Dna11:",Dna11
          ))
          print(paste0(
            " alpha:",alpha," beta:",beta," pab:",pab," Ef1:",Ef1," Ef2:",Ef2,
            " Ed1:",Ed1," Ed2:",Ed2," Es1:",Es1, " Es2:",Es2," Eu1:",Eu1," Eu2:",Eu2," eta:",eta
          ))
          stop("ERROR: probability or moments out of bounds")
        }
      }
      
      # Introduce population dynamics
      N <- x00 + x01 + x10 + x11
      rRelease <- Es * (rMob + rDup)
      #rLoss <- rUpt * N + rClear
      
      Dna_total <- Dna01 + Dna10 + Dna11
      rUpt_reg <- rUpt * Dna_total / (Dna_total + Dna_cst) # limit DNA uptake when DNA low
      
      # eDNA dynamics
      if(Dna_total > 0){
        dDna01 <- rRelease * (x01 + eta1 * x11) - rUpt_reg * N * Dna01 / Dna_total - rClear * Dna01
        dDna10 <- rRelease * (x10 + eta1 * x11) - rUpt_reg * N * Dna10 / Dna_total - rClear * Dna10
        dDna11 <- rRelease * (eta * x11)        - rUpt_reg * N * Dna11 / Dna_total - rClear * Dna11
        # Gene gains
        tau00_01 <- rUpt_reg * x00 * Dna01 / Dna_total
        tau00_10 <- rUpt_reg * x00 * Dna10 / Dna_total
        tau00_11 <- rUpt_reg * x00 * Dna11 / Dna_total
        tau01_11 <- rUpt_reg * x01 * Dna10 / Dna_total
        tau10_11 <- rUpt_reg * x10 * Dna01 / Dna_total
      } else {
        dDna01 <- rRelease * (x01 + eta1 * x11)
        dDna10 <- rRelease * (x10 + eta1 * x11)
        dDna11 <- rRelease * (eta * x11)
        tau00_01 <- 0
        tau00_10 <- 0
        tau00_11 <- 0
        tau01_11 <- 0
        tau10_11 <- 0
      }
      
      # Losses
      l <- Es * rMob
      tau01_00 <- l * x01
      tau10_00 <- l * x10
      tau11_00 <- l * eta  * x11
      tau11_01 <- l * eta1 * x11
      tau11_10 <- l * eta1 * x11
      
      # Variable death rates
      if(envId==1){ # constant antibiotic
        #dR <- deathRateFunc(Time)
        dR <- c(deathBasal+deathIncrease,deathBasal+deathIncrease,deathBasal+deathIncrease,deathBasal)
        damp <- 1 - Ntotal / K1 # logistic growth rate
      } else { # no antibiotic
        dR <- c(deathBasal,deathBasal,deathBasal,deathBasal)
        damp <- 1 - Ntotal / K2 # logistic growth rate
      }
      
      r00 <- (b00 * damp - dR[1])
      r01 <- (b01 * damp - dR[2])
      r10 <- (b10 * damp - dR[3])
      r11 <- (b11 * damp - dR[4])
      
      # Populations
      dx00 <- r00 * x00 + (tau01_00 + tau10_00 + tau11_00) - (tau00_01 + tau00_10 + tau00_11)
      dx01 <- r01 * x01 + (tau00_01 + tau11_01) - (tau01_11 + tau01_00)
      dx10 <- r10 * x10 + (tau00_10 + tau11_10) - (tau10_11 + tau10_00)
      dx11 <- r11 * x11 + (tau00_11 + tau01_11 + tau10_11) - (tau11_00 + tau11_01 + tau11_10)
      
      # Component derivatives (intermediate forms have underscore suffix)
      dNf_ <- r11 * x11 # unchanged arc distance
      dNd_ <- tau00_11 - tau11_00 # double gene events (gain & loss)
      dNs_ <- -(tau11_01 + tau11_10) # double gene to single gene
      dNu_ <- +(tau01_11 + tau10_11) # single gene to double gene
      
      dNEf1 <- dNf_ * Ef1 + dNd_ * Ed1 + dNs_ * Es1 + dNu_ * Eu1
      dNEf2 <- dNf_ * Ef2 + dNd_ * Ed2 + dNs_ * Es2 + dNu_ * Eu2
      
      if(envId==1){
        # Cell migration
        dx00 <- dx00 + G * K1/K2 * x00_env2 - G * x00
        dx01 <- dx01 + G * K1/K2 * x01_env2 - G * x01
        dx10 <- dx10 + G * K1/K2 * x10_env2 - G * x10
        dx11 <- dx11 + G * K1/K2 * x11_env2 - G * x11
        # Arc distance modifications
        if(x11_env2 > 0) {
          Ef1_env2 <- NEf1_env2 / x11_env2 # Moments 1 and 2 (about zero) and parameters
          Ef2_env2 <- NEf2_env2 / x11_env2
        } else { # uniform arc distance
          Ef1_env2 <- 1/4
          Ef2_env2 <- 1/12
        }
        #dNEf1 <- dNEf1 + G * (x11_env2 * Ef1_env2 - x11 * Ef1)
        #dNEf2 <- dNEf2 + G * (x11_env2 * Ef2_env2 - x11 * Ef2)
        dNEf1 <- dNEf1 + G * K1/K2 * x11_env2 * Ef1_env2 - G * x11 * Ef1
        dNEf2 <- dNEf2 + G * K1/K2 * x11_env2 * Ef2_env2 - G * x11 * Ef2
      } # end env1 specs
      if(envId==2){
        # Cell migration
        dx00 <- dx00 - G * K1/K2 * x00 + G * x00_env1
        dx01 <- dx01 - G * K1/K2 * x01 + G * x01_env1
        dx10 <- dx10 - G * K1/K2 * x10 + G * x10_env1
        dx11 <- dx11 - G * K1/K2 * x11 + G * x11_env1
        # Arc distance modifications
        if(x11_env1 > 0) {
          Ef1_env1 <- NEf1_env1 / x11_env1 # Moments 1 and 2 (about zero) and parameters
          Ef2_env1 <- NEf2_env1 / x11_env1
        } else { # uniform arc distance
          Ef1_env1 <- 1/4
          Ef2_env1 <- 1/12
        }
        dNEf1 <- dNEf1 + G * x11_env1 * Ef1_env1 - G * K1/K2 * x11 * Ef1
        dNEf2 <- dNEf2 + G * x11_env1 * Ef2_env1 - G * K1/K2 * x11 * Ef2
      } # end env2 specs
      
      state_env <- c(dx00, dx01, dx10, dx11, dNEf1, dNEf2, dDna01, dDna10, dDna11)
      names(state_env) <- paste0(c("dx00", "dx01", "dx10", "dx11", "dNEf1", "dNEf2", "dDna01", "dDna10", "dDna11"), "_env", envId)
      assign(paste0("state_env",envId), state_env)
    } # end environment loop
    to_return <- c(state_env1, state_env2)
    
    return(list(
      to_return,
      #c(dx00, dx01, dx10, dx11, dNEf1, dNEf2, dDna01, dDna10, dDna11),
      # a = a, b = b,
      eta = eta#,
      # Ef1 = Ef1, Ed1 = Ed1, Es1 = Es1#,
      # x11gain = (tau00_11 + tau01_11 + tau10_11),
      # x11loss = (tau11_00 + tau11_01 + tau11_10),
      # x11conta = tau01_11 + tau10_11,
      # r01 = r01,
      # dx01 = dx01
    )) # end return
  }) # end with()
}

### Lysis & Loss / no excision ###
LVcompet_lysis <- function(Time, State, Pars, DEBUG=FALSE) { # State & Pars are named vectors with population i suffix _i
  # Number of species (pop)
  nSpecies <- Pars["nSpecies"]
  # Total population (all species)
  Ntotal <- sum(State[c(
    paste0("x00_", 1:nSpecies),
    paste0("x01_", 1:nSpecies),
    paste0("x10_", 1:nSpecies),
    paste0("x11_", 1:nSpecies))])
  
  stepState <- State # to store the output of the simulation
  stepState <- sapply(stepState, function(x) 0)
  
  for(i in 1:nSpecies) {
    # Argument preparation
    State_pop_i <- State[names(stepState) %like% paste0("_",i)] # variable of species i
    names(State_pop_i) <-  gsub(paste0("_",i), "", names(State_pop_i)) # remove pop suffix for populationStep()
    Pars_pop_i <- Pars[names(Pars) %like% paste0("_",i)] # parameters of species i
    names(Pars_pop_i) <-  gsub(paste0("_",i), "", names(Pars_pop_i)) # remove pop suffix for populationStep()
    # Simulation step
    stepState[names(stepState) %like% paste0("_",i)] <- populationStep_lysis(
      Time,
      State_pop_i,
      c(Pars_pop_i, Ntotal = Ntotal),
      DEBUG=DEBUG
    )[[1]]
  }
  
  return(list(
    stepState
  ))
}

populationStep_lysis <- function(Time, State, Pars, DEBUG=FALSE) { # population step for one species
  with(as.list(c(State, Pars)), { # with(data, expression)
    
    if(DEBUG){
      #if(x00<0 | x01<0 | x10<0 | x11<0 | Dna01<0 | Dna10<0 | Dna11<0){ # prevent negative populations
      if(x00<(-0.1) | x01<(-0.1) | x10<(-0.1) | x11<(-0.1) | Dna01<(-0.01) | Dna10<(-0.01) | Dna11<(-0.01)){ # prevnt negative population (with small error marge)
        print(paste0(
          "x00:",x00," x01:",x01," x10:",x10," x11:",x11,
          " Dna01:",Dna01," Dna10:",Dna10," Dna11:",Dna11,
          " Ntotal:",Ntotal
        ))
        stop("ERROR: cell or Dna < 0")
      }
    }
    
    # Moments 1 and 2 (about zero) and parameters
    if(x11 > 1E-6) { # using x11>0 may lead to simulation failure when x11 tends to 0
      # ! arc distance results are not valid for small x11 population
      Ef1 <- NEf1 / x11
      Ef2 <- NEf2 / x11
    } else { # uniform arc distance
      Ef1 <- 1/4
      Ef2 <- 1/12
    }
    
    # Find arc distance distribution parameters
    ab <- as.double(absolve(Ef1, Ef2 - Ef1^2))
    alpha <- ab[1]
    beta  <- ab[2]
    
    # Probability of comobilization
    # comobg: unconditional probability of cotransfer
    pab <- comobg(gamma, alpha, beta)
    
    # Conditional moments: double-gene event
    # i.e. arc distance mean and variance conditional to a HGT event involving both genes
    Md <- as.double(muvarg(gamma, alpha, beta))
    Ed1 <- Md[1]
    Ed2 <- Md[2] + Ed1^2
    
    # Conditional moments: single-gene deletion, P(!(AB)) (!OBSOLETE!)
    #Es1 <- (Ef1 - pab * Ed1) / (1 - pab)
    #Es2 <- (Ef2 - pab * Ed2) / (1 - pab)
    
    # Conditional moments: single-gene insertion (uniform arc distance)
    Eu1 <- 1/4
    Eu2 <- 1/12
    
    # Intermediate values
    # Es: probability P(A)=P(B) that an HGT event involves A and/or B
    # eta: probability P(A|B) comobilization conditional to HGT event involves A and/or B
    # pab: probability comobilization when HGT event
    eta <- pab / Es # conditional probability comobilization
    eta1 <- 1 - eta # conditional probability single mobilization
    
    # CORRECTION 21-07-30, P(!A,B)
    # Conditional moments: single-gene deletion
    Es1 <- (Es * Ef1 - pab * Ed1) / (Es * eta1)
    Es2 <- (Es * Ef2 - pab * Ed2) / (Es * eta1)
    
    if(DEBUG){
      if(Ef1<0 | Ef1>0.2501 | Ef2<0 | Ef2>0.8334 | # >1/4 & >1/12 + error marge allowed
         Ed1<0 | Ed1>0.2501 | Ed2<0 | Ed2>0.8334 |
         #Es1<0 | Es1>1/4 | Es2<0 | Es2>1/12 |
         #Eu1<0 | Eu1>1/4 | Eu2<0 | Eu2>1/12 |
         pab<0 | pab>1 | eta<0 | eta>1){
        print(paste0(
          " time:",Time," x00:",x00," x01:",x01," x10:",x10," x11:",x11," Dna01:",Dna01," Dna10:",Dna10," Dna11:",Dna11
        ))
        print(paste0(
          " alpha:",alpha," beta:",beta," pab:",pab," Ef1:",Ef1," Ef2:",Ef2,
          " Ed1:",Ed1," Ed2:",Ed2," Es1:",Es1, " Es2:",Es2," Eu1:",Eu1," Eu2:",Eu2," eta:",eta
        ))
        stop("ERROR: probability or moments out of bounds")
      }
    }
    
    # Add x00 niche cells (constant x00 population)
    x00_tot <- x00 + x00_niche
    
    # Variable death rates
    dR <- deathRateFunc_tris(Time, db=deathBasal, di=deathIncrease, dp=deathPeriod, trt=treatmentType, lp=lifePeriod)
    
    # Introduce population dynamics
    N <- x00_tot + x01 + x10 + x11 # total cells of this species
    #rRelease <- Es * (rMob + rDup) # DNA release rate
    rRelease_01 <- Es * rDup + dR[2] # DNA release rate - no excision / release on cell lysis
    rRelease_10 <- Es * rDup + dR[3] # DNA release rate - no excision / release on cell lysis
    rRelease_11 <- Es * rDup + dR[4] # DNA release rate - no excision / release on cell lysis
    #rLoss <- rUpt * N + rClear # DNA removal rate
    Dna_total <- Dna01 + Dna10 + Dna11
    rUpt_reg <- rUpt * Dna_total / (Dna_total + Dna_cst) # limit DNA uptake when DNA low
    
    # eDNA dynamics
    if(Dna_total > 0){
      dDna01 <- rRelease_01 * (x01 + eta1 * x11) - rUpt_reg * N * Dna01 / Dna_total - rClear * Dna01
      dDna10 <- rRelease_10 * (x10 + eta1 * x11) - rUpt_reg * N * Dna10 / Dna_total - rClear * Dna10
      dDna11 <- rRelease_11 * (eta * x11)        - rUpt_reg * N * Dna11 / Dna_total - rClear * Dna11
      # Gene gains
      tau00_01 <- rUpt_reg * x00_tot * Dna01 / Dna_total
      tau00_10 <- rUpt_reg * x00_tot * Dna10 / Dna_total
      tau00_11 <- rUpt_reg * x00_tot * Dna11 / Dna_total
      tau01_11 <- rUpt_reg * x01 * Dna10 / Dna_total
      tau10_11 <- rUpt_reg * x10 * Dna01 / Dna_total
    } else {
      dDna01 <- rRelease_01 * (x01 + eta1 * x11)
      dDna10 <- rRelease_10 * (x10 + eta1 * x11)
      dDna11 <- rRelease_11 * (eta * x11)
      tau00_01 <- 0
      tau00_10 <- 0
      tau00_11 <- 0
      tau01_11 <- 0
      tau10_11 <- 0
    }
    
    # Gene losses
    l <- Es * rMob # rMob = MGE loss but eDNA is not released
    tau01_00 <- l * x01
    tau10_00 <- l * x10
    tau11_00 <- l * eta  * x11 # double-gene loss
    tau11_01 <- l * eta1 * x11 # single gene loss
    tau11_10 <- l * eta1 * x11 # single gene loss
    
    if(DEBUG){
      if(is.na(dR[1]) | is.na(dR[2]) | is.na(dR[3]) | is.na(dR[4])){
        stop("ERROR: death rate function output NA")
      }
    }
    
    # Bacterial growth
    damp <- 1 - Ntotal / K # logistic growth rate
    r00 <- (b00 * damp - dR[1]) # original formula
    r01 <- (b01 * damp - dR[2]) # original formula
    r10 <- (b10 * damp - dR[3]) # original formula
    r11 <- (b11 * damp - dR[4]) # original formula
    
    # Cell dynamics
    dx00 <- r00 * x00 + (tau01_00 + tau10_00 + tau11_00) - (tau00_01 + tau00_10 + tau00_11) + WTin * (K - x00)
    dx01 <- r01 * x01 + (tau00_01 + tau11_01) - (tau01_11 + tau01_00) - WTin * x01
    dx10 <- r10 * x10 + (tau00_10 + tau11_10) - (tau10_11 + tau10_00) - WTin * x10
    dx11 <- r11 * x11 + (tau00_11 + tau01_11 + tau10_11) - (tau11_00 + tau11_01 + tau11_10) - WTin * x11
    
    # Component derivatives (intermediate forms have underscore suffix)
    dNf_ <- (r11 - WTin) * x11 # unchanged arc distance
    dNd_ <- tau00_11 - tau11_00 # double gene gain and loss
    dNs_ <- -(tau11_01 + tau11_10) # double gene to single gene
    dNu_ <- +(tau01_11 + tau10_11) # single gene to double gene
    
    # Arc distance dynamics
    dNEf1 <- dNf_ * Ef1 + dNd_ * Ed1 + dNs_ * Es1 + dNu_ * Eu1
    dNEf2 <- dNf_ * Ef2 + dNd_ * Ed2 + dNs_ * Es2 + dNu_ * Eu2
    
    return(list(
      c(dx00, dx01, dx10, dx11, dNEf1, dNEf2, dDna01, dDna10, dDna11),
      # a = a, b = b,
      eta = eta#,
      # Ef1 = Ef1, Ed1 = Ed1, Es1 = Es1#,
      # x11gain = (tau00_11 + tau01_11 + tau10_11),
      # x11loss = (tau11_00 + tau11_01 + tau11_10),
      # x11conta = tau01_11 + tau10_11,
      # r01 = r01,
      # dx01 = dx01
    )) # end return
  }) # end with()
}

### Loss / loss does not affect m distribution ###
LVcompet_loss <- function(Time, State, Pars, DEBUG=FALSE) { # State & Pars are named vectors with population i suffix _i
  # Number of species (pop)
  nSpecies <- Pars["nSpecies"]
  # Total population (all species)
  Ntotal <- sum(State[c(
    paste0("x00_", 1:nSpecies),
    paste0("x01_", 1:nSpecies),
    paste0("x10_", 1:nSpecies),
    paste0("x11_", 1:nSpecies))])
  
  stepState <- State # to store the output of the simulation
  stepState <- sapply(stepState, function(x) 0)
  
  for(i in 1:nSpecies) {
    # Argument preparation
    State_pop_i <- State[names(stepState) %like% paste0("_",i)] # variable of species i
    names(State_pop_i) <-  gsub(paste0("_",i), "", names(State_pop_i)) # remove pop suffix for populationStep()
    Pars_pop_i <- Pars[names(Pars) %like% paste0("_",i)] # parameters of species i
    names(Pars_pop_i) <-  gsub(paste0("_",i), "", names(Pars_pop_i)) # remove pop suffix for populationStep()
    # Simulation step
    stepState[names(stepState) %like% paste0("_",i)] <- populationStep_loss(
      Time,
      State_pop_i,
      c(Pars_pop_i, Ntotal = Ntotal),
      DEBUG=DEBUG
    )[[1]]
  }
  
  return(list(
    stepState
  ))
}

populationStep_loss <- function(Time, State, Pars, DEBUG=FALSE) { # population step for one species
  with(as.list(c(State, Pars)), { # with(data, expression)
    
    if(DEBUG){
      #if(x00<0 | x01<0 | x10<0 | x11<0 | Dna01<0 | Dna10<0 | Dna11<0){ # prevent negative populations
      if(x00<(-0.1) | x01<(-0.1) | x10<(-0.1) | x11<(-0.1) | Dna01<(-0.01) | Dna10<(-0.01) | Dna11<(-0.01)){ # prevnt negative population (with small error marge)
        print(paste0(
          "x00:",x00," x01:",x01," x10:",x10," x11:",x11,
          " Dna01:",Dna01," Dna10:",Dna10," Dna11:",Dna11,
          " Ntotal:",Ntotal
        ))
        stop("ERROR: cell or Dna < 0")
      }
    }
    
    # Moments 1 and 2 (about zero) and parameters
    if(x11 > 1E-6) { # using x11>0 may lead to simulation failure when x11 tends to 0
      # ! arc distance results are not valid for small x11 population
      Ef1 <- NEf1 / x11
      Ef2 <- NEf2 / x11
    } else { # uniform arc distance
      Ef1 <- 1/4
      Ef2 <- 1/12
    }
    
    # Find arc distance distribution parameters
    ab <- as.double(absolve(Ef1, Ef2 - Ef1^2))
    alpha <- ab[1]
    beta  <- ab[2]
    
    # Probability of comobilization
    # comobg: unconditional probability of cotransfer
    pab <- comobg(gamma, alpha, beta)
    
    # Conditional moments: double-gene event
    # i.e. arc distance mean and variance conditional to a HGT event involving both genes
    Md <- as.double(muvarg(gamma, alpha, beta))
    Ed1 <- Md[1]
    Ed2 <- Md[2] + Ed1^2
    
    # Conditional moments: single-gene deletion, P(!(AB)) (!OBSOLETE!)
    #Es1 <- (Ef1 - pab * Ed1) / (1 - pab)
    #Es2 <- (Ef2 - pab * Ed2) / (1 - pab)
    
    # Conditional moments: single-gene insertion (uniform arc distance)
    Eu1 <- 1/4
    Eu2 <- 1/12
    
    # Intermediate values
    # Es: probability P(A)=P(B) that an HGT event involves A and/or B
    # eta: probability P(A|B) comobilization conditional to HGT event involves A and/or B
    # pab: probability comobilization when HGT event
    eta <- pab / Es # conditional probability comobilization
    eta1 <- 1 - eta # conditional probability single mobilization
    
    # CORRECTION 21-07-30, P(!A,B)
    # Conditional moments: single-gene deletion
    Es1 <- (Es * Ef1 - pab * Ed1) / (Es * eta1)
    Es2 <- (Es * Ef2 - pab * Ed2) / (Es * eta1)
    
    if(DEBUG){
      if(Ef1<0 | Ef1>0.2501 | Ef2<0 | Ef2>0.8334 | # >1/4 & >1/12 + error marge allowed
         Ed1<0 | Ed1>0.2501 | Ed2<0 | Ed2>0.8334 |
         #Es1<0 | Es1>1/4 | Es2<0 | Es2>1/12 |
         #Eu1<0 | Eu1>1/4 | Eu2<0 | Eu2>1/12 |
         pab<0 | pab>1 | eta<0 | eta>1){
        print(paste0(
          " time:",Time," x00:",x00," x01:",x01," x10:",x10," x11:",x11," Dna01:",Dna01," Dna10:",Dna10," Dna11:",Dna11
        ))
        print(paste0(
          " alpha:",alpha," beta:",beta," pab:",pab," Ef1:",Ef1," Ef2:",Ef2,
          " Ed1:",Ed1," Ed2:",Ed2," Es1:",Es1, " Es2:",Es2," Eu1:",Eu1," Eu2:",Eu2," eta:",eta
        ))
        stop("ERROR: probability or moments out of bounds")
      }
    }
    
    # Add x00 niche cells (constant x00 population)
    x00_tot <- x00 + x00_niche
    
    # Variable death rates
    dR <- deathRateFunc_tris(Time, db=deathBasal, di=deathIncrease, dp=deathPeriod, trt=treatmentType, lp=lifePeriod)
    
    # Introduce population dynamics
    N <- x00_tot + x01 + x10 + x11 # total cells of this species
    #rRelease <- Es * (rMob + rDup) # DNA release rate
    rRelease <- Es * rDup # rMob does not release DNA
    #rLoss <- rUpt * N + rClear # DNA removal rate
    Dna_total <- Dna01 + Dna10 + Dna11
    rUpt_reg <- rUpt * Dna_total / (Dna_total + Dna_cst) # limit DNA uptake when DNA low + limit max uptake to rUpt
    
    # eDNA dynamics
    if(Dna_total > 0){
      dDna01 <- rRelease * (x01 + eta1 * x11) - rUpt_reg * N * Dna01 / Dna_total - rClear * Dna01
      dDna10 <- rRelease * (x10 + eta1 * x11) - rUpt_reg * N * Dna10 / Dna_total - rClear * Dna10
      dDna11 <- rRelease * (eta * x11)        - rUpt_reg * N * Dna11 / Dna_total - rClear * Dna11
      # Gene gains
      tau00_01 <- rUpt_reg * x00_tot * Dna01 / Dna_total
      tau00_10 <- rUpt_reg * x00_tot * Dna10 / Dna_total
      tau00_11 <- rUpt_reg * x00_tot * Dna11 / Dna_total
      tau01_11 <- rUpt_reg * x01 * Dna10 / Dna_total
      tau10_11 <- rUpt_reg * x10 * Dna01 / Dna_total
    } else {
      dDna01 <- rRelease * (x01 + eta1 * x11)
      dDna10 <- rRelease * (x10 + eta1 * x11)
      dDna11 <- rRelease * (eta * x11)
      tau00_01 <- 0
      tau00_10 <- 0
      tau00_11 <- 0
      tau01_11 <- 0
      tau10_11 <- 0
    }
    
    # Gene losses
    l <- Es * rMob # rMob = MGE loss but eDNA is not released
    tau01_00 <- l * x01
    tau10_00 <- l * x10
    tau11_00 <- l * eta  * x11 # double-gene loss
    tau11_01 <- l * eta1 * x11 # single gene loss
    tau11_10 <- l * eta1 * x11 # single gene loss
    
    if(DEBUG){
      if(is.na(dR[1]) | is.na(dR[2]) | is.na(dR[3]) | is.na(dR[4])){
        stop("ERROR: death rate function output NA")
      }
    }
    
    # Bacterial growth
    damp <- 1 - Ntotal / K # logistic growth rate
    r00 <- (b00 * damp - dR[1]) # original formula
    r01 <- (b01 * damp - dR[2]) # original formula
    r10 <- (b10 * damp - dR[3]) # original formula
    r11 <- (b11 * damp - dR[4]) # original formula
    
    # Cell dynamics
    dx00 <- r00 * x00 + (tau01_00 + tau10_00 + tau11_00) - (tau00_01 + tau00_10 + tau00_11) + WTin * (K - x00)
    dx01 <- r01 * x01 + (tau00_01 + tau11_01) - (tau01_11 + tau01_00) - WTin * x01
    dx10 <- r10 * x10 + (tau00_10 + tau11_10) - (tau10_11 + tau10_00) - WTin * x10
    dx11 <- r11 * x11 + (tau00_11 + tau01_11 + tau10_11) - (tau11_00 + tau11_01 + tau11_10) - WTin * x11
    
    # Component derivatives (intermediate forms have underscore suffix)
    dNf_ <- (r11 - WTin) * x11 # unchanged arc distance
    dNd_ <- tau00_11 - tau11_00 # double gene gain and loss
    dNs_ <- -(tau11_01 + tau11_10) # double gene to single gene
    dNu_ <- +(tau01_11 + tau10_11) # single gene to double gene
    
    # Gene excision does not affect
    # dNf_ <- r11 * x11 - (tau11_00 + tau11_01 + tau11_10) # unchanged arc distance / growth and gene loss
    # dNd_ <- tau00_11 # double gene gain (reduce arc distance)
    # dNs_ <- 0 # single gene loss
    # dNu_ <- tau01_11 + tau10_11 # single gene to double gene (uniform arc distance)
    
    # Arc distance dynamics
    dNEf1 <- dNf_ * Ef1 + dNd_ * Ed1 + dNs_ * Es1 + dNu_ * Eu1
    dNEf2 <- dNf_ * Ef2 + dNd_ * Ed2 + dNs_ * Es2 + dNu_ * Eu2
    
    return(list(
      c(dx00, dx01, dx10, dx11, dNEf1, dNEf2, dDna01, dDna10, dDna11),
      # a = a, b = b,
      eta = eta#,
      # Ef1 = Ef1, Ed1 = Ed1, Es1 = Es1#,
      # x11gain = (tau00_11 + tau01_11 + tau10_11),
      # x11loss = (tau11_00 + tau11_01 + tau11_10),
      # x11conta = tau01_11 + tau10_11,
      # r01 = r01,
      # dx01 = dx01
    )) # end return
  }) # end with()
}

### 2 patch model - rUpt reg - migration equal ###
# migration env1->env2 and env2->env1 equal even if K1 and K2 different
# rMob -> loss without release
LVcompet_2patchReg_loss <- function(Time, State, Pars, DEBUG = FALSE) { # State & Pars are named vectors with population i suffix _i
  # Number of species (pop)
  nSpecies <- Pars["nSpecies"]
  # Total population (all species)
  Ntotal_env1 <- sum(State[c(
    paste0("x00_env1_sp", 1:nSpecies),
    paste0("x01_env1_sp", 1:nSpecies),
    paste0("x10_env1_sp", 1:nSpecies),
    paste0("x11_env1_sp", 1:nSpecies))])
  Ntotal_env2 <- sum(State[c(
    paste0("x00_env2_sp", 1:nSpecies),
    paste0("x01_env2_sp", 1:nSpecies),
    paste0("x10_env2_sp", 1:nSpecies),
    paste0("x11_env2_sp", 1:nSpecies))])
  
  stepState <- State # to store the output of the simulation
  stepState <- sapply(stepState, function(x) 0)
  
  for(i in 1:nSpecies) {
    # Argument preparation
    State_pop_i <- State[names(stepState) %like% paste0("_sp",i)] # variable of species i
    names(State_pop_i) <-  gsub(paste0("_sp",i), "", names(State_pop_i)) # remove pop suffix for populationStep()
    Pars_pop_i <- Pars[names(Pars) %like% paste0("_sp",i)] # parameters of species i
    names(Pars_pop_i) <-  gsub(paste0("_sp",i), "", names(Pars_pop_i)) # remove pop suffix for populationStep()
    # Simulation step
    stepState[names(stepState) %like% paste0("_sp",i)] <- populationStep_2patchReg_loss(
      Time,
      State_pop_i,
      c(Pars_pop_i, Ntotal_env1 = Ntotal_env1, Ntotal_env2 = Ntotal_env2),
      DEBUG = DEBUG
    )[[1]]
  }
  
  return(list(
    stepState
  ))
}

populationStep_2patchReg_loss <- function(Time, State, Pars, DEBUG = FALSE) { # population step for one species
  with(as.list(c(State, Pars)), { # with(data, expression)
    
    for(envId in 1:2){ # iterate through environments
      # Get state of environment envId
      x00 <- get(paste0("x00_env", envId))
      x01 <- get(paste0("x01_env", envId))
      x10 <- get(paste0("x10_env", envId))
      x11 <- get(paste0("x11_env", envId))
      NEf1 <- get(paste0("NEf1_env", envId))
      NEf2 <- get(paste0("NEf2_env", envId))
      Dna01 <- get(paste0("Dna01_env", envId))
      Dna10 <- get(paste0("Dna10_env", envId))
      Dna11 <- get(paste0("Dna11_env", envId))
      Ntotal <- get(paste0("Ntotal_env", envId))
      
      if(DEBUG){
        #if(x00<0 | x01<0 | x10<0 | x11<0 | Dna01<0 | Dna10<0 | Dna11<0){ # prevent negative populations
        if(x00<(-0.1) | x01<(-0.1) | x10<(-0.1) | x11<(-0.1) | Dna01<(-0.01) | Dna10<(-0.01) | Dna11<(-0.01)){ # prevent negative population (with small error marge)
          print(paste0(
            " time:",Time," env:",envId," x00:",x00," x01:",x01," x10:",x10," x11:",x11," Dna01:",Dna01," Dna10:",Dna10," Dna11:",Dna11
          ))
          stop("ERROR: cell or Dna < 0")
        }
      }
      
      # Moments 1 and 2 (about zero) and parameters
      if(x11 > 1E-6) { # >0 may lead to simulation failure when x11 tends to 0
        Ef1 <- NEf1 / x11
        Ef2 <- NEf2 / x11
      } else { # uniform arc distance
        Ef1 <- 1/4
        Ef2 <- 1/12
      }
      
      # Find arc distance distribution parameters
      ab <- as.double(absolve(Ef1, Ef2 - Ef1^2))
      # print(ab)
      alpha <- ab[1]
      beta  <- ab[2]
      
      # Probability of comobilization
      pab <- comobg(gamma, alpha, beta)
      
      # Conditional moments: double-gene event
      Md <- as.double(muvarg(gamma, alpha, beta))
      Ed1 <- Md[1]
      Ed2 <- Md[2] + Ed1^2
      
      # Conditional moments: single-gene deletion (OBSOLETE)
      # Es1 <- (Ef1 - pab * Ed1) / (1 - pab)
      # Es2 <- (Ef2 - pab * Ed2) / (1 - pab)
      
      # Conditional moments: single-gene insertion (uniform arc distance)
      Eu1 <- 1/4
      Eu2 <- 1/12
      
      # Intermediate values
      eta <- pab / Es # /Es is here to compensate Es in rRelease term
      eta1 <- 1 - eta
      
      # TENTATIVE CORRECTION 21-07-30, P(!A,B)
      # Conditional moments: single-gene deletion
      Es1 <- (Es * Ef1 - pab * Ed1) / (Es * eta1)
      Es2 <- (Es * Ef2 - pab * Ed2) / (Es * eta1)
      
      if(DEBUG){
        if(Ef1<0 | Ef1>0.2501 | Ef2<0 | Ef2>0.8334 | # >1/4 & >1/12 + error marge allowed
           Ed1<0 | Ed1>0.2501 | Ed2<0 | Ed2>0.8334 |
           #Es1<0 | Es1>1/4 | Es2<0 | Es2>1/12 |
           #Eu1<0 | Eu1>1/4 | Eu2<0 | Eu2>1/12 |
           pab<0 | pab>1 | eta<0 | eta>1){
          print(paste0(
            " time:",Time," env:",envId," x00:",x00," x01:",x01," x10:",x10," x11:",x11," Dna01:",Dna01," Dna10:",Dna10," Dna11:",Dna11
          ))
          print(paste0(
            " alpha:",alpha," beta:",beta," pab:",pab," Ef1:",Ef1," Ef2:",Ef2,
            " Ed1:",Ed1," Ed2:",Ed2," Es1:",Es1, " Es2:",Es2," Eu1:",Eu1," Eu2:",Eu2," eta:",eta
          ))
          stop("ERROR: probability or moments out of bounds")
        }
      }
      
      # Introduce population dynamics
      N <- x00 + x01 + x10 + x11
      #rRelease <- Es * (rMob + rDup)
      rRelease <- Es * rDup
      #rLoss <- rUpt * N + rClear
      
      Dna_total <- Dna01 + Dna10 + Dna11
      rUpt_reg <- rUpt * Dna_total / (Dna_total + Dna_cst) # limit DNA uptake when DNA low
      
      # eDNA dynamics
      if(Dna_total > 0){
        dDna01 <- rRelease * (x01 + eta1 * x11 /2) - rUpt_reg * N * Dna01 / Dna_total - rClear * Dna01
        dDna10 <- rRelease * (x10 + eta1 * x11 /2) - rUpt_reg * N * Dna10 / Dna_total - rClear * Dna10
        dDna11 <- rRelease * (eta * x11)        - rUpt_reg * N * Dna11 / Dna_total - rClear * Dna11
        # Gene gains
        tau00_01 <- rUpt_reg * x00 * Dna01 / Dna_total
        tau00_10 <- rUpt_reg * x00 * Dna10 / Dna_total
        tau00_11 <- rUpt_reg * x00 * Dna11 / Dna_total
        tau01_11 <- rUpt_reg * x01 * Dna10 / Dna_total
        tau10_11 <- rUpt_reg * x10 * Dna01 / Dna_total
      } else {
        dDna01 <- rRelease * (x01 + eta1 * x11 /2) # x11 releases bot 01 and 10 genes with probability 0.5
        dDna10 <- rRelease * (x10 + eta1 * x11 /2)
        dDna11 <- rRelease * (eta * x11)
        tau00_01 <- 0
        tau00_10 <- 0
        tau00_11 <- 0
        tau01_11 <- 0
        tau10_11 <- 0
      }
      
      # Losses
      l <- Es * rMob
      tau01_00 <- l * x01
      tau10_00 <- l * x10
      tau11_00 <- l * eta  * x11
      tau11_01 <- l * eta1 * x11
      tau11_10 <- l * eta1 * x11
      
      # Variable death rates
      if(envId==1){ # constant antibiotic
        #dR <- deathRateFunc(Time)
        dR <- c(deathBasal+deathIncrease,deathBasal+deathIncrease,deathBasal+deathIncrease,deathBasal)
        damp <- 1 - Ntotal / K1 # logistic growth rate
      } else { # no antibiotic
        dR <- c(deathBasal,deathBasal,deathBasal,deathBasal)
        damp <- 1 - Ntotal / K2 # logistic growth rate
      }
      
      r00 <- (b00 * damp - dR[1])
      r01 <- (b01 * damp - dR[2])
      r10 <- (b10 * damp - dR[3])
      r11 <- (b11 * damp - dR[4])
      
      # Populations
      dx00 <- r00 * x00 + (tau01_00 + tau10_00 + tau11_00) - (tau00_01 + tau00_10 + tau00_11)
      dx01 <- r01 * x01 + (tau00_01 + tau11_01) - (tau01_11 + tau01_00)
      dx10 <- r10 * x10 + (tau00_10 + tau11_10) - (tau10_11 + tau10_00)
      dx11 <- r11 * x11 + (tau00_11 + tau01_11 + tau10_11) - (tau11_00 + tau11_01 + tau11_10)
      
      # Component derivatives (intermediate forms have underscore suffix)
      dNf_ <- r11 * x11 # unchanged arc distance
      dNd_ <- tau00_11 - tau11_00 # double gene events (gain & loss)
      dNs_ <- -(tau11_01 + tau11_10) # double gene to single gene
      dNu_ <- +(tau01_11 + tau10_11) # single gene to double gene
      
      # dNf_ <- r11 * x11 - (tau11_00 + tau11_01 + tau11_10) # unchanged arc distance / growth and gene loss
      # dNd_ <- tau00_11 # double gene gain (reduce arc distance)
      # dNs_ <- 0 # single gene loss
      # dNu_ <- tau01_11 + tau10_11 # single gene to double gene (uniform arc distance)
      
      dNEf1 <- dNf_ * Ef1 + dNd_ * Ed1 + dNs_ * Es1 + dNu_ * Eu1
      dNEf2 <- dNf_ * Ef2 + dNd_ * Ed2 + dNs_ * Es2 + dNu_ * Eu2
      
      if(envId==1){
        # Cell migration
        dx00 <- dx00 + G * K1/K2 * x00_env2 - G * x00
        dx01 <- dx01 + G * K1/K2 * x01_env2 - G * x01
        dx10 <- dx10 + G * K1/K2 * x10_env2 - G * x10
        dx11 <- dx11 + G * K1/K2 * x11_env2 - G * x11
        # Arc distance modifications
        if(x11_env2 > 0) {
          Ef1_env2 <- NEf1_env2 / x11_env2 # Moments 1 and 2 (about zero) and parameters
          Ef2_env2 <- NEf2_env2 / x11_env2
        } else { # uniform arc distance
          Ef1_env2 <- 1/4
          Ef2_env2 <- 1/12
        }
        #dNEf1 <- dNEf1 + G * (x11_env2 * Ef1_env2 - x11 * Ef1)
        #dNEf2 <- dNEf2 + G * (x11_env2 * Ef2_env2 - x11 * Ef2)
        dNEf1 <- dNEf1 + G * K1/K2 * x11_env2 * Ef1_env2 - G * x11 * Ef1
        dNEf2 <- dNEf2 + G * K1/K2 * x11_env2 * Ef2_env2 - G * x11 * Ef2
      } # end env1 specs
      if(envId==2){
        # Cell migration
        dx00 <- dx00 - G * K1/K2 * x00 + G * x00_env1
        dx01 <- dx01 - G * K1/K2 * x01 + G * x01_env1
        dx10 <- dx10 - G * K1/K2 * x10 + G * x10_env1
        dx11 <- dx11 - G * K1/K2 * x11 + G * x11_env1
        # Arc distance modifications
        if(x11_env1 > 0) {
          Ef1_env1 <- NEf1_env1 / x11_env1 # Moments 1 and 2 (about zero) and parameters
          Ef2_env1 <- NEf2_env1 / x11_env1
        } else { # uniform arc distance
          Ef1_env1 <- 1/4
          Ef2_env1 <- 1/12
        }
        dNEf1 <- dNEf1 + G * x11_env1 * Ef1_env1 - G * K1/K2 * x11 * Ef1
        dNEf2 <- dNEf2 + G * x11_env1 * Ef2_env1 - G * K1/K2 * x11 * Ef2
      } # end env2 specs
      
      state_env <- c(dx00, dx01, dx10, dx11, dNEf1, dNEf2, dDna01, dDna10, dDna11)
      names(state_env) <- paste0(c("dx00", "dx01", "dx10", "dx11", "dNEf1", "dNEf2", "dDna01", "dDna10", "dDna11"), "_env", envId)
      assign(paste0("state_env",envId), state_env)
    } # end environment loop
    to_return <- c(state_env1, state_env2)
    
    return(list(
      to_return,
      #c(dx00, dx01, dx10, dx11, dNEf1, dNEf2, dDna01, dDna10, dDna11),
      # a = a, b = b,
      eta = eta#,
      # Ef1 = Ef1, Ed1 = Ed1, Es1 = Es1#,
      # x11gain = (tau00_11 + tau01_11 + tau10_11),
      # x11loss = (tau11_00 + tau11_01 + tau11_10),
      # x11conta = tau01_11 + tau10_11,
      # r01 = r01,
      # dx01 = dx01
    )) # end return
  }) # end with()
}

### 2 patch model - rUpt reg - migration equal ###
# migration env1->env2 and env2->env1 equal even if K1 and K2 different
# rMob -> loss without release
# spaceBis: antibiotic1 in patch1 and antibiotic2 in patch2
LVcompet_2patchReg_loss_spaceBis <- function(Time, State, Pars, DEBUG = FALSE) { # State & Pars are named vectors with population i suffix _i
  # Number of species (pop)
  nSpecies <- Pars["nSpecies"]
  # Total population (all species)
  Ntotal_env1 <- sum(State[c(
    paste0("x00_env1_sp", 1:nSpecies),
    paste0("x01_env1_sp", 1:nSpecies),
    paste0("x10_env1_sp", 1:nSpecies),
    paste0("x11_env1_sp", 1:nSpecies))])
  Ntotal_env2 <- sum(State[c(
    paste0("x00_env2_sp", 1:nSpecies),
    paste0("x01_env2_sp", 1:nSpecies),
    paste0("x10_env2_sp", 1:nSpecies),
    paste0("x11_env2_sp", 1:nSpecies))])
  
  stepState <- State # to store the output of the simulation
  stepState <- sapply(stepState, function(x) 0)
  
  for(i in 1:nSpecies) {
    # Argument preparation
    State_pop_i <- State[names(stepState) %like% paste0("_sp",i)] # variable of species i
    names(State_pop_i) <-  gsub(paste0("_sp",i), "", names(State_pop_i)) # remove pop suffix for populationStep()
    Pars_pop_i <- Pars[names(Pars) %like% paste0("_sp",i)] # parameters of species i
    names(Pars_pop_i) <-  gsub(paste0("_sp",i), "", names(Pars_pop_i)) # remove pop suffix for populationStep()
    # Simulation step
    stepState[names(stepState) %like% paste0("_sp",i)] <- populationStep_2patchReg_loss_spaceBis(
      Time,
      State_pop_i,
      c(Pars_pop_i, Ntotal_env1 = Ntotal_env1, Ntotal_env2 = Ntotal_env2),
      DEBUG = DEBUG
    )[[1]]
  }
  
  return(list(
    stepState
  ))
}

populationStep_2patchReg_loss_spaceBis <- function(Time, State, Pars, DEBUG = FALSE) { # population step for one species
  with(as.list(c(State, Pars)), { # with(data, expression)
    
    for(envId in 1:2){ # iterate through environments
      # Get state of environment envId
      x00 <- get(paste0("x00_env", envId))
      x01 <- get(paste0("x01_env", envId))
      x10 <- get(paste0("x10_env", envId))
      x11 <- get(paste0("x11_env", envId))
      NEf1 <- get(paste0("NEf1_env", envId))
      NEf2 <- get(paste0("NEf2_env", envId))
      Dna01 <- get(paste0("Dna01_env", envId))
      Dna10 <- get(paste0("Dna10_env", envId))
      Dna11 <- get(paste0("Dna11_env", envId))
      Ntotal <- get(paste0("Ntotal_env", envId))
      
      if(DEBUG){
        #if(x00<0 | x01<0 | x10<0 | x11<0 | Dna01<0 | Dna10<0 | Dna11<0){ # prevent negative populations
        if(x00<(-0.1) | x01<(-0.1) | x10<(-0.1) | x11<(-0.1) | Dna01<(-0.01) | Dna10<(-0.01) | Dna11<(-0.01)){ # prevent negative population (with small error marge)
          print(paste0(
            " time:",Time," env:",envId," x00:",x00," x01:",x01," x10:",x10," x11:",x11," Dna01:",Dna01," Dna10:",Dna10," Dna11:",Dna11
          ))
          stop("ERROR: cell or Dna < 0")
        }
      }
      
      # Moments 1 and 2 (about zero) and parameters
      if(x11 > 1E-6) { # >0 may lead to simulation failure when x11 tends to 0
        Ef1 <- NEf1 / x11
        Ef2 <- NEf2 / x11
      } else { # uniform arc distance
        Ef1 <- 1/4
        Ef2 <- 1/12
      }
      
      # Find arc distance distribution parameters
      ab <- as.double(absolve(Ef1, Ef2 - Ef1^2))
      # print(ab)
      alpha <- ab[1]
      beta  <- ab[2]
      
      # Probability of comobilization
      pab <- comobg(gamma, alpha, beta)
      
      # Conditional moments: double-gene event
      Md <- as.double(muvarg(gamma, alpha, beta))
      Ed1 <- Md[1]
      Ed2 <- Md[2] + Ed1^2
      
      # Conditional moments: single-gene deletion (OBSOLETE)
      # Es1 <- (Ef1 - pab * Ed1) / (1 - pab)
      # Es2 <- (Ef2 - pab * Ed2) / (1 - pab)
      
      # Conditional moments: single-gene insertion (uniform arc distance)
      Eu1 <- 1/4
      Eu2 <- 1/12
      
      # Intermediate values
      eta <- pab / Es # /Es is here to compensate Es in rRelease term
      eta1 <- 1 - eta
      
      # TENTATIVE CORRECTION 21-07-30, P(!A,B)
      # Conditional moments: single-gene deletion
      Es1 <- (Es * Ef1 - pab * Ed1) / (Es * eta1)
      Es2 <- (Es * Ef2 - pab * Ed2) / (Es * eta1)
      
      if(DEBUG){
        if(Ef1<0 | Ef1>0.2501 | Ef2<0 | Ef2>0.8334 | # >1/4 & >1/12 + error marge allowed
           Ed1<0 | Ed1>0.2501 | Ed2<0 | Ed2>0.8334 |
           #Es1<0 | Es1>1/4 | Es2<0 | Es2>1/12 |
           #Eu1<0 | Eu1>1/4 | Eu2<0 | Eu2>1/12 |
           pab<0 | pab>1 | eta<0 | eta>1){
          print(paste0(
            " time:",Time," env:",envId," x00:",x00," x01:",x01," x10:",x10," x11:",x11," Dna01:",Dna01," Dna10:",Dna10," Dna11:",Dna11
          ))
          print(paste0(
            " alpha:",alpha," beta:",beta," pab:",pab," Ef1:",Ef1," Ef2:",Ef2,
            " Ed1:",Ed1," Ed2:",Ed2," Es1:",Es1, " Es2:",Es2," Eu1:",Eu1," Eu2:",Eu2," eta:",eta
          ))
          stop("ERROR: probability or moments out of bounds")
        }
      }
      
      # Introduce population dynamics
      N <- x00 + x01 + x10 + x11
      #rRelease <- Es * (rMob + rDup)
      rRelease <- Es * rDup
      #rLoss <- rUpt * N + rClear
      
      Dna_total <- Dna01 + Dna10 + Dna11
      rUpt_reg <- rUpt * Dna_total / (Dna_total + Dna_cst) # limit DNA uptake when DNA low
      
      # eDNA dynamics
      if(Dna_total > 0){
        dDna01 <- rRelease * (x01 + eta1 * x11 /2) - rUpt_reg * N * Dna01 / Dna_total - rClear * Dna01
        dDna10 <- rRelease * (x10 + eta1 * x11 /2) - rUpt_reg * N * Dna10 / Dna_total - rClear * Dna10
        dDna11 <- rRelease * (eta * x11)        - rUpt_reg * N * Dna11 / Dna_total - rClear * Dna11
        # Gene gains
        tau00_01 <- rUpt_reg * x00 * Dna01 / Dna_total
        tau00_10 <- rUpt_reg * x00 * Dna10 / Dna_total
        tau00_11 <- rUpt_reg * x00 * Dna11 / Dna_total
        tau01_11 <- rUpt_reg * x01 * Dna10 / Dna_total
        tau10_11 <- rUpt_reg * x10 * Dna01 / Dna_total
      } else {
        dDna01 <- rRelease * (x01 + eta1 * x11 /2) # x11 releases bot 01 and 10 genes with probability 0.5
        dDna10 <- rRelease * (x10 + eta1 * x11 /2)
        dDna11 <- rRelease * (eta * x11)
        tau00_01 <- 0
        tau00_10 <- 0
        tau00_11 <- 0
        tau01_11 <- 0
        tau10_11 <- 0
      }
      
      # Losses
      l <- Es * rMob
      tau01_00 <- l * x01
      tau10_00 <- l * x10
      tau11_00 <- l * eta  * x11
      tau11_01 <- l * eta1 * x11
      tau11_10 <- l * eta1 * x11
      
      # Variable death rates
      if(envId==1){ # antibiotic 1
        #dR <- deathRateFunc(Time)
        dR <- c(deathBasal+deathIncrease, deathBasal+deathIncrease, deathBasal, deathBasal)
        damp <- 1 - Ntotal / K1 # logistic growth rate
      } else { # antibiotic 2
        dR <- c(deathBasal+deathIncrease, deathBasal, deathBasal+deathIncrease, deathBasal)
        damp <- 1 - Ntotal / K2 # logistic growth rate
      }
      
      r00 <- (b00 * damp - dR[1])
      r01 <- (b01 * damp - dR[2])
      r10 <- (b10 * damp - dR[3])
      r11 <- (b11 * damp - dR[4])
      
      # Populations
      dx00 <- r00 * x00 + (tau01_00 + tau10_00 + tau11_00) - (tau00_01 + tau00_10 + tau00_11)
      dx01 <- r01 * x01 + (tau00_01 + tau11_01) - (tau01_11 + tau01_00)
      dx10 <- r10 * x10 + (tau00_10 + tau11_10) - (tau10_11 + tau10_00)
      dx11 <- r11 * x11 + (tau00_11 + tau01_11 + tau10_11) - (tau11_00 + tau11_01 + tau11_10)
      
      # Component derivatives (intermediate forms have underscore suffix)
      dNf_ <- r11 * x11 # unchanged arc distance
      dNd_ <- tau00_11 - tau11_00 # double gene events (gain & loss)
      dNs_ <- -(tau11_01 + tau11_10) # double gene to single gene
      dNu_ <- +(tau01_11 + tau10_11) # single gene to double gene
      
      # dNf_ <- r11 * x11 - (tau11_00 + tau11_01 + tau11_10) # unchanged arc distance / growth and gene loss
      # dNd_ <- tau00_11 # double gene gain (reduce arc distance)
      # dNs_ <- 0 # single gene loss
      # dNu_ <- tau01_11 + tau10_11 # single gene to double gene (uniform arc distance)
      
      dNEf1 <- dNf_ * Ef1 + dNd_ * Ed1 + dNs_ * Es1 + dNu_ * Eu1
      dNEf2 <- dNf_ * Ef2 + dNd_ * Ed2 + dNs_ * Es2 + dNu_ * Eu2
      
      if(envId==1){
        # Cell migration
        dx00 <- dx00 + G * K1/K2 * x00_env2 - G * x00
        dx01 <- dx01 + G * K1/K2 * x01_env2 - G * x01
        dx10 <- dx10 + G * K1/K2 * x10_env2 - G * x10
        dx11 <- dx11 + G * K1/K2 * x11_env2 - G * x11
        # Arc distance modifications
        if(x11_env2 > 0) {
          Ef1_env2 <- NEf1_env2 / x11_env2 # Moments 1 and 2 (about zero) and parameters
          Ef2_env2 <- NEf2_env2 / x11_env2
        } else { # uniform arc distance
          Ef1_env2 <- 1/4
          Ef2_env2 <- 1/12
        }
        #dNEf1 <- dNEf1 + G * (x11_env2 * Ef1_env2 - x11 * Ef1)
        #dNEf2 <- dNEf2 + G * (x11_env2 * Ef2_env2 - x11 * Ef2)
        dNEf1 <- dNEf1 + G * K1/K2 * x11_env2 * Ef1_env2 - G * x11 * Ef1
        dNEf2 <- dNEf2 + G * K1/K2 * x11_env2 * Ef2_env2 - G * x11 * Ef2
      } # end env1 specs
      if(envId==2){
        # Cell migration
        dx00 <- dx00 - G * K1/K2 * x00 + G * x00_env1
        dx01 <- dx01 - G * K1/K2 * x01 + G * x01_env1
        dx10 <- dx10 - G * K1/K2 * x10 + G * x10_env1
        dx11 <- dx11 - G * K1/K2 * x11 + G * x11_env1
        # Arc distance modifications
        if(x11_env1 > 0) {
          Ef1_env1 <- NEf1_env1 / x11_env1 # Moments 1 and 2 (about zero) and parameters
          Ef2_env1 <- NEf2_env1 / x11_env1
        } else { # uniform arc distance
          Ef1_env1 <- 1/4
          Ef2_env1 <- 1/12
        }
        dNEf1 <- dNEf1 + G * x11_env1 * Ef1_env1 - G * K1/K2 * x11 * Ef1
        dNEf2 <- dNEf2 + G * x11_env1 * Ef2_env1 - G * K1/K2 * x11 * Ef2
      } # end env2 specs
      
      state_env <- c(dx00, dx01, dx10, dx11, dNEf1, dNEf2, dDna01, dDna10, dDna11)
      names(state_env) <- paste0(c("dx00", "dx01", "dx10", "dx11", "dNEf1", "dNEf2", "dDna01", "dDna10", "dDna11"), "_env", envId)
      assign(paste0("state_env",envId), state_env)
    } # end environment loop
    to_return <- c(state_env1, state_env2)
    
    return(list(
      to_return,
      #c(dx00, dx01, dx10, dx11, dNEf1, dNEf2, dDna01, dDna10, dDna11),
      # a = a, b = b,
      eta = eta#,
      # Ef1 = Ef1, Ed1 = Ed1, Es1 = Es1#,
      # x11gain = (tau00_11 + tau01_11 + tau10_11),
      # x11loss = (tau11_00 + tau11_01 + tau11_10),
      # x11conta = tau01_11 + tau10_11,
      # r01 = r01,
      # dx01 = dx01
    )) # end return
  }) # end with()
}



#### TESTING ####
# State_pop_i <- State[names(stepState) %like% paste0("_",i)] # variable of species i
# names(State_pop_i) <-  gsub(paste0("_",i), "", names(State_pop_i)) # remove pop suffix for populationStep()
# Pars_pop_i <- Pars[names(Pars) %like% paste0("_",i)] # parameters of species i
# names(Pars_pop_i) <-  gsub(paste0("_",i), "", names(Pars_pop_i)) # remove pop suffix for populationStep()
# # names vectors to variable to test equtions inside function
# for(it in 1:length(State_pop_i)){
#   assign(names(State_pop_i)[it], State_pop_i[it])
# }
# for(it in 1:length(Pars_pop_i)){
#   assign(names(Pars_pop_i)[it], Pars_pop_i[it])
# }


########################################################
## ---------------- MANDENA MICROCEBE --------------- ##
## ---- CAPTURE-MARK_RECAPTURE SURVIVAL ANALYSIS ---- ##
########################################################

rm(list = ls())


## ------ LIBRARIES ------

library(coda)
library(nimble)
library(xtable)
library(data.table)
library(ggplot2)



## -----------------------------------------------------------------------------
## ------ I. FIT NIMBLE MODEL ------

## ------   1. LOAD MANDENA CMR DATA ------

load("./data/Microcebe_Mandena_data.RData")



## ------   2. WRITE NIMBLE MODEL ------

nimModel <- nimbleCode({
  
  ## DEMOGRAPHIC PROCESS
  psi ~ dunif(0,1)
  
  for(d in 1:2){
    for(s in 1:2){
      z.time[d,s] ~ dbern(psi)
      z.temp[d,s] ~ dbern(psi)
      z.transloc[d,s] ~ dbern(psi)
      logit.phi0[d,s] ~ dnorm(0,0.01)     ## Baseline survival[site,sex]
      beta.time[d,s] ~ dnorm(0,0.01)      ## Temporal effect[site,sex]
      beta.temp[d,s] ~ dnorm(0,0.01)      ## Temperature effect[site,sex]
      beta.transloc[d,s] ~ dnorm(0,0.01)  ## Translocation effect[site,sex]
      for(m in 1:n.months){
        logit(PHI[d,s,1,m]) <- logit.phi0[d,s] +
          beta.temp[d,s] * z.temp[d,s] * temp[m] +
          beta.time[d,s] * z.time[d,s] * month[m]
        logit(PHI[d,s,2,m]) <- logit.phi0[d,s] +
          beta.temp[d,s] * z.temp[d,s] * temp[m] +
          beta.time[d,s] * z.time[d,s] * month[m] +
          beta.transloc[d,s] * z.transloc[d,s]
      }#m
    }#ss
  }#s
  
  
  ## Individual state
  for(i in 1:n.individuals){
    z[i,1] ~ dbern(1)
    for(t in 1:n.intervals[i]){
      phi[i,t] <- prod(PHI[status[i],sex[i],transloc[i],start.int[i,t]:end.int[i,t]])
      z[i,t+1] ~ dbern(z[i,t] * phi[i,t])
    }#t
  }#i
  
  
  ## DETECTION PROCESS
  for(f in 1:n.sites){
    for(s in 1:2){
      lambda0[f,s] ~ dnorm(0,0.01)    ## Detection Hazard rate female/disturbed
      gamma.temp[f,s] ~ dnorm(0,0.01) ## Detection Hazard rate female/disturbed
      z.det[f,s] ~ dbern(psi)
    }
  }
  
  ## Individual detection
  for(i in 1:n.individuals){
    for(t in 2:n.sessions[i]){
      log(lambda[i,t]) <- lambda0[site[i],sex[i]] +
        gamma.temp[site[i],sex[i]] * temp2[i,t] * z.det[site[i],sex[i]]
      p[i,t] <- 1-exp(-lambda[i,t] * duration[i,t])
      y[i,t] ~ dbern(p[i,t] * z[i,t])
    }#t
  }#i
}) 



## ------   3. FIT MODEL TO DATA ------
##-- Create a NIMBLE model object
Rmodel <- nimbleModel( code = nimModel,
                       constants = nimConstants,
                       data = nimData,
                       inits = nimInits,
                       calculate = F)
Rmodel$calculate()

##-- Configure and Build RJ-MCMC objects
conf <- configureMCMC(Rmodel,
                      monitors = nimParams,
                      print = FALSE)

configureRJ(conf,
            targetNodes = c("beta.temp"),
            indicatorNodes = c('z.temp'),
            control = list(mean = 0, scale = .2))


configureRJ(conf,
            targetNodes = c( "beta.time"),
            indicatorNodes = c('z.time'),
            control = list(mean = 0, scale = .2))

configureRJ(conf,
            targetNodes = c("beta.transloc"),
            indicatorNodes = c('z.transloc'),
            control = list(mean = 0, scale = .2))

configureRJ(conf,
            targetNodes = c("gamma.temp"),
            indicatorNodes = c('z.det'),
            control = list(mean = 0, scale = .2))

##-- Check the assigned samplers
conf$printSamplers(nimParams)

##-- Build MCMC
Rmcmc <- buildMCMC(conf)

##-- Compile and Run MCMC
## Finally, we compile both the model and MCMC objects and execute the
## compiled MCMC for 50 000 iterations and 4 chains.
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
MCMC_runtime <- system.time(
  nimOutput <- runMCMC( Cmcmc,
                        niter = 1000,
                        nburnin = 0,
                        nchains = 3,
                        thin = 1,
                        samplesAsCodaMCMC = T))

##-- Traceplots
plot(nimOutput)


## -----------------------------------------------------------------------------

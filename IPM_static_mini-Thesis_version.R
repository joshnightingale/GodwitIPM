# IPM in parallel 

#### load packages ####
library(nimble)
library(nimbleEcology)
library(magrittr)
library(doParallel) # also attaches 'foreach' and 'parallel'



#### Load data ####
load("IPM_CONSTANTS.RData")
list2env(IPM_CONSTANTS, globalenv())
IPM_CONSTANTS$constraint_data_s <- 1 # only need one

load("IPM_DATA.RData")

load("IPM_INITS.RData")
# names(IPM_INITS)[[5]] <- "R"

#### run settings ####
nt <- 50
nt2 <- 100
nb <- 100000
ni <- 100000 + nb
nc <- 1


#### Model code ####
IPM_code <- nimbleCode({
  
  #### CJS model for ringing data ####
  
  for (j in 1:N.age){ # separate estimates for each age
    ### survival
    surv.exp[j] <- expit(mu.surv[j])
  } # close j
  
  for (j in 2:N.age){ # separate estimates for each age
    # Regularising priors on independent betas for geographical effects
    beta.LAT[j] ~ dnorm(0, 0.1) # age-dep latitude effect on survival
    beta.OLD[j] ~ dnorm(0, 0.1) # age-dep colonisation effect on survival
    
  } # close j
  
  mu.surv.rec ~ dgamma(0.5, 5) # penalty term for birds marked as chicks
  mu.surv[1] <- mu.surv[2] - mu.surv.rec
  
  mu.surv[2] ~ dnorm(0, sd=1) # Vague priors for mean age-dependent survival
  mu.surv[3] ~ dnorm(0, sd=1) # Vague priors for mean age-dependent survival
  constraint_data_s ~ dconstraint(mu.surv[2] <= mu.surv[3])
  
  # betas for birds marked as chicks are equal to those ringed in 1st winter
  beta.LAT[1] <- beta.LAT[2]
  beta.OLD[1] <- beta.OLD[2]
  
  
  
  
  
  
  ### hyperprior for detection
  p.const ~ dnorm(0, 1)  # Vague prior for mean p
  mu.p <- expit(p.const) # logit scale reporting
  p.sd ~ dunif(0, 5)
  
  
  for (i in 1:N.ind){
    ## detection - individual random effect
    logit(p[i]) ~ dnorm(p.const, sd=p.sd)
    
    for (t in FIRST[i]:(N.years-1)){
      # age-dependent survival with latitude and colonisation effects
      phi[i,t] <- expit( mu.surv[AGE[i,t]] + ( KNOWN[i, t] * ( 
        (beta.LAT[AGE[i,t]] * LAT[i,t]) +
          (beta.OLD[AGE[i,t]] * OLD[i,t]) ) )
      )
    } #t
    
    # Likelihood
    y[i, FIRST[i]:N.years] ~ dCJS_vs(probSurvive = (phi[i, FIRST[i]:(N.years)]), 
                                     probCapture = p[i], 
                                     len = length(FIRST[i]:N.years))
    
  } #i
  
  
  
  #### Gaussian state-space / DM model  for Count Data ####
  
  #### priors for sites
  for (i in 1:Nsite) {  
    # sd of observation process
    obs.sd[i] ~ dgamma(3,1) 
    
    # sd of first N
    n1.sd[i] ~dgamma(1,1)
    
    ### site-specific survival and  recruitment by decade 
    
    ### site-specific age-specific survival
    s.ad.site[i] <- expit( mu.surv[3] + beta.LAT[3]*LAT.c[i] + beta.OLD[3]*OLD.c[i] )
    s.j.site[i] <- expit( mu.surv[2] + beta.LAT[2]*LAT.c[i] + beta.OLD[2]*OLD.c[i] )
    
    ### site-specific recruitment
    r.site.mu.norm[i] ~ dnorm(r.mean.norm, r.site.sd)
    r.site.mu[i] <- expit( r.site.mu.norm[i] )
    
    r.site[i] <- expit ( r.site.mu.norm[i] + r.beta.lat*LAT.c[i] + r.beta.col*OLD.c[i] )
    
  } # close i
  
  
  ##### hyperpriors
  
  ### recruitment
  
  ## recruitment per decade
  r.mean.norm ~ dnorm(0, 1) # overall mean recruitment prob
  r.mean <- expit(r.mean.norm)
  
  r.site.sd ~ dgamma(1, 1) # sd for site estimates
  
  r.proc.sd ~ dgamma(0.1, 1) # sd for recruitment process
  s.proc.sd ~ dgamma(1, 1) # sd for survival process
  
  
  # recruitment coefs
  r.beta.lat ~ dnorm(0, sd=0.1)
  r.beta.col ~ dnorm(0, sd=0.1)
  
  
  
  #### Biological process
  for(i in 1:Nsite) {
    
    ## first count
    N[i, 1] ~ T(dnorm(Y.first[i], sd=n1.sd[i]), 0, )
    
    # fill params from first year to avoid NA error
    R[i, 1] <- N[i,1] / 1
    S[i, 1] <- N[i,1] / 2
    
    
    ## subsequent counts
    for (t in 2:N.years) {
      
      ## survival process
      S[i, t] ~ T(dnorm( mean = ( N[i, (t-1)] * s.ad.site[i] ), 
                         sd = ( s.proc.sd * (N[i, (t-1)] + 1) ) ), 0, ) 
      
      ## recruitment process (inc. juv. survival)
      R[i, t] ~ T(dnorm( mean= (N[i, (t-1)] * r.site[i] * s.j.site[i]), 
                         sd = ( r.proc.sd * (N[i, (t-1)] + 1) ) ), 0, ) 
      
      
      ## sum to give total population
      N[i, t] <- S[i, t] + R[i, t]
      
      ## R and S cannot be negative!
      # constraint_data_r[i,t] ~ dconstraint( R[i,t] >= 0 & S[i,t] >= 0 )
    } # close i
  } # close t
  
  
  #### Observation process
  for(i in 1:Nsite) {
    for (t in 1:N.years) {
      for (j in 1:J[i,t]) { # use replicate counts where available
        Y.c[i, t, j] ~ dnorm(N[i, t], sd=obs.sd[i])
      } # close j
    } # close i
  } # close t
})



#### make cluster ####

ncore <- 3     # Number to use
cl <- makeCluster(ncore)
registerDoParallel(cl)

seeds <- 1:ncore


old.time <- Sys.time()
result <- foreach(x = seeds, .packages=c("nimble", "nimbleEcology")) %dopar% {
  set.seed(x)
  
  #### Create model ####
  IPM_mod <- nimbleModel(IPM_code,
                         constants = IPM_CONSTANTS,
                         data=IPM_DATA,
                         inits=append(IPM_INITS,
                                      list(mu.surv=c(runif(1, -1, 0), runif(1, 0.01, 1.5), runif(1, 1.51, 3)),
                                           p.const = runif(1, -1, 0.5), p.sd=runif(1, 0.1, 1),
                                           beta.LAT = rep(runif(2, -0.05, 0.05), times=c(2,1)),
                                           beta.OLD = rep(runif(2, -0.05, 0.05), times=c(2,1)),
                                           r.mean=runif(1, -0.5, 0.5)) ),
                         dimensions = list(phi = c(N.ind, (N.years)), Y.c=dim(IPM_DATA$Y.c))
                         # , check = T, calculate =T
                         , check = F, calculate =F
  )
  IPM_mod$initializeInfo()
  
  
  
  #### parameters to monitor ####
  parameters <- c("surv.exp", "mu.p", "p.sd", "beta.LAT", "beta.OLD", # CJS
                  'r.beta.lat', 'r.beta.col', 'r.mean', 'r.site.sd', 'r.proc.sd', 's.proc.sd') # Count
  parameters2 <- c("phi", "p", "N", "Y.c", "R", "S", 'r.site.mu') # latents for extra thinning
  
  
  
  #### MCMC settings and samplers ####
  conf_ipm <- configureMCMC(IPM_mod, monitors = parameters, monitors2 = parameters2,
                            autoBlock = F) # use blocking for production run
  conf_ipm # print settings
  
  mcmc <- buildMCMC(conf_ipm)
  
  
  #### Compilation ####
  Cmod <- compileNimble(IPM_mod, showCompilerOutput = T)
  Cmcmc <- compileNimble(mcmc, showCompilerOutput = T)
  
  
  
  #### Sample! ####
  samps <- runMCMC(Cmcmc, 
                   niter = ni, thin = nt
                   , thin2=nt2
                   , nburnin = nb, nchains = 1, 
                   samplesAsCodaMCMC = T)
  
  return(samps)
  
} # end result
Sys.time() - old.time


save(result, file="IPM_static_mini.RData")

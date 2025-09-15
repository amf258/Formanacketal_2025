################################################################################
### Bring in necessary packages
library(rjags)
library(jagsUI)

## Set working directory
setwd("/scratch/amf258/livedead/final_top_species")

################################################################################
### Read in argument

# Retrieve the second argument passed as an argument
args <- as.numeric(commandArgs(trailingOnly = TRUE))
k <- args[1]
status <-args[2]

################################################################################
### Read in files
source("./scripts/MCMCrestart_v2.r")

# Define a subset of species to analyze
spc<-c("abla", "auch", "piab", "potr", "abal",
       "pien", "pisy", "quru")
# spc<-c("auch", "potr", "quru")

if (status==1){ ## living
  ## Retrive argument
  k <- args[1]
  
  ## Data object for JAGS
  load(paste0("./data/L_", spc[k], "_datalist.Rdata"))
  
  ## Path to JAGS model script
  # model1 <-  "./scripts/Model_abal_v1.R"
  model2 <-  "./scripts/Model_abal_v3.R"
  
  ## Load inits
  load(file = paste0("./results/inits/L_", spc[k], "_inits4.Rdata"))
  
  ################################################################################
  ### Setting parameters and paths for JAGS run
  
  # update the model to track initial variables
  init.variables = c("delta","mu.mu.a0","mu.mu.aAge", "mu.mu.aR1",
                     "mu.a","sig.a","mu.aa","sig.aa","sig.mu.a0",
                     "sig.mu.aAge","sig.mu.aR1","mu.mu.AR1",
                     "sig.muAR1","sig","sig.AR1","sig.aR1",
                     "sig.aAge","sig.a0")
  # 
  # jm5 <- update(jm4, parameters.to.save=init.variables, n.adapt=NULL,
  #               n.iter=20)
  # 
  # 
  # newinits <-  initfind(jm5$samples)
  # 
  # inits1 = newinits$initials
  # 
  # variable.names = c("mu.a0", "mu.aR1", "mu.aAge", "mu.a",
  #                    "mu.aa", "sig.a0", "sig.aR1", "a", "aa",
  #                    "sig.a", "sig.aa", "sig.AR1", "sig.aAge",
  #                    "sig", "weight", "mu.mu.a0","mu.mu.aAge",
  #                    "sig.mu.aAge", "sig.mu.a0", "deviance")
  
  # ## Fit model 
  jm1 <- jagsUI::jags(data = dat.list,
                      inits = inits4,
                      model.file = model2,
                      parameters.to.save = init.variables,
                      # Rhat.limit=1.2,
                      # max.iter=200000,
                      # iter.increment=10000,
                      parallel = TRUE,
                      n.chains = 3,
                      n.iter = 2000,
                      n.burnin = 1000,
                      # n.thin = 2,
                      DIC = TRUE)
  
  
  # save(jm1, file = paste0("./results/mod1/L_", spc[k], "_mod1.Rdata"))
  # 
  newinits <-  initfind(jm1$samples)

  inits2 = newinits$initials
  
  variable.names = c("mu.a0", "mu.aR1", "mu.aAge", "mu.a",
                     "mu.aa", "sig.a0", "sig.aR1", "a", "aa",
                     "sig.a", "sig.aa", "sig.AR1", "sig.aAge",
                     "sig", "weight", "mu.mu.a0","mu.mu.aAge",
                     "sig.mu.aAge", "sig.mu.a0", "deviance")
  
  ## Fit model 
  jm5 <- jagsUI::autojags(data = dat.list,
                          inits = inits2,
                          model.file = model2,
                          parameters.to.save = variable.names,
                          Rhat.limit=1.2,
                          max.iter=50000,
                          iter.increment=8000,
                          parallel = TRUE,
                          n.chains = 3,
                          # n.iter = 5000,
                          n.burnin = 1000,
                          n.thin = 2,
                          DIC = TRUE)
  
  save(jm5, file = paste0("./results/mod5/L_", spc[k], "_mod5.Rdata"))
  
  # update variable names to include net sensitivities
  # variable.names = c("mu.a", "mu.aa", "sig.a", "sig.aa", 
  #                    "sig",  "deviance", "dmudWP", 
  #                    "dmudWT", "dmudSP", "dmudST")
  # 
  # net_sens <- update(jm4, parameters.to.save=variable.names, n.adapt=NULL, 
  #                    n.iter=5000)
  # 
  # save(net_sens, file = paste0("./results/cp/L_post_", spc[k], "_net_sens.Rdata"))
  
}

if (status==2){ ## Dead
  ## Retrive argument
  k <- args[1]
  
  ## Data object for JAGS
  load(paste0("./data/D_", spc[k], "_datalist.Rdata"))
  
  ## Path to JAGS model script
  # model1 <-  "./scripts/Model_abal_v1.R"
  model2 <-  "./scripts/Model_abal_v3.R"
  
  ## Load inits
  load(file = paste0("./results/inits/D_", spc[k], "_inits4.Rdata"))
  
  ################################################################################
  ### Setting parameters and paths for JAGS run
  
  # update the model to track initial variables
  init.variables = c("delta","mu.mu.a0","mu.mu.aAge", "mu.mu.aR1",
                     "mu.a","sig.a","mu.aa","sig.aa","sig.mu.a0",
                     "sig.mu.aAge","sig.mu.aR1","mu.mu.AR1",
                     "sig.muAR1","sig","sig.AR1","sig.aR1",
                     "sig.aAge","sig.a0")
  
  # jm5 <- update(jm4, parameters.to.save=init.variables, n.adapt=NULL,
  #               n.iter=20)
  # 
  # 
  # newinits <-  initfind(jm5$samples)
  # 
  # inits1 = newinits$initials
  
  # variable.names = c("mu.a0", "mu.aR1", "mu.aAge", "mu.a",
  #                    "mu.aa", "sig.a0", "sig.aR1", "a", "aa",
  #                    "sig.a", "sig.aa", "sig.AR1", "sig.aAge",
  #                    "sig", "weight", "mu.mu.a0","mu.mu.aAge",
  #                    "sig.mu.aAge", "sig.mu.a0", "deviance")
  
  ## Fit model 
  jm1 <- jagsUI::jags(data = dat.list,
                      inits = inits4,
                      model.file = model2,
                      parameters.to.save = init.variables,
                      # Rhat.limit=1.2,
                      # max.iter=200000,
                      # iter.increment=10000,
                      parallel = TRUE,
                      n.chains = 3,
                      n.iter = 2000,
                      n.burnin = 1000,
                      # n.thin = 2,
                      DIC = TRUE)

  
  # save(jm1, file = paste0("./results/mod1/D_", spc[k], "_mod1.Rdata"))
  # 
  newinits <-  initfind(jm1$samples)

  inits2 = newinits$initials
  
  variable.names = c("mu.a0", "mu.aR1", "mu.aAge", "mu.a",
                     "mu.aa", "sig.a0", "sig.aR1", "a", "aa",
                     "sig.a", "sig.aa", "sig.AR1", "sig.aAge",
                     "sig", "weight", "mu.mu.a0","mu.mu.aAge",
                     "sig.mu.aAge", "sig.mu.a0", "deviance")
  
  ## Fit model 
  jm5 <- jagsUI::autojags(data = dat.list,
                          inits = inits2,
                          model.file = model2,
                          parameters.to.save = variable.names,
                          Rhat.limit=1.2,
                          max.iter=50000,
                          iter.increment=8000,
                          parallel = TRUE,
                          n.chains = 3,
                          # n.iter = 5000,
                          n.burnin = 1000,
                          n.thin = 2,
                          DIC = TRUE)
  
  save(jm5, file = paste0("./results/mod5/D_", spc[k], "_mod5.Rdata"))
  
  modfit <- update(jmPE, parameters.to.save=c("mu.LogWidth", "LogWidth"), n.adapt=NULL,
                   n.iter=1000)
  
  # update variable names to include net sensitivities
  # variable.names = c("mu.a", "mu.aa", "sig.a", "sig.aa", 
  #                    "sig",  "deviance", "dmudWP", 
  #                    "dmudWT", "dmudSP", "dmudST")
  # 
  # net_sens <- update(jm4, parameters.to.save=variable.names, n.adapt=NULL, 
  #                    n.iter=5000)
  # 
  # save(net_sens, file = paste0("./results/cp/L_post_", spc[k], "_net_sens.Rdata"))
  
}

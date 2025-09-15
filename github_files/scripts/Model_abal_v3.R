## This model is a hierarchical SAM that estimates seasonal ----
## climate effects up to 4 years prior to ring formation 


# data{
  # Standardizing climate variables by site (added to data prep script)
  #   for(y in 1:Nyears.climate){
  #     for(s in 1:Nsites){
  #       ppt.season[y,s,1] = (ppt.winter[y,s]-mean(ppt.winter[,s]))/sd((ppt.winter[,s]))
  #       ppt.season[y,s,2] = (ppt.summer[y,s]-mean(ppt.summer[,s]))/sd((ppt.summer[,s]))
  #       temp.season[y,s,1] = (temp.winter[y,s]-mean(temp.winter[,s]))/sd((temp.winter[,s]))
  #       temp.season[y,s,2] = (temp.summer[y,s]-mean(temp.summer[,s]))/sd((temp.summer[,s]))
  #     }
  #   }
  # Use this to assign equal weights to delta to try to get better
  # initial values for model parameters
  # # for(u in 1:Nsites){
  # for(v in 1:Nvars){
  #   for(s in 1:Nseasons){ # season
  #     for(j in 1:Nlags){
  #       delta[v,s,j] <- 1
  #     }
  #   }
  # }
  # # }
  #   # gamma dist'n arguments for folded-Cauchy priors, with
  #   # degrees of freedom (v) and scale parameter (A):
  #   # v = 2
  #   # A = 0.1
  #   # arg1 <- v/2
  #   # arg2 <- v*A*A/2
# }

model{
  
  for(r in 1:Nring){
    # Likelihood for log ring-width data:			
    LogWidth[r] ~ dnorm(mu.LogWidth[r], tau)
    # Replicated data for evaluating model fit:
    LogWidth.rep[r] ~ dnorm(mu.LogWidth[r], tau)
    # Mean model for age-detrended, tree (or core) level data. With core-level coefficients.
    # antX[year,V,season] = antecedent climate variable for variable V = 1 for precip, V = 2 for temp.
    # season = 1 for "cool" ("winter") season, season = 2 for "warm" ("summer") season.
    mu.LogWidth[r] <- a0[CoreID[r]] + aAge[CoreID[r]]*Age[r] +
      main.effects[r] + inter.effects[r] + aR1[CoreID[r]]*AR1[r]

    
    # Main effects of antecedent climate
    for(v in 1:Nvars){
      for(s in 1:Nseasons){
        main.effects1[r,v,s] <- a[site[r],v,s]*antX[Year[r],site[r],v,s]
      }
      main.effects2[r,v] <- sum(main.effects1[r,v,])
    }
    main.effects[r] <- sum(main.effects2[r,])
    
    
    # 2-way antecedent climate interactions
    for(j in 1:Ninteractions){
      inter.effects1[r,j] <- aa[site[r],j]*antX[Year[r],site[r],v1[j],s1[j]]*antX[Year[r],site[r],v2[j],s2[j]]
    }
    inter.effects[r] <- sum(inter.effects1[r,])
    
    # For missing "prior" ring width indices (i.e., initial ring in series)
    AR1[r]~dnorm(mu.AR1[site[r],Year[r]], tau.AR1)
  }
  
  # # Net senstivities site level
  # for(y in y.start:y.end){
  #   for(k in 1:Nsites){
  #         dydWP[y,k] = a[k,1,1] + aa[k,1]*antX[y,k,1,2] + ...
  #   }
  # }
  
  # Net senstivities site level
  for(y in 1:Nyears){ 
    dmudWP[y] <- mu.a[1,1] + mu.aa[1]*mean(antX[y,,1,2]) + mu.aa[2]*mean(antX[y,,2,1]) + 
      mu.aa[3]*mean(antX[y,,2,2])
    
    dmudWT[y] <- mu.a[2,1] + mu.aa[2]*mean(antX[y,,1,1]) + mu.aa[4]*mean(antX[y,,1,2]) + 
      mu.aa[6]*mean(antX[y,,2,2]) 
    
    dmudSP[y] <- mu.a[1,2] + mu.aa[1]*mean(antX[y,,1,1]) + mu.aa[4]*mean(antX[y,,2,1]) + 
      mu.aa[5]*mean(antX[y,,2,2]) 
    
    dmudST[y] <- mu.a[2,2] + mu.aa[3]*mean(antX[y,,1,1]) + mu.aa[5]*mean(antX[y,,1,2]) + 
      mu.aa[6]*mean(antX[y,,2,1]) 
  }
  
  # Compute antecedent climate variables for climate variable v and time "block" into 
  # the past t (t = 1 is the current year):
  #vary these by status as well
  for(v in 1:Nvars){
    for(s in 1:Nseasons){ # season
      # for(u in 1:Nsites){
      for(j in 1:Nlags){
        # Assign a dirichlet(1,1,...,1) prior to the importance weights using the "delta-trick"
        delta[v,s,j] ~ dgamma(1,1)
        #delta[v,s,j] ~ dunif(1,3)
        weight[v,s,j] <- delta[v,s,j]/sumD1[v,s]
        
        for(i in year.start:Nyears.climate){ 
          for(u in 1:Nsites){
            # Here, v = 1 represents precipitation, and v = 2 is temperature.
            # and, s = 1 is winter, and s = 2 is summer.
            antX1[i,u,v,s,j] <- weight[v,s,j]*(equals(v,1)*ppt.season[(i)-j+1,u,s] + 
                                                 equals(v,2)*temp.season[(i)-j+1,u,s])
            
          }
        }
      }
      # # sum of the unnormalized weights:
      # sumD1[v,s] <- sum(delta[v,s,])
      
      # compute the antecedent climate variables by summing the weighted values 
      # across all annual lag periods:
      for(i in year.start:Nyears.climate){ 
        for(u in 1:Nsites){ 
          # This should be set-up/coded so that the first row (year) in antX 
          # corresponds to Year ID = 1 in the ring width data file.
          antX[i - year.start + 1,u,v,s] <- sum(antX1[i,u,v,s,])
        }
      }
      
      # sum of the unnormalized weights:
      sumD1[v,s] <- sum(delta[v,s,])
      
    }
  }
  
  
  
  # Assign hierarchical priors to the core-level effects in the mean model,
  # and assign relatively non-informative priors to population-level
  # parameters:
  # Add parameter expansion:
  for(r in 1:Ncores){
    # intercept, Age, and AR1 effect:
    a0[r] ~ dnorm(mu.a0[siteID[r]],tau.a0)
    aAge[r] ~ dnorm(mu.aAge[siteID[r]],tau.aAge)
    aR1[r] ~ dnorm(mu.aR1[siteID[r]],tau.aR1)
    
  }
  
  #### SITE-level parameters
  for(q in 1:Nsites){
    mu.a0[q] ~ dnorm(mu.mu.a0,tau.mu.a0)
    mu.aAge[q] ~ dnorm(mu.mu.aAge,tau.mu.aAge)
    mu.aR1[q] ~ dnorm(mu.mu.aR1,tau.mu.aR1)
    for(v in 1:Nvars){
      for(s in 1:Nseasons){
        a[q,v,s] ~ dnorm(mu.a[v,s],tau.a[v,s])
      }
    }
    for(j in 1:Ninteractions){
      aa[q,j] ~ dnorm(mu.aa[j],tau.aa[j])
    }
  }
  
  # Priors for population-level parameters (across all sites):
  mu.mu.a0 ~ dnorm(0,0.0001)
  mu.mu.aAge ~ dnorm(0,0.0001)
  mu.mu.aR1 ~ dnorm(0,0.0001)
  for(v in 1:Nvars){
    for(s in 1:Nseasons){
      mu.a[v,s] ~ dnorm(0,0.0001)
      tau.a[v,s] <- pow(sig.a[v,s],-2)
      sig.a[v,s] ~ dunif(0,100)
    }
  }
  for(j in 1:Ninteractions){
    mu.aa[j] ~ dnorm(0,0.0001)
    tau.aa[j] <- pow(sig.aa[j],-2)
    sig.aa[j] ~ dunif(0,100)
  }
  
  sig.mu.a0 ~ dunif(0,100)
  sig.mu.aAge ~ dunif(0,100)
  sig.mu.aR1 ~ dunif(0,100)  
  
  tau.mu.a0 <- pow(sig.mu.a0,-2)# same as prior for tau.a0
  tau.mu.aAge <- pow(sig.mu.aAge,-2)
  tau.mu.aR1 <- pow(sig.mu.aR1,-2)
  
  
  
  # Priors for the mean associated with the AR1 likelihood (for potential missing data):
  for (b in 1:Nsites){
    mu.mu.AR1[b] ~ dunif(-2,2)
    for(y in 1:Nyears){
      mu.AR1[b,y] ~ dnorm(mu.mu.AR1[b], tau.muAR1)
    }
  }
  sig.muAR1 ~ dunif(0,100)
  tau.muAR1 <- pow(sig.muAR1,-2)
  
  # KO: Prior for overall mean; may want to change the lower/upper limits on the uniform
  # prior to better reflect the possible range of values for the AR1 variable.
  
  
  # Priors for standard deviations / precisions associated with the likelihoods:
  # KO: Will need to provide inits for sig.AR1 (use similar to what is use for sig)
  sig ~ dunif(0,100)
  tau <- pow(sig,-2)	
  sig.AR1 ~ dunif(0,100)
  sig.aR1 ~ dunif(0,100)
  sig.aAge ~ dunif(0,100)
  sig.a0 ~ dunif(0,100)
  tau.AR1 <- pow(sig.AR1,-2)	
  tau.aR1 <- pow(sig.aR1,-2)	
  tau.aAge <- pow(sig.aAge,-2)	
  tau.a0 <- pow(sig.a0,-2)	
}





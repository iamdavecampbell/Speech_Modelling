############################
# functions in this file:
# RUN_PT_SS           = run vanilla MCMC or use parallel tempering to sample from the posterior of a STATESPACE model.  This is the main function
# logprior_SS         + evaluate the log prior for the parameters
# model_propagate_SS  + one step ahead prediction for the model (constructed using matrix notation) top produce the mean, and sampled observation
# model_predict_SS    + predict ahead looping over all time steps 
# loglikelihood_SS    + calculate the COMPLETE DATA log likelihood for the model
# logpost_SS          + calculate the log posterior
# make_fake_data_SS   + generate a fake dataset with missing values etc
## dataload           + Load the dataset and also keep the "alphas" indicating an observation has being made
############################


## Run Parallel Tempering for the State Space Model Formulation
Run_PT_SS = function(niter,  #MCMC iterations
                  temperatures=1, # temperatures for PT
                  data, # raw data with potential missing values
                  P = 1, # max lag
                  nstops = 10,
                  Goal_acceptance_rate = c(A = .23,
                                           B = .23,
                                           C = .23,
                                           sigma2_delta = .44,
                                           sigma2_e = .44,
                                           pa = .44,
                                           X = .23), # target sampling acceptance rates
                  filename = NULL,
                  step_var = rep(0.1,7),# initial guess at transition variance can be user specified.
                  CCor = diag(1, ncol(data)*(ncol(data)-1)),# potential correlation structure for C if known 
                  sigma2_pars = c(delta = .01, e = .0001) # prior parameters for the sigmas.
          ){
# Perform parallel tempering or vanilla MCMC to obtain samples from the model:
#
# DTP_{jt} = \alpha_{jt} X_{jt} + e_{jt}
# X_{jt} = a_j+b_{j} * X_{jt-1} + \sum_{k\in \{1,...,J\} \setminus \{j\}}c_{jk} * X_{kt-1}+ \delta_{jt}$$
#
# INPUTS
# niter = number of MCMC iterations
# temperatures = used for running parallel tempering.  If you want vanilla MCMC, then set top 1.  For more chains use something like c(.25,.5,.75,.875,.95,1).    0 = sample from the prior, 1 = sample from the posterior.  
# data = just the raw data.  No extra columns, no dates.  Missing values are fine.  They are sampled.
# P = the maximum autoregressive lag (it probably works but is untested for anything bigger than 1)
# nstops = number of times to pause the MCMC and calculate the acceptance rate.  If the iterations are still in the first half, then the transition variance will be adjusted
#         consequently YOU MUST DISCARD THE FIRST HALF OF THE ITERATIONS since they probably won't have the correct limiting distribution
#
# Goal_acceptance_rate = target Metropolis Hastings acceptance rate for the sampling of {A, B, C, sigma2, DTP0}.  Typically want .44 for 1 dimension, shrinking down to .23 for 5 or more dimensions
#
# filename = for saving interim results NOTE: IF THE FILENAME BEING WITH "LOG" THEN IT IS ASSUMED THAT THE LOG DATASET IS BEING USED
#
# step_var = the Metropolis Hastings transition variance initial guess
# CCor = the covariance structure for the metropolis proposals for the C block.
#
#
#
# Outputs a list with:
#
# swappers = iteration by iteration summary of which chains were proposed to swap and the outcome.  For Parallel Tempering
# swap_success_rate = an upper triangle matrix showing the success rate for proposed swaps between two different temperature values
# accepts      = Counts of the number of times a proposed accepted MH move was accepted
# rate         = Metropolis Hastings acceptance rate ( accepted MH moves / proposed MH moves )
# step_var     = final transition variance for proposing MH moves 
# CCor         = final estimated correlation structure for proposing moves within the C block
# Goal_acceptance_rate = same as input
# filename      = same as input
# temperatures   = same as input
# log_posterior = values of the log un normalized posterior at the end of each iteration
# theta         = MCMC samples of model parmameters {A,B,C, sigma2_delta, sigma2_e, pa}
# sigma2        = MCMC samples of the model variance term
# DTP_missing   = MCMC samples of the missing observations 
# theta_labels  = names of parameters for the values in theta
# A_index       = indices to extract just the A block from the theta matrix
# B_index       = indices to extract just the B block from the theta matrix
# C_index       = indices to extract just the C block from the theta matrix
# DTP0_index    = indices to extract just the DTP0 block from the theta matrix
#
#
  
  Nparblocks = 7 #{A,B,C,sigma2_delta, sigma2_e, p_a, X} 
  Ncountries = J = ncol(data) #Ncountries
  Number_of_a = J = ncol(data) # one intercept per country
  Number_of_b = J*P            # cols of B are lag.
  Number_of_c = (J-1)*J        # We don't go into lags larger than 1 for the other country effects
  Number_of_pa  = 1           # It's always one, but might expand later to make it time varying.
  Number_of_sigma2 = 2         # process, then observation
  Number_of_X0 = J             # Hold onto the initial states (the rest too...)
  Ndata = length(data)
  # construct name labels:
  bgrid = expand.grid(1:J,1:P)
  cgrid = expand.grid(1:J,1:J,1); cgrid =cgrid |> filter(Var1 !=Var2);
  theta_labels = c(#a
    paste0("A",1:J),
    #"b1lag1",...,"b1lagP",...,"bJlag1",..., "bJlagP"
    sort(apply(bgrid,1,function(x){paste0("B",x[1],"lag",x[2])})),
    #c
    sort(apply(cgrid,1,function(x){paste0("C",x[1],"Xc",x[2],"lag",x[3])})),
    c("sigma2_delta", "sigma2_e"),
    "pa",
    paste0("X0_",1:J)
  )
  A_index     = theta_labels |> grep(pattern = "^A")
  B_index     = theta_labels |> grep(pattern = "^B")
  C_index     = theta_labels |> grep(pattern = "^C")
  sigma2_delta_index = theta_labels |> grep(pattern = "sigma2_delta")
  sigma2_e_index     = theta_labels |> grep(pattern = "sigma2_e")
  pa_index    = theta_labels |> grep(pattern = "^pa")
  X0_index    = theta_labels |> grep(pattern = "^X0")

  alphas = data*0
  alphas[data !=0 | !is.na(data)] = 1
  alphas[data ==0 |  is.na(data)] = 0
  
  # prep and preallocate:
  # start with all positive values...
  thetastart = runif(Number_of_a+Number_of_b+Number_of_c+Number_of_pa+Number_of_sigma2+Number_of_X0, min = 0, max = .1)     
  theta                 = rep(list(rbind(thetastart,
                                         matrix(NA, 
                                                ncol = Number_of_a+Number_of_b+Number_of_c+Number_of_pa+Number_of_sigma2+Number_of_X0,
                                                nrow = niter-1,
                                                dimnames = list(NULL,theta_labels)))), length(temperatures))
  if(length(step_var)==1){step_var          = rep(step_var, Nparblocks)}
  if(!is.list(step_var)){step_var           = rep(list(step_var), length(temperatures))}
  if(!is.list(CCor)){CCor                   = rep(list(CCor),     length(temperatures))}
  X0Cor = rep(list(diag(1,Ncountries)),     length(temperatures)) # initial guess at the latent state correlation
  # set up a starting point, even if it's the raw data with imputed mean.
  Xstart = data
  Xstart[1,Xstart[1,]==0] = mean(Xstart[1,Xstart[1,]!=0])
  for(ro in 2:nrow(Xstart)){
    Xstart[ro,which(is.na(Xstart[ro,]))] = Xstart[ro-1,which(is.na(Xstart[ro,]))] 
    Xstart[ro,Xstart[ro,]==0]            = Xstart[ro-1,Xstart[ro,]==0] 
  }
  # build the X matrix by flattening the data matrix, reconstruct the X matrix using:
  # matrix(Xmat[[1]][1,], ncol = J):
  Xmat            = rep(list(rbind(c(Xstart),
                            matrix(NA, nrow = niter-1, ncol = Ndata))),              length(temperatures)) 
  
  log_posterior   = rep(list(rep(NA,niter)),                                       length(temperatures))
  accepts         = rep(list(matrix(0, nrow = nstops, ncol = Nparblocks)),         length(temperatures))   # track acceptance count for each stop and for the components {A,B,C,sigma2}
  rate            = rep(list(matrix(NA,nrow = nstops, ncol = Nparblocks)),         length(temperatures))   # track acceptance rate  for each stop and for the components {A,B,C,sigma2}  
  
  current_stop = 1
  data_full = list()
  
  # goals for the acceptance rate:
  Ubound = Goal_acceptance_rate + .05
  Lbound = Goal_acceptance_rate - .05
  
  # temp object for the correlation adaptation
  CCor2 = list()
  X0Cor2 = list()
  # keep track of swap acceptance rates:
  swappers      = matrix(NA, ncol = 3, nrow = niter)
  # set.seed(1234)
  for(iter in 1:(niter-1)){
    # iter=0
    # iter = iter + 1  
    ##### Check on swap #####
    maybeswap = sample(0:(length(temperatures)+1),2,replace=TRUE)
    if(all(maybeswap[1]!=maybeswap[2],min(maybeswap)>0,max(maybeswap)<=length(temperatures),iter>2)){
      #propose a parameter swap
      post1pars1 = logpost_SS(theta = theta[[maybeswap[1]]][iter,], 
                              theta_labels,
                              P = P,
                              alphas,
                              data,
                              X = Xmat[[maybeswap[1]]][iter,]|> matrix(ncol = Ncountries),
                              temperature = temperatures[maybeswap[1]],
                              filename = filename,
                              sigma2_pars = c(delta = .01, e = .0001))
        
      post2pars2 = logpost_SS(theta = theta[[maybeswap[2]]][iter,], 
                              theta_labels,
                              P = P,
                              alphas,
                              data,
                              X = Xmat[[maybeswap[2]]][iter,]|> matrix(ncol = Ncountries),
                              temperature = temperatures[maybeswap[2]],
                              filename = filename,
                              sigma2_pars = c(delta = .01, e = .0001))
      
      post1pars2 = logpost_SS(theta = theta[[maybeswap[2]]][iter,], 
                              theta_labels,
                              P = P,
                              alphas,
                              data,
                              X = Xmat[[maybeswap[2]]][iter,]|> matrix(ncol = Ncountries),
                              temperature = temperatures[maybeswap[1]],
                              filename = filename,
                              sigma2_pars = c(delta = .01, e = .0001))
      post2pars1 = logpost_SS(theta = theta[[maybeswap[1]]][iter,], 
                              theta_labels,
                              P = P,
                              alphas,
                              data,
                              X = Xmat[[maybeswap[1]]][iter,]|> matrix(ncol = Ncountries),
                              temperature = temperatures[maybeswap[2]],
                              filename = filename,
                              sigma2_pars = c(delta = .01, e = .0001))
      
      
      if(runif(1)< exp(post2pars1+post1pars2   -post1pars1 - post2pars2)){
        #accept the swap
        theta[[maybeswap[1]]][iter+1,]   = theta[[maybeswap[2]]][iter,]
        theta[[maybeswap[2]]][iter+1,]   = theta[[maybeswap[1]]][iter,]
        
        Xmat[[maybeswap[1]]][iter+1,] = Xmat[[maybeswap[2]]][iter,]
        Xmat[[maybeswap[2]]][iter+1,] = Xmat[[maybeswap[1]]][iter,]
        
        log_posterior[[maybeswap[1]]][iter+1] = post1pars2
        log_posterior[[maybeswap[2]]][iter+1] = post2pars1
        
        swappers[iter,]=c(sort(maybeswap),1)
        
        remaining_chains = 1:length(temperatures)
        remaining_chains = remaining_chains[-maybeswap]

      }else{
        # no swap accepted
        theta[[maybeswap[1]]][iter+1,]   = theta[[maybeswap[1]]][iter,]
        theta[[maybeswap[2]]][iter+1,]   = theta[[maybeswap[2]]][iter,]
        
        Xmat[[maybeswap[1]]][iter+1,] = Xmat[[maybeswap[1]]][iter,]
        Xmat[[maybeswap[2]]][iter+1,] = Xmat[[maybeswap[2]]][iter,]
        
        log_posterior[[maybeswap[1]]][iter+1] = post1pars1
        log_posterior[[maybeswap[2]]][iter+1] = post2pars2
        
        swappers[iter,]  = c(sort(maybeswap),0)
        
        remaining_chains = 1:length(temperatures)
        remaining_chains = remaining_chains[-maybeswap]
      }
    }else{
      # no swap proposed:
      remaining_chains = 1:length(temperatures)
      swappers[iter,]  = 0
    }
    
    ###
    # Now go in and save all the sampled states (X)
    #
    #
    ###
    
    ## pick up from here.
    # - sample $p_\alpha$ since that's easy and only requires the presence/absense.  But later this will be turned into something time varying and therefore more insightful.
    # - **propose** $X_{j0}$ and then propagate $X_{jt}\mid X_{j0}$. Then make a decision for this MH step based on the collective set of $X_{j\cdot}$.
    # - **sample parameters in blocks** $a,b,c$ using MH.  Propose values, then propagate the model forward, assessing the data fit as you go.  
    # - **sample the variance terms** $\sigma^2_delta$ and $\sigma^2_e$ using Metropolis Hastings, though later this can be done using a more efficient and direct Gibbs step.
    #
    # sample A / decide
    # sample B / decide
    # sample C / decide
    # sample sigma2 / decide
    # sample p_a / decide
    #
    # decide on X using
    for(chain in remaining_chains){
      theta_use  = theta[[chain]][iter,]
      XUse       = Xmat[[chain]][iter,] |> matrix(ncol = Ncountries) # latent states
      
      
      ##### A #####
      # propose a value from an easy distribution
      Aprop               = rnorm(n = Number_of_a, mean = theta_use[A_index], sd = step_var[[chain]][1]);
      theta_prop          = theta_use
      theta_prop[A_index] = Aprop
      # the ratio of un-normalized posteriors.  Note that my proposal
      # distribution is symmetric so Q_{ij}=Q_{ji}
      
      
      log_post_prop = logpost_SS(theta = theta_prop, 
                            theta_labels = theta_labels,
                            data = data,
                            P = 1,
                            alphas = alphas,
                            X = XUse,
                            temperature = temperatures[chain],
                            filename = filename,
                            sigma2_pars = c(delta = .01, e = .0001))  # proposed
      log_post_old  = logpost_SS(theta = theta_use, 
                                 theta_labels = theta_labels,
                                 data = data,
                                 P = 1,
                                 alphas = alphas,
                                 X = XUse,
                                 temperature = temperatures[chain],
                                 filename = filename,
                                 sigma2_pars = c(delta = .01, e = .0001)) # last accepted
      logalpha =   log_post_prop - log_post_old
      
      if(!is.na(logalpha) && runif(1) < exp(logalpha)){
        accepts[[chain]][current_stop,1] = accepts[[chain]][current_stop,1]+1;
        theta_use     = theta_prop
        log_post_old  = log_post_prop
        
      }
      ##### B #####
      # propose a value from an easy distribution
      Bprop               = rnorm(n = Number_of_b, mean = theta_use[B_index], sd = step_var[[chain]][2]);
      theta_prop          = theta_use
      theta_prop[B_index] = Bprop
      # the ratio of un-normalized posteriors.  Note that my proposal
      # distribution is symmetric so Q_{ij}=Q_{ji}
      log_post_prop = logpost_SS(theta = theta_prop, 
                                 theta_labels = theta_labels,
                                 data = data,
                                 P = 1,
                                 alphas = alphas,
                                 X = XUse,
                                 temperature = temperatures[chain],
                                 filename = filename,
                                 sigma2_pars = c(delta = .01, e = .0001))  # proposed
      logalpha =   log_post_prop - log_post_old
      
      if(!is.na(logalpha) && runif(1) < exp(logalpha)){
        accepts[[chain]][current_stop,2] = accepts[[chain]][current_stop,2]+1;
        theta_use     = theta_prop
        log_post_old  = log_post_prop
      }
      ##### C #####
      # propose a value from an easy distribution
      # Cprop               = rnorm(n = Number_of_c, mean = theta_use[C_index], sd = step_var[[chain]][3]);
      Cprop               = rmvnorm(1, mean = theta_use[C_index], 
                                    sigma = diag(step_var[[chain]][3], length(C_index))%*%
                                      CCor[[chain]] %*% 
                                      diag(step_var[[chain]][3], length(C_index)))
      
      theta_prop          = theta_use
      theta_prop[C_index] = Cprop
      # the ratio of un-normalized posteriors.  Note that my proposal
      # distribution is symmetric so Q_{ij}=Q_{ji}
      log_post_prop = logpost_SS(theta = theta_prop, 
                                 theta_labels = theta_labels,
                                 data = data,
                                 P = 1,
                                 alphas = alphas,
                                 X = XUse,
                                 temperature = temperatures[chain],
                                 filename = filename,
                                 sigma2_pars = c(delta = .01, e = .0001))  # proposed
      logalpha =   log_post_prop - log_post_old
      
      if(!is.na(logalpha) && runif(1) < exp(logalpha)){
        accepts[[chain]][current_stop,3] = accepts[[chain]][current_stop,3]+1;
        theta_use     = theta_prop
        log_post_old  = log_post_prop
      }
      
      ##### X0 #####
      # This will be super slow to burn in since the X values can't move very far 
      # having all the data move far as well.  Hopefully PT will help, realistically sequential monte carlo will be better.
      # propose a value from an easy distribution
      X0prop               = rmvnorm(1, mean = theta_use[X0_index], 
                                       sigma = diag(step_var[[chain]][7], length(X0_index))%*%
                                         X0Cor[[chain]] %*% 
                                         diag(step_var[[chain]][7], length(X0_index)))
        
      theta_prop = theta_use
      theta_prop[X0_index] = X0prop
      # the ratio of un-normalized posteriors.  Note that my proposal
      # distribution is symmetric so Q_{ij}=Q_{ji}
      log_post_prop = logpost_SS(theta = theta_prop, 
                                 theta_labels = theta_labels,
                                 data = data,
                                 P = 1,
                                 alphas = alphas,
                                 X = XUse,
                                 temperature = temperatures[chain],
                                 filename = filename,
                                 sigma2_pars = c(delta = .01, e = .0001))  # proposed
        logalpha =   log_post_prop - log_post_old 
        
        if(!is.na(logalpha) && runif(1) < exp(logalpha)){
          accepts[[chain]][current_stop,7] = accepts[[chain]][current_stop,7]+1;
          theta_use     = theta_prop
          log_post_old  = log_post_prop
        }  
      
        ##### X* #####
        # This will be super duper slow since we now have a for loop that moves around every since piece of X
        # however at least it will allow things to progress.
        # using the model to rebuild the full trajectory will have terrible sampling properties.
        
        for(xindex in 2:nrow(XUse)){
           Xstarprop = XUse
           Xstarprop[xindex,]   = rmvnorm(1, mean = XUse[xindex,], 
                                          sigma = diag(step_var[[chain]][7], length(X0_index))%*%
                                            X0Cor[[chain]] %*% 
                                            diag(step_var[[chain]][7], length(X0_index)))
           # the ratio of un-normalized posteriors.  Note that my proposal
           # distribution is symmetric so Q_{ij}=Q_{ji}
           log_post_prop = logpost_SS(theta = theta_prop, 
                                      theta_labels = theta_labels,
                                      data = data,
                                      P = 1,
                                      alphas = alphas,
                                      X = Xstarprop,
                                      temperature = temperatures[chain],
                                      filename = filename,
                                      sigma2_pars = c(delta = .01, e = .0001))  # proposed
           logalpha =   log_post_prop - log_post_old 
           
           if(!is.na(logalpha) && runif(1) < exp(logalpha)){
             # not tracking acceptance rates here.  It's a lot to carry around.
             # accepts[[chain]][current_stop,7] = accepts[[chain]][current_stop,7]+1;
             XUse          = Xstarprop
             log_post_old  = log_post_prop
           }  
        }
        
      ##### pa #####
      # propose a value from an easy distribution
      paprop               = rtruncnorm(n = Number_of_pa, 
                                        a=0, 
                                        b=1, 
                                        mean = theta_use[pa_index], 
                                        sd = step_var[[chain]][6])
      theta_prop           = theta_use
      theta_prop[pa_index] = paprop
      # the ratio of un-normalized posteriors.  Note that my proposal
      # distribution is NOT symmetric
      
      
      log_post_prop = logpost_SS(theta = theta_prop, 
                                 theta_labels = theta_labels,
                                 data = data,
                                 P = 1,
                                 alphas = alphas,
                                 X = XUse,
                                 temperature = temperatures[chain],
                                 filename = filename,
                                 sigma2_pars = c(delta = .01, e = .0001))  # proposed
      logalpha =   log_post_prop - 
        log_post_old + 
        sum(log(dtruncnorm(theta_use[pa_index],  
                           a=0, b=1, 
                           mean = paprop, 
                           sd = step_var[[chain]][6]))) -
        sum(log(dtruncnorm(paprop, 
                           a=0, b=1, 
                           mean = theta_use[pa_index],  
                           sd = step_var[[chain]][6])))
      
      
      if(!is.na(logalpha) && runif(1) < exp(logalpha)){
        accepts[[chain]][current_stop,6] = accepts[[chain]][current_stop,6]+1;
        theta_use     = theta_prop
        log_post_old  = log_post_prop
        
      }
    
      
      
      
     
     
      ##### Sigma2_delta #####
      #
      # This can probably be done directly from the conditional distribution, 
      # This is not wrong, just inefficient
      #
      # propose a value from an easy distribution: truncated Normal.  
      sigma2_delta_prop          = rtruncnorm(n = 1, a=0, b=Inf, 
                                        mean = theta_use[sigma2_delta_index],
                                        sd = step_var[[chain]][4])
      theta_prop                     = theta_use
      theta_prop[sigma2_delta_index] = sigma2_delta_prop
      
      # the ratio of un-normalized posteriors.  Note that my proposal
      # distribution is NOT symmetric 
      
      log_post_prop = logpost_SS(theta = theta_prop, 
                                 theta_labels = theta_labels,
                                 data = data,
                                 P = 1,
                                 alphas = alphas,
                                 X = XUse,
                                 temperature = temperatures[chain],
                                 filename = filename,
                                 sigma2_pars = c(delta = .01, e = .0001))  # proposed
      
      logalpha =   log_post_prop - log_post_old + 
              log(dtruncnorm(theta_use[sigma2_delta_index],  
                                 a=0, b=Inf, 
                                 mean = sigma2_delta_prop, 
                                 sd = step_var[[chain]][4])) -
              log(dtruncnorm(sigma2_delta_prop, 
                                   a=0, b=Inf, 
                                   mean = theta_use[sigma2_delta_index],  
                                   sd = step_var[[chain]][4]))
              
      if(!is.na(logalpha) && runif(1) < exp(logalpha)){
        accepts[[chain]][current_stop,4] = accepts[[chain]][current_stop,4]+1;
        theta_use               = theta_prop
        log_post_old            = log_post_prop
      }
      
      
      
      ##### Sigma2_e #####
      #
      # This can probably be done directly from the conditional distribution, 
      # This is not wrong, just inefficient
      #
      # propose a value from an easy distribution: truncated Normal.  
      sigma2_e_prop          = rtruncnorm(n = 1, a=0, b=Inf, 
                                              mean = theta_use[sigma2_e_index],
                                              sd = step_var[[chain]][5])
      theta_prop                     = theta_use
      theta_prop[sigma2_e_index] = sigma2_e_prop
      
      # the ratio of un-normalized posteriors.  Note that my proposal
      # distribution is NOT symmetric 
      
      log_post_prop = logpost_SS(theta = theta_prop, 
                                 theta_labels = theta_labels,
                                 data = data,
                                 P = 1,
                                 alphas = alphas,
                                 X = XUse,
                                 temperature = temperatures[chain],
                                 filename = filename,
                                 sigma2_pars = c(delta = .01, e = .0001))  # proposed
      
      logalpha =   log_post_prop - log_post_old + 
        log(dtruncnorm(theta_use[sigma2_e_index],  
                       a=0, b=Inf, 
                       mean = sigma2_e_prop, 
                       sd = step_var[[chain]][5])) -
        log(dtruncnorm(sigma2_e_prop, 
                       a=0, b=Inf, 
                       mean = theta_use[sigma2_e_index],  
                       sd = step_var[[chain]][5]))
      
      if(!is.na(logalpha) && runif(1) < exp(logalpha)){
        accepts[[chain]][current_stop,5] = accepts[[chain]][current_stop,5]+1;
        theta_use               = theta_prop
        log_post_old            = log_post_prop
      }
      
      
      
      
      
      
      ##### Lock in the Output ##### 
      theta[[chain]][iter+1,]        = theta_use
      Xmat[[chain]][iter+1,]         = XUse
      log_posterior[[chain]][iter+1] = log_post_old
      
      }# end loop over NChain
    
    ##### Track sampling performance indicators  #####
    if(iter==floor(current_stop*niter/nstops)){  
      # calculate the acceptance rate.  Note that I use 1+ accepts / 2+ iterations.
      # If I don't that and we don't accept any proposals the MCMC will crash and fail.
      
      for(chain in 1:length(temperatures)){
        denom = niter/nstops - # potential MH swaps within chain - potential chain to chain swaps
          length(which(
            swappers[((current_stop-1)*niter/nstops+1):(((current_stop)*niter)/nstops),1:2] == chain
          ))
        # cat(paste("iter", iter, " chain ", chain, " denom ", denom, " acc ", past  e(acce  pts[[chain]][current_stop,], collapse = "_"), "\n"))
        rate[[chain]][current_stop,]= (1+accepts[[chain]][current_stop,])/(2+denom);
        # cat(paste("iter = ",iter," chain = ", chain, " Acceptance rate = ", past  e(roun  d(rate[[chain]][current_stop,], 4), collapse = " , " ), "\n"))
      }
      cat(paste("iter:", iter, " rate ", paste(round(rate[[chain]][current_stop-1,],4), collapse = "_"),"\n"))
      current_stop = current_stop+1;
      
      if(!is.null(filename)){save.image(filename)} # optional, but helpful}
      for(chain in 1:length(temperatures)){
        for(par_index in 1:Nparblocks){
          if((rate[[chain]][current_stop-1,par_index] > Ubound[par_index] ||  
              rate[[chain]][current_stop-1,par_index] < Lbound[par_index])  && 
             iter<= niter/2){
            # do the adaptation but stop by the halfway point
            step_var[[chain]][par_index] = step_var[[chain]][par_index] * rate[[chain]][current_stop-1,par_index] / Goal_acceptance_rate[par_index]; 
          }
        }
        # Cmat correlation structure:
        # modify if it has been run for a while, but stop before it's been running for too long.  & Also make sure that there have been some accepted values:
        if( current_stop>2 &&
            iter<= niter/2 && 
            nrow( unique(theta[[chain]][floor(iter/2):iter,C_index]))>2){ # less often make an adjustment to the correlation structure: 
          CCor2[[chain]] = cor(theta[[chain]][floor(iter/2):iter,C_index])*.8 # shrink correlations down since I'm sure they're wrong
          diag(CCor2[[chain]]) = 1
          CCor[[chain]] = (CCor[[chain]] + CCor2[[chain]]) / 2 # average them so that it moves in a less twitchy manner.
        }
        # X0 correlation structure:
        if( current_stop>2 &&
            iter<= niter/2 && 
            nrow( unique(theta[[chain]][floor(iter/2):iter,X0_index]))>2){ # less often make an adjustment to the correlation structure: 
          X0Cor2[[chain]] = cor(theta[[chain]][floor(iter/2):iter,X0_index])*.8 # shrink correlations down since I'm sure they're wrong
          diag(X0Cor2[[chain]]) = 1
          X0Cor[[chain]] = (X0Cor[[chain]] + X0Cor2[[chain]]) / 2 # average them so that it moves in a less twitchy manner.
        }
      }
    }
    
  }# end loop over iterations
  
  
  ACCTable = matrix(0, nrow = length(temperatures), ncol = length(temperatures))
  PosTable = matrix(0, nrow = length(temperatures), ncol = length(temperatures))
  if(length(temperatures)>1){
    SWAPS = swappers[rowSums(swappers)!=0,]
    for(ro in 1:(nrow(SWAPS))){
      PosTable[SWAPS[ro,1],SWAPS[ro,2]] =  ACCTable[SWAPS[ro,1],SWAPS[ro,2]] + 1  
      ACCTable[SWAPS[ro,1],SWAPS[ro,2]] =  ACCTable[SWAPS[ro,1],SWAPS[ro,2]] + SWAPS[ro,3]
    }
    swap_success_rate = ACCTable/PosTable
  }else{swap_success_rate = SWAPS = NULL}
  
  return(list(
    swappers             = swappers,
    swap_success_rate    = swap_success_rate,
    accepts              = accepts,
    rate                 = rate,
    step_var             = step_var,
    CCor                 = CCor,
    Goal_acceptance_rate = Goal_acceptance_rate,
    filename             = filename,
    temperatures         = temperatures,
    log_posterior        = log_posterior,
    theta                = theta,
    theta_labels         = theta_labels,
    A_index     = A_index,
    B_index     = B_index,
    C_index     = C_index,
    sigma2_delta_index = sigma2_delta_index,
    sigma2_e_index     = sigma2_e_index,
    pa_index    = pa_index,
    X0_index    = X0_index,
    Xmat = Xmat,
    alphas = alphas
  ))
}


## dataload
dataload = function(datafile, topic_of_interest, countries2use, start_date = ymd("2008-01-01")){
  
  dataset = read_csv(datafile)
  data = dataset |> 
    filter(anchor == topic_of_interest)|>
    pivot_wider( names_from = country, 
                 values_from = avg_proportion) |>
    arrange(date)|>
    select(contains(c(countries2use, "date", "anchor")))|> 
    filter(date >= start_date)
  
  dates_in_use =   data |> pull(date)
  
  data =   data |> select(-c("date", "anchor")) |> as.matrix()
  
  alphas = matrix(0, nrow = nrow(data), ncol = ncol(data))
  alphas[!is.na(data)] = 1
  
  return(list(dates_in_use = dates_in_use,
              data   =   data,
              alphas = alphas))
  
}








## logprior_SS
# theta,
# theta_labels,
# alphas,
# filename
# sigma2_pars
logprior_SS = function(theta,theta_labels, alphas,filename = NULL, sigma2_pars = c(delta = .01, e = .0001)){
  # note the order of the sigma2 prior parameters must match that of theta
  # $$\sigma_{\delta}^2\sim exponential(.01)$$ state process noise (converts to prior mean for the SD of .1)
  # $$\sigma_{e}^2\sim exponential(.0001)$$     observation noise (converts to prior mean for the SD of .01)
  A_index     = theta_labels |> grep(pattern = "^A")
  B_index     = theta_labels |> grep(pattern = "^B")
  C_index     = theta_labels |> grep(pattern = "^C")
  sigma_index = theta_labels |> grep(pattern = "^sigma")
  pa_index    = theta_labels |> grep(pattern = "^pa")
  X0_index    = theta_labels |> grep(pattern = "^X0")
  
  if(any(theta[sigma_index] <= 0)){ #case 1 negative variance
    # negative variance is bad
    logprior = -Inf
  }else{#positive variance
      logprior = sum(dnorm(theta[c(A_index,B_index,C_index,X0_index)], log = TRUE)) +
        # theta[pa_index] # uniform prior
        -sum(theta[sigma_index]/sigma2_pars)
    }
  return(logprior)
}

## model_propagate_SS
# ONE STEP AHEAD ONLY. X Only
model_propagate_SS = function(theta,theta_labels, Xstep_minus_P,P = 1){
  # Xstep_minus_P = matrix(X_full[(time_point-P):(time_point-1),], nrow = P)
  # P = lag depth
  # predict the model ahead by one time increment
  #$$X_{jt} = a_j+b_{j} * X_{jt-1} + \sum_{k\in \{1,...,J\} \setminus \{j\}}c_{jk} * X_{kt-1}+ \delta_{jt}$$
  J     = ncol(Xstep_minus_P) #Ncountries
  A_index     = theta_labels |> grep(pattern = "^A")
  B_index     = theta_labels |> grep(pattern = "^B")
  C_index     = theta_labels |> grep(pattern = "^C")
  sigma2_delta_index = theta_labels |> grep(pattern = "sigma2_delta")
  pa_index    = theta_labels |> grep(pattern = "^pa")
  X0_index    = theta_labels |> grep(pattern = "^X0")
  
  a     = theta[A_index] |> matrix(nrow = P)
  # cols of B are lag.  Rows of b are country
  b     = theta[B_index] |> matrix(nrow = P, byrow = TRUE)
  
  ctemp = theta[C_index]|> matrix(nrow = J, byrow = TRUE) 
  # for the sake of matrix multiplication further down, make this a square with diagonal zeros.
  cmat = matrix(NA, ncol = J, nrow = J)
  # Fill upper triangle from  (excluding diagonal)
  cmat[upper.tri(cmat)] = ctemp[upper.tri(ctemp, diag = TRUE)]
  
  # Fill lower triangle from M (excluding diagonal)
  cmat[lower.tri(cmat)] = ctemp[lower.tri(ctemp, diag = FALSE)]
  
  # Fill diagonal with D
  diag(cmat) = 0
  
  #$$X_{jt} = a_j+b_{j} * X_{jt-1} + \sum_{k\in \{1,...,J\} \setminus \{j\}}c_{jk} * X_{kt-1}+ \delta_{jt}$$
  # so that we take in the most recently evaluated Xstep_minus_P and put the lagged version in the right place
  error = rnorm(J,mean = 0,sd = sqrt(theta[sigma2_delta_index])) |> matrix(nrow = P)
  predicted_mean = a + 
    apply(b*Xstep_minus_P,2,sum) + # sum down columns if P>1
    t(cmat %*% Xstep_minus_P[nrow(Xstep_minus_P),])
  return(list(predicted_mean = predicted_mean,
              error = error,
              prediction = predicted_mean+error))
}

## model_predict
model_predict_SS = function(theta,theta_labels,data, P = 1){
  
  A_index     = theta_labels |> grep(pattern = "^A")
  B_index     = theta_labels |> grep(pattern = "^B")
  C_index     = theta_labels |> grep(pattern = "^C")
  sigma_index = theta_labels |> grep(pattern = "^sigma")
  pa_index    = theta_labels |> grep(pattern = "^pa")
  X0_index    = theta_labels |> grep(pattern = "^X0")
  
  X_full = data_full = data # just to preallocate; these are replaced below.
  # Using the past, predict one step ahead to infill the present. 
  X_full[1,] = theta[X0_index]
  alphas = matrix(sample(x=c(0,1), size = length(data), replace = TRUE, prob = c(theta[pa_index],1-theta[pa_index])),
                  nrow = nrow(data), ncol = ncol(data))

  for(time_point in (P+1):nrow(data_full)){
    modepred                = model_propagate_SS(theta,
                                                 theta_labels,
                                                 matrix(X_full[(time_point-P):(time_point-1),], nrow = P), 
                                                 P = P)
    X_full[time_point,]     = modepred$prediction
  }
  # convert latent states to observations:
  data_full = alphas * (X_full + matrix(rnorm(n = length(data), mean = 0, sd = sqrt(theta["sigma2_e"])),
                                       nrow = nrow(data), ncol = ncol(data)))
  return(list(X_full     = X_full,
              data_full  = data_full,
              alphas = alphas))
}




## loglikelihood_SS: Complete data likelihood.
# theta       
# theta_labels
# alphas
# data  
# X     
# P   
# Log likelihood for the State Space is calculated for the observed data conditional on the infilled missing data.
loglikelihood_SS = function(theta,
                            theta_labels,
                            alphas,
                            data, 
                            X, 
                            P = 1){
  A_index        = theta_labels |> grep(pattern = "^A")
  B_index        = theta_labels |> grep(pattern = "^B")
  C_index        = theta_labels |> grep(pattern = "^C")
  sigma2_delta_index    = theta_labels |> grep(pattern = "sigma2_delta")
  sigma2_e_index = theta_labels |> grep(pattern = "sigma2_e")
  pa_index       = theta_labels |> grep(pattern = "^pa")
  X0_index       = theta_labels |> grep(pattern = "^X0")
  
  loglik =  
    #P(alpha | p_\alpha)
    sum(alphas==1) * log(theta[pa_index]) + sum(alphas==0) * log(1-theta[pa_index]) 
    #P(X_t|X_{t-1}, theta)
    for(time_point in (P+1):nrow(data)){
      modepred                = model_propagate_SS(theta,theta_labels,
                                                   matrix(
                                                     X[(time_point-P):(time_point-1),], nrow = P), P = P)
      loglik = loglik + 
        sum(
          dnorm(mean = modepred$predicted_mean,
                sd   = sqrt(theta[sigma2_delta_index]), 
                X[time_point,],
                log = TRUE)
        )
    }
    # observations conditional on states (X)
  loglik =  loglik +
    sum(dnorm(mean = X, 
          sd   = sqrt(theta[sigma2_e_index]), 
          data, 
          log  = TRUE)[alphas ==1] # exclude missing values in the likelihood calculation since those are handled elsewhere.
      ) 
  return(loglik)
}






## logpost_SS
logpost_SS = function(theta, 
                      theta_labels,
                      P = 1,
                      alphas,
                      data,
                      X,
                      temperature = 1,
                      filename = NULL,
                      sigma2_pars = c(delta = .01, e = .0001)){
  temperature*sum(loglikelihood_SS(theta        = theta,
                                   theta_labels = theta_labels,
                                   alphas = alphas,
                                   data   = data,
                                   X      = X,
                                   P      = P
                                   )) +
    logprior_SS(theta,
                theta_labels,
                alphas,
                filename,
                sigma2_pars = sigma2_pars)
  
  
}







## make_fake_data
make_fake_data_SS = function(Nmissing = 25, Number_of_DTP0 = 1, Ncountries = 5, Ntimes = 200){
  
  data = matrix(1, ncol = Ncountries, nrow = Ntimes)   #preallocate
  data[1,]   = runif(Ncountries)*.25                   #actual starting point with one missing observation
  DTP0_missing_index = sample(size = Number_of_DTP0, 1:Ncountries)
  data[1,DTP0_missing_index]   = NA
  miss       = sample(
    setdiff(1:(Ntimes*Ncountries), 1+Ntimes*0:(Ncountries-1)), 
    Nmissing) # set up missing values (can't be starting point)
  data[miss] = NA             #place eventual missing values
  missing_index = is.na(data) # the most useful format is logicals for locations.
  # NOTE THAT DTP0 is handled differently so remove it from the missing index
  missing_index[1,] = FALSE
  data[-1,] = NA              # we'll simulate the whole dataset from here down below
  
  J      = Ncountries #Ncountries
  P      = 1 # max lag
  Number_of_a = J
  Number_of_b = J*P
  Number_of_c = (J-1)*J 
  # for now we don't have any lags beyond 1 within the country (C) matrix because I'd need to think about
  # how to store and handle them in the propagation model
  Number_of_pa  = 1           # It's always one, but might expand later to make it time varying.
  Number_of_sigma2 = 2         # process, then observation
  Number_of_X0 = J             # Hold onto the initial states (the rest too...)
  
  bgrid = expand.grid(1:J,1:P)
  cgrid = expand.grid(1:J,1:J,1); cgrid =cgrid |> filter(Var1 !=Var2);
  theta_labels = c(#a
    paste0("A",1:J),
    #"b1lag1",...,"b1lagP",...,"bJlag1",..., "bJlagP"
    sort(apply(bgrid,1,function(x){paste0("B",x[1],"lag",x[2])})),
    #c
    sort(apply(cgrid,1,function(x){paste0("C",x[1],"Xc",x[2],"lag",x[3])})),
    c("sigma2_delta", "sigma2_e"),
    "pa",
    paste0("X0_",1:J)
  )
  A_index     = theta_labels |> grep(pattern = "^A")
  B_index     = theta_labels |> grep(pattern = "^B")
  C_index     = theta_labels |> grep(pattern = "^C")
  sigma_index = theta_labels |> grep(pattern = "^sigma")
  pa_index    = theta_labels |> grep(pattern = "^pa")
  X0_index    = theta_labels |> grep(pattern = "^X0")
  
  # deal with starting time point missing values:
  X0   = runif(n = Number_of_X0)
  
  theta_true  =  c(rnorm((Number_of_a+Number_of_b+Number_of_c), sd = .1), 
                   runif(1, 0, .1 ),# sigma2_delta (process noise)
                   runif(1, 0, .01 ),# sigma2_e observational
                   runif(Number_of_pa,0,.1),
                   X0)
  names(theta_true) = theta_labels
  
  
  
  fake_data = model_predict_SS(theta_true,theta_labels,data, P = 1)

  return(list(
    theta_true = theta_true,
    X_full     = fake_data$X_full,  
    data_full  = fake_data$data_full,
    alphas     = fake_data$alphas
      )
  )
}





















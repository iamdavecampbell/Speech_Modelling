############################
# functions in this file:
# RUN_PT           = run vanilla MCMC or use parallel tempering to sample from the posterior.  This is the main function
# logprior         = evaluate the log prior for the parameters
# model_propagate  = one step ahead prediction for the model (constructed using matrix notation) top produce the mean, and sampled observation
# model_predict    = predict ahead looping over all time steps 
# infill missing   = reconstruct the data but with single imputed values.  These are MCMC samples anyways so it's fine.
# loglikelihood    = calculate the log likelihood for the model
# logpost          = calculate the log posterior
# make_fake_data   = generate a fake dataset with missing values etc
# dataload         = load the real dataset.
############################


Run_PT = function(niter,  #MCMC iterations
                  temperatures=1, # temperatures for PT
                  data, # raw data with potential missing values
                  P= 1, # max lag
                  nstops = 10,
                  Goal_acceptance_rate = c(.23,.23,.23,.44,.44), # target sampling acceptance rate for {A, B, C, sigma2, DTP0}
                  filename = NULL,
                  step_var = rep(0.1, 5),# initial guess at transition variance can be user specified.
                  CCor = diag(1, ncol(data)*(ncol(data)-1)) # potential correlation structure for C if known 
){
  
  ################
  # Perform parallel tempering or vanilla MCMC to obtain samples from the model:
  #
  # DTP_{ijt} = a_j + \sum_pb_{jp} * DTP_{ijt-p} + \sum_{k\in \{1,...,J\} \setminus \{j\}}c_{jk} * DTP_{ikt-1}  + e_{i,j,t}
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
  # theta         = MCMC samples of model parmameters {A,B,C, DTP0}
  # sigma2        = MCMC samples of the model variance term
  # DTP_missing   = MCMC samples of the missing observations 
  # theta_labels  = names of parameters for the values in theta
  # A_index       = indices to extract just the A block from the theta matrix
  # B_index       = indices to extract just the B block from the theta matrix
  # C_index       = indices to extract just the C block from the theta matrix
  # DTP0_index    = indices to extract just the DTP0 block from the theta matrix
  #
  ##### 
  
  missing_index      = is.na(data) # the most useful format is logicals for locations.
  DTP0_missing_index = which(is.na(data[1,]))
  Number_of_DTP0 = length(DTP0_missing_index)
  Nparblocks = if(Number_of_DTP0>0){5}else{4}#{A,B,C,sigma2, DTP0?}
  
  Ncountries = J = ncol(data) #Ncountries
  Number_of_a = J = ncol(data) # one intercept per country
  Number_of_b = J*P            # cols of B are lag.
  Number_of_c = (J-1)*J        # We don't go into lags larger than 1 for the other country effects
  
  # construct name labels:
  bgrid = expand.grid(1:J,1:P)
  cgrid = expand.grid(1:J,1:J,1); cgrid =cgrid |> filter(Var1 !=Var2);
  theta_labels = c(#a
    paste0("A",1:J),
    #"b1lag1",...,"b1lagP",...,"bJlag1",..., "bJlagP"
    sort(apply(bgrid,1,function(x){paste0("B",x[1],"lag",x[2])})),
    #c
    sort(apply(cgrid,1,function(x){paste0("C",x[1],"Xc",x[2],"lag",x[3])})),
    if(Number_of_DTP0>0){paste0("DTP0_",DTP0_missing_index)}else{NULL}
  )
  A_index = theta_labels |> grep(pattern = "^A")
  B_index = theta_labels |> grep(pattern = "^B")
  C_index = theta_labels |> grep(pattern = "^C")
  DTP0_index = theta_labels |> grep(pattern = "^DTP")
  
  
  
  
  # prep and preallocate:
  
  thetastart = c(runif((Number_of_a+Number_of_b+Number_of_c)),
                 if(Number_of_DTP0>0){runif(Number_of_DTP0)*.25}else{NULL})# trying to be smart about the initial value
  theta                 = rep(list(rbind(thetastart,
                                         matrix(NA, 
                                                ncol = Number_of_a+Number_of_b+Number_of_c+Number_of_DTP0, 
                                                nrow = niter-1,
                                                dimnames = list(NULL,theta_labels)))), length(temperatures))
  if(!is.list(step_var)){step_var           = rep(list(step_var), length(temperatures))}
  if(!is.list(CCor)){CCor                   = rep(list(CCor),     length(temperatures))}
  sigma2          = rep(list(c(1,rep(NA, niter-1))),                                      length(temperatures))
  DTP_missing     = rep(list(matrix(NA, ncol = sum(missing_index), nrow = niter)), length(temperatures))
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
  
  # keep track of swap acceptance rates:
  swappers      = matrix(NA, ncol = 3, nrow = niter)
  # set.seed(1234)
  for(iter in 1:(niter-1)){
    # iter=0
    # iter = iter + 1  
    # Check on swap:
    maybeswap = sample(0:(length(temperatures)+1),2,replace=TRUE)
    if(all(maybeswap[1]!=maybeswap[2],min(maybeswap)>0,max(maybeswap)<=length(temperatures),iter>2)){
      #propose a parameter swap
      post1pars1 = logpost(theta = theta[[maybeswap[1]]][iter,], 
                           sigma2 = sigma2[[maybeswap[1]]][iter], 
                           data = data, 
                           P = P,  missing_index,
                           temperature = temperatures[maybeswap[1]],
                           DTP0=theta[[maybeswap[1]]][iter,DTP0_index],
                           filename = filename)
      post2pars2 = logpost(theta = theta[[maybeswap[2]]][iter,], 
                           sigma2 = sigma2[[maybeswap[2]]][iter], 
                           data = data, 
                           P = P,  missing_index,
                           temperature = temperatures[maybeswap[2]],
                           DTP0=theta[[maybeswap[2]]][iter,DTP0_index],
                           filename = filename)
      
      post1pars2 = logpost(theta = theta[[maybeswap[2]]][iter,], 
                           sigma2 = sigma2[[maybeswap[2]]][iter], 
                           data = data, 
                           P = P,  missing_index,
                           temperature = temperatures[maybeswap[1]],
                           DTP0=theta[[maybeswap[2]]][iter,DTP0_index],
                           filename = filename)
      post2pars1 = logpost(theta = theta[[maybeswap[1]]][iter,], 
                           sigma2 = sigma2[[maybeswap[1]]][iter], 
                           data = data, 
                           P = P,  missing_index,
                           temperature = temperatures[maybeswap[2]],
                           DTP0=theta[[maybeswap[1]]][iter,DTP0_index],
                           filename = filename)
      
      
      if(runif(1)< exp(post2pars1+post1pars2   -post1pars1 - post2pars2)){
        #accept the swap
        theta[[maybeswap[1]]][iter+1,]   = theta[[maybeswap[2]]][iter,]
        theta[[maybeswap[2]]][iter+1,]   = theta[[maybeswap[1]]][iter,]
        sigma2[[maybeswap[1]]][iter+1]   = sigma2[[maybeswap[2]]][iter]
        sigma2[[maybeswap[2]]][iter+1]   = sigma2[[maybeswap[1]]][iter]
        DTP_missing[[maybeswap[1]]][iter+1,] = DTP_missing[[maybeswap[2]]][iter,]
        DTP_missing[[maybeswap[2]]][iter+1,] = DTP_missing[[maybeswap[1]]][iter,]
        log_posterior[[maybeswap[1]]][iter+1] = post1pars2
        log_posterior[[maybeswap[2]]][iter+1] = post2pars1
        swappers[iter,]=c(sort(maybeswap),1)
        remaining_chains = 1:length(temperatures)
        remaining_chains = remaining_chains[-maybeswap]
        
        dt1 = data_full[[maybeswap[1]]]
        dt2 = data_full[[maybeswap[2]]]
        data_full[[maybeswap[1]]] = dt2
        data_full[[maybeswap[1]]] = dt1
        
      }else{
        # no swap accepted
        theta[[maybeswap[1]]][iter+1,]   = theta[[maybeswap[1]]][iter,]
        theta[[maybeswap[2]]][iter+1,]   = theta[[maybeswap[2]]][iter,]
        sigma2[[maybeswap[1]]][iter+1]   = sigma2[[maybeswap[1]]][iter]
        sigma2[[maybeswap[2]]][iter+1]   = sigma2[[maybeswap[2]]][iter]
        log_posterior[[maybeswap[1]]][iter+1] = post1pars1
        log_posterior[[maybeswap[2]]][iter+1] = post2pars2
        DTP_missing[[maybeswap[1]]][iter+1,] = DTP_missing[[maybeswap[1]]][iter,]
        DTP_missing[[maybeswap[2]]][iter+1,] = DTP_missing[[maybeswap[2]]][iter,]
        swappers[iter,]  = c(sort(maybeswap),0)
        
        remaining_chains = 1:length(temperatures)
        remaining_chains = remaining_chains[-maybeswap]
      }
    }else{
      # no swap proposed:
      remaining_chains = 1:length(temperatures)
      swappers[iter,]  = 0
    }
    
    for(chain in remaining_chains){
      theta_use  = theta[[chain]][iter,]
      sigma2_use = sigma2[[chain]][iter]
      
      # start with the data
      data_full[[chain]] = data
      # fill in the sampled starting point
      data_full[[chain]][1,DTP0_missing_index] = theta_use[DTP0_index]
      # sample the missing values
      data_full[[chain]]          = 
        infill_missing(theta = theta_use,
                       sigma2 = sigma2_use,
                       data = data,
                       P = P)
      DTP_missing[[chain]][iter+1,] = data_full[[chain]][missing_index]
      
      
      
      ##### A #####
      # propose a value from an easy distribution
      Aprop               = rnorm(n = Number_of_a, mean = theta_use[A_index], sd = step_var[[chain]][1]);
      theta_prop          = theta_use
      theta_prop[A_index] = Aprop
      # the ratio of un-normalized posteriors.  Note that my proposal
      # distribution is symmetric so Q_{ij}=Q_{ji}
      log_post_prop = logpost(theta = theta_prop, 
                              sigma2 = sigma2_use, 
                              data = data, 
                              P = P,  missing_index,
                              DTP0 = theta_prop[DTP0_index],
                              temperature = temperatures[chain],
                              filename = filename)  # proposed
      log_post_old  = logpost(theta = theta_use,  
                              sigma2 = sigma2_use, 
                              data = data,
                              P = P,  missing_index,
                              DTP0 = theta_use[DTP0_index],
                              temperature = temperatures[chain],
                              filename = filename) # last accepted
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
      log_post_prop = logpost(theta = theta_prop, 
                              sigma2 = sigma2_use, 
                              data = data,
                              P = P,  
                              missing_index,
                              DTP0 = theta_prop[DTP0_index],
                              temperature = temperatures[chain],
                              filename = filename)  # proposed
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
      log_post_prop = logpost(theta = theta_prop, 
                              sigma2 = sigma2_use, 
                              data = data,
                              P = P,  
                              missing_index,
                              DTP0 = theta_prop[DTP0_index],
                              temperature = temperatures[chain],
                              filename = filename)  # proposed
      logalpha =   log_post_prop - log_post_old
      
      if(!is.na(logalpha) && runif(1) < exp(logalpha)){
        accepts[[chain]][current_stop,3] = accepts[[chain]][current_stop,3]+1;
        theta_use     = theta_prop
        log_post_old  = log_post_prop
      }
      
      ##### DTP0 #####. 
      # This will be super slow to burn in since DTP0 can't move very far without 
      # having all the data move far as well.  Probably PT will help a lot here.
      # propose a value from an easy distribution
      if(Number_of_DTP0>0){
        DTP0prop               = rtruncnorm(n = length(DTP0_index), 
                                            a=0, 
                                            b=1, 
                                            mean = theta_use[DTP0_index], 
                                            sd = step_var[[chain]][5])
        theta_prop = theta_use
        theta_prop[DTP0_index] = DTP0prop
        data_full_prop = data_full[[chain]]
        data_full_prop[1,DTP0_missing_index] = DTP0prop
        # the ratio of un-normalized posteriors.  Note that my proposal
        # distribution is not symmetric since I'm using a truncated normal
        log_post_prop = logpost(theta = theta_prop, 
                                sigma2 = sigma2_use, 
                                data = data, 
                                P = P,  missing_index,
                                DTP0 = theta_prop[DTP0_index],
                                temperature = temperatures[chain],
                                filename = filename)  # proposed
        logalpha =   log_post_prop - 
          log_post_old + 
          sum(log(dtruncnorm(theta_use[DTP0_index],  
                             a=0, b=Inf, 
                             mean = DTP0prop, 
                             sd = step_var[[chain]][5]))) -
          sum(log(dtruncnorm(DTP0prop, 
                             a=0, b=Inf, 
                             mean = theta_use[DTP0_index],  
                             sd = step_var[[chain]][5])))
        
        if(!is.na(logalpha) && runif(1) < exp(logalpha)){
          accepts[[chain]][current_stop,Nparblocks] = accepts[[chain]][current_stop,Nparblocks]+1;
          theta_use     = theta_prop
          log_post_old  = log_post_prop
        }  
      }
      
      #### Done those Gibbs {A,B,C} steps, lock in the resulting theta vector ####
      theta[[chain]][iter+1,] = theta_use
      
      
      #### Sample from sigma^2 | everything else
      #
      # This can probably be done directly from the conditional distribution, 
      # This is not wrong, but maybe inefficient
      #
      
      ##### sigma2 #####
      # propose a value from an easy distribution: truncated Normal.  
      sigma2_prop          = rtruncnorm(n = 1, a=0, b=Inf, mean = sigma2_use, sd = step_var[[chain]][4])
      # the ratio of un-normalized posteriors.  Note that my proposal
      # distribution is symmetric so Q_{ij}=Q_{ji}
      
      log_post_prop = logpost(theta = theta_use, 
                              sigma2 = sigma2_prop, 
                              data = data, 
                              P = P,  
                              missing_index,
                              DTP0 = theta_use[DTP0_index],
                              temperature = temperatures[chain],
                              filename = filename)  # proposed
      logalpha =   log_post_prop - 
        log_post_old + 
        log(dtruncnorm(sigma2_use,  a=0, b=Inf, mean = sigma2_prop, sd = step_var[[chain]][4])) -
        log(dtruncnorm(sigma2_prop, a=0, b=Inf, mean = sigma2_use,  sd = step_var[[chain]][4]))
      if(!is.na(logalpha) && runif(1) < exp(logalpha)){
        accepts[[chain]][current_stop,4] = accepts[[chain]][current_stop,4]+1;
        sigma2[[chain]][iter+1] = sigma2_prop
        log_post_old            = log_post_prop
      }else{
        sigma2[[chain]][iter+1] = sigma2[[chain]][iter]
      }
      log_posterior[[chain]][iter+1] = log_post_old
    }# end loop over NChain
    
    # unlist(lapply(log_posterior, function(x){x[iter+1]}))
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
        # modify if it has been run for a while, but stop before it's been running for too long.  & Also make sure that there have been some accepted values:
        if( current_stop>2 &&
            iter<= niter/2 && 
            nrow( unique(theta[[chain]][floor(iter/2):iter,C_index]))>2){ # less often make an adjustment to the correlation structure: 
          CCor2[[chain]] = cor(theta[[chain]][floor(iter/2):iter,C_index])*.8 # shrink correlations down since I'm sure they're wrong
          diag(CCor2[[chain]]) = 1
          CCor[[chain]] = (CCor[[chain]] + CCor2[[chain]]) / 2 # average them so that it moves in a less twitchy manner.
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
    sigma2               = sigma2,
    DTP_missing          = DTP_missing,
    theta_labels         = theta_labels,
    A_index              = A_index,
    B_index              = B_index,
    C_index              = C_index,
    DTP0_index           = DTP0_index
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
  
  
  return(list(dates_in_use =   data |> pull(date),
              data =   data |> select(-c("date", "anchor")) |> as.matrix()))
}







# ## logprior - prior on sigma2 is  exponential(1)
# logprior = function(theta,sigma2,DTP0=NULL, filename = NULL){
#   # $$\sigma\overset{iid}{\sim }exponential(1)$$  
#   if(sigma2 <= 0){ #case 1 negative variance
#     # negative variance is bad
#     logprior = -Inf
#   }else{
#     if((is.null(DTP0) | length(DTP0)==0) ){ #positive variance, no missing values at time zero
#       # no starting point to deal with,
#       logprior = sum(dnorm(theta, log = TRUE)) -sigma2
#     }else{#positive variance, with missing values at time zero
#       if(ifelse(is.null(filename), TRUE, !grep(filename, pattern = "^log", ignore.case = TRUE))){#positive variance, with missing values at time zero, and not taking the log of observations
#         # null filename or not starting with "log" flag:
#         # there is a starting point DTP0 to deal with...
#         # if we are using the log of the data (removes the negative constraint on starting point)
#         # do this by checking the filename for the starting term "log"
#         if(min(DTP0)>=0 && max(DTP0)<=1 ){#positive variance, with missing values at time zero, and not taking the log of observations, and starting point is between zero and 1
#           # prior for time zero evaluation is uniform(0,1)
#           logprior = sum(dnorm(theta, log = TRUE)) - sigma2
#         }else{#positive variance, with missing values at time zero, and not taking the log of observations, and starting point is NOT between zero and 1
#           # DTP0 is too large or too small
#           logprior = -Inf
#         }
#       }else{#positive variance, with missing values at time zero, and  taking the log of observations
#         # there is a starting point DTP0 to deal with
#         logprior = sum(dnorm(theta, log = TRUE)) +
#           sum(dnorm(DTP0, log = TRUE, mean = -3, sd = 4)) - ##### mean is ~log(.05) and mean-2sd is ~log(.0009), mean+2sd is ~log(3)=1.1
#           sigma2
#       }
#     }
#   }
#   return(logprior)
# }
# 
# # 

## logprior
logprior = function(theta,sigma2,DTP0=NULL, filename = NULL){
  # \sigma2 ~ Exponential with mean .01
  if(sigma2 <= 0){ #case 1 negative variance
    # negative variance is bad
    logprior = -Inf
  }else{
    if((is.null(DTP0) | length(DTP0)==0) ){ #positive variance, no missing values at time zero
      # no starting point to deal with,
      logprior = sum(dnorm(theta, log = TRUE)) -
        sigma2/.01
    }else{#positive variance, with missing values at time zero
      if(ifelse(is.null(filename), TRUE, !grep(filename, pattern = "^log", ignore.case = TRUE))){#positive variance, with missing values at time zero, and not taking the log of observations
      # null filename or not starting with "log" flag:
      # there is a starting point DTP0 to deal with...
      # if we are using the log of the data (removes the negative constraint on starting point)
      # do this by checking the filename for the starting term "log"
      if(min(DTP0)>=0 && max(DTP0)<=1 ){#positive variance, with missing values at time zero, and not taking the log of observations, and starting point is between zero and 1
        # prior for time zero evaluation is uniform(0,1)
        logprior = sum(dnorm(theta, log = TRUE)) -
          sigma2/.01
      }else{#positive variance, with missing values at time zero, and not taking the log of observations, and starting point is NOT between zero and 1
        # DTP0 is too large or too small
        logprior = -Inf
       }
        }else{#positive variance, with missing values at time zero, and  taking the log of observations
            # there is a starting point DTP0 to deal with
            logprior = sum(dnorm(theta, log = TRUE)) +
              sum(dnorm(DTP0, log = TRUE, mean = -3, sd = 4)) - ##### mean is ~log(.05) and mean-2sd is ~log(.0009), mean+2sd is ~log(3)=1.1
              sigma2/.01
        }
        }
      }
  return(logprior)
}
############################
## model_propagate
model_propagate = function(theta,sigma2,data,P = 1){
  
  # if missing values are in too early then the format passed in may not be a matrix
  if(!is.matrix(data)){
    data = matrix(data, ncol = length(data))
  }
  # P = lag depth
  # predict the model ahead by one time increment
  #$$DTP_{ijt} = a_j + \sum_pb_{jp} * DTP_{ijt-p} + \sum_{k\in \{1,...,J\} \setminus \{j\}}c_{jk} * DTP_{ikt-1}   + d * X_{i,j,t}  + e_{i,j,t}#
  J     = ncol(data) #Ncountries
  error = rnorm(J,mean = 0,sd = sqrt(sigma2))
  Number_of_a = J
  Number_of_b = J*P
  Number_of_c = (J-1)*J# We don't go into lags larger than 1 for the other countries
  a     = theta[1:Number_of_a] |> matrix(nrow = Number_of_a)
  # cols of B are lag.  Rows of b are country
  b     = theta[(Number_of_a+1): (Number_of_a+Number_of_b)] |> matrix(nrow = Number_of_a, byrow = TRUE)
  
  ctemp = theta[(Number_of_a+Number_of_b+1):(Number_of_a+Number_of_b+Number_of_c)]|> matrix(nrow = Number_of_a, byrow = TRUE) 
  # for the sake of matrix multiplication further down, make this a square with diagonal zeros.
  cmat = matrix(NA, ncol = J, nrow = J)
  # Fill upper triangle from  (excluding diagonal)
  cmat[upper.tri(cmat)] = ctemp[upper.tri(ctemp, diag = TRUE)]
  
  # Fill lower triangle from M (excluding diagonal)
  cmat[lower.tri(cmat)] = ctemp[lower.tri(ctemp, diag = FALSE)]
  
  # Fill diagonal with D
  diag(cmat) = 0
  
  
  
  #$$DTP_{ijt} = a_j + \sum_pb_{jp} * DTP_{ijt-p} + \sum_{k\in \{1,...,J\} \setminus \{j\}}c_{jk} * DTP_{ikt-1}   + d * X_{i,j,t}  + e_{i,j,t}#
  # so that we take in the most recently evaluated data and put the lagged version in the right place
  data_index = nrow(data):(nrow(data)-(P-1))
  predicted_mean = a + 
    apply(b*(data[data_index,]),1,sum) +
    cmat%*%data[nrow(data),]
  return(list(predicted_mean = predicted_mean,
              error = error,
              prediction = predicted_mean+error))
}

## model_predict
model_predict = function(theta,sigma2,data,P = 1,DTP0 = NULL){

  data_full = data_full_mean = data
  # Using the past, predict one step ahead to infill the present. 
  if(sum(is.na(data_full[1,]))>0){data_full[1,is.na(data_full[1,])] = DTP0}
  for(time_point in (P+1):nrow(data_full)){
    modepred = model_propagate(theta,sigma2,data_full[1:(time_point-1),],P = P)
    data_full[time_point,]           = modepred$prediction
    data_full_mean[time_point,]      = modepred$predicted_mean
  }
  return(list(data_full_prediction = data_full,
              data_full_mean       = data_full_mean))
}

## infill_missing
infill_missing = function(theta,sigma2,data,P = 1){
  
  # pore-allocate and find the missing values
  data_full = data
  missing_index_time = is.na(data) |> apply(1,any) |> which()
  
  # Using the past, predict one step ahead to infill the present. 
  for(time_point in missing_index_time){
    country_missing = which(is.na(data_full[time_point,]))
    
    data_full[time_point,country_missing] = model_propagate(theta,sigma2,data_full[1:(time_point-1),],P = 1)$prediction[country_missing]
    
  }
  return(data_full)
}





## loglikelihood
# Log likelihood is calculated for the observed data conditional on the infilled missing data.
loglikelihood = function(theta,sigma2,data,P = 1, missing_index,DTP0 = NULL){
  time = nrow(data)
  loglik = 0
  prop_model = model_predict(theta,sigma2,data,P=P,DTP0=DTP0)
    loglik = loglik + sum(
      dnorm(mean = prop_model$data_full_prediction[!missing_index], 
            sd   = sqrt(sigma2), 
            data[!missing_index], 
            log  = TRUE) # exclude missing values in the likelihood calculation since those are handled elsewhere.
    )
  return(loglik)
}






## logpost
logpost = function(theta,sigma2,data,P = 1,
                   missing_index,
                   DTP0=NULL,
                   temperature = 1,
                   filename = NULL){
  
  temperature*sum(loglikelihood(theta,sigma2,data,P = P, missing_index,DTP0)) + logprior(theta,sigma2,DTP0, filename)
  
  
}









## make_fake_data
make_fake_data = function(Nmissing = 25, Number_of_DTP0 = 1, Ncountries = 5, Ntimes = 200){

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
    
    # deal with starting time point missing values:

    DTP0 = runif(n = Number_of_DTP0)
    theta_true  =  c(rnorm((Number_of_a+Number_of_b+Number_of_c), sd = .1),
                 DTP0)# 
    sigma2_true = .5

    data[1,DTP0_missing_index] = theta_true[(Number_of_a+Number_of_b+Number_of_c+1):(Number_of_a+Number_of_b+Number_of_c+Number_of_DTP0)]

    data = infill_missing(theta_true,sigma2_true,data,P = 1)
    data[missing_index] = NA
    data[1,DTP0_missing_index] = NA
    
    return(list(
      Number_of_a = Number_of_a,
      Number_of_b = Number_of_b,
      Number_of_c = Number_of_c,
      Number_of_DTP0   = Number_of_DTP0,
      DTP0_missing_index = DTP0_missing_index,
      data = data,
      missing_index = missing_index,
      Nmissing = Nmissing,
      Ncountries = Ncountries,
      Ntimes = Ntimes,
      sigma2_true = sigma2_true,
      theta_true = theta_true
    )
    )
}


















## ----version----------------------
#| code-fold: true
version
library(tidyverse)
library(lubridate)
library(rstan)
options(mc.cores = 4)
rstan_options(auto_write = TRUE)

source("MCMC_SS_functions.R")     # Just used for data loading



## ----user_defined_options---------
# 1. location and filename
datafolder = "cbspeeches-main/inst/data-misc/"
datafile   = paste0(datafolder,"nmf-time-series-g20.csv")

#2. topic of interest
candidate_topics = c("inflat" ,"brexit", "covid", "cbdc", "ukrain" )
topic_of_interest = candidate_topics[1]

#3 countries to use.  
# possible option
G7_countries = c ("Canada", "France", "Germany", "Italy", "Japan",  "United Kingdom","United States")
# possible option
CUSA = c ("Canada","United States")
# actual decision:
countries2use = G7_countries# unique(dataset$country)

#4 first month to use, 
# start of dataset
start_date = ymd("2008-11-01")
# otherwise the month before the first mention in the dataset.
if(topic_of_interest=="brexit"){
  start_date = ymd("2016-02-01")
}
if(topic_of_interest=="covid"){
  start_date = ymd("2020-01-01")
}
if(topic_of_interest=="cbdc"){
  start_date = ymd("2016-02-01")
}
if(topic_of_interest=="ukrain"){
  start_date = ymd("2014-03-01")
}

#5  prefix for naming output filenames
prefix = "stan_" 



## ----dataload---------------------
#| warning: false

#. just loading and maybe transforming the data:
data_full = dataload(datafile, topic_of_interest, countries2use, start_date)
dates_in_use = data_full$dates_in_use
if(prefix == "log"){
  data         = log(data_full$data)
}else{
  data = data_full$data
}


## ----preamble---------------------
#| code-fold: false
if(!exists("niter")){ niter = 50000 }    # total number of iterations



## ----MCMC-------------------------
#| code-fold: true
#| eval: true


filename = paste0(prefix,"dataSS_",topic_of_interest,"_results_",paste0(countries2use, collapse = "_"),"_",niter,".Rdata")
print(filename)




T1 = Sys.time()

# Set up the data so that NAs and zeros are treated the same way.
data[is.na(data)] = 0
colnames(data) = NULL

# Prepare list for Stan
stan_data <- list(
  T = nrow(data),
  K = ncol(data),
  I = diag(ncol(data)),
  y_obs = data
)

# Compile and run
fit <- stan(
  file = "state_space.stan",
  data = stan_data,
  iter = niter, 
  init = "vb",
  chains = 4
)



elapsed = Sys.time() - T1 
saveRDS(PT_results$CCor, paste0("CCor_",gsub(filename, pattern = "Rdata", replacement = "rds")))
cat(paste("total compute time: ", round(as.numeric(elapsed, units = "secs"),2)," seconds"))
save.image(filename)




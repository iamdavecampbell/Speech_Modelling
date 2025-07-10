
## VAR model with missing observations

We model the topic proportion for topic **i**, from country **j** at time **t**, $DTP_{ijt}$, as a function of lagged information, $DTP_{ijt-1}$, information from countries other than country $j\in 1,...,J$, $DTP_{i(-j)t-1}$, and exogenous variables, $X_{i,j,t}$,

<!-- $$DTP_{ijt} = a + b * DTP_{ijt-1} + c * DTP_{i(-j)t-1}   + d * X_{i,j,t} + f * Country \times Month_{j,t} + e$$ -->
$$DTP_{ijt} = a_j + \sum_pb_{jp} * DTP_{ijt-p} + \sum_{k\in \{1,...,J\} \setminus \{j\}}c_{jk} * DTP_{ikt-1}     + e_{i,j,t}, \ \ \ where\ \ \ e\sim N(0,\sigma^2_{i})$$

- MCMC_functions.R = functions for performing MCMC on a VAR model.  Can also simulate  data
- Data_PTModelling.qmd = runs MCMC on data


## State Space Formulation:

Observation process


$$DTP_{ijt} = \alpha_{ijt} \left[X_{ijt} + e_{i,j,t}\right], \ \ \ where\ \ \ e_{ijt}\sim N(0,\sigma^2_{e,i}), \ \ \ and \ \ \ \alpha_{ijt} \sim Bernoulli(p_a)$$

 Un-observable Transition process:

Dropping the reliance on the topic index _i_, the model transitions ahead based on the  stochastic process: 

$$X_{jt} = a_j+b_{j} * X_{jt-1} + \sum_{k\in \{1,...,J\} \setminus \{j\}}c_{jk} * X_{kt-1}+ \delta_{jt}.$$

- MCMC_SS_functions.R = functions for performing MCMC on a state space model.  
- Data_PTModelling_state_space.qmd = runs MCMC on data

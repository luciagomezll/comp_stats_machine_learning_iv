#################################################################################
### SIMULATION OF THE STUDY OF RETURNS TO COLLEGE
#################################################################################

setwd("C:/Users/user/Documents/GitHub/comp_stats_machine_learning_iv")
rm(list=ls())

# Load libraries
if (!require("haven")) {install.packages("haven")}
library(haven)
set.seed(12345)

#=========================================================
# a. Generate the simulation of the data set
#=========================================================

local_level_dta <- read_dta('localvariables.dta')
basic_level_dta <- read_dta('basicvariables.dta')

# parameters for simulations
obs <- 1000
p <- 5
mean_mu <- c(0.0)
sigma_mu <- matrix(c(1,0.6,0.6,1),2)
obs_local <- obs/10

# DGP function 
dgp <- function(individual_obs,local_obs,individual_data,local_data,continuous_iv=FALSE,weak=FALSE,long=FALSE){
        if (continuous_iv== TRUE) {
                X <- matrix(rnorm(obs*p,0,1),obs,p)
                u <- mvrnorm(obs,mean_mu,sigma_mu)
                X.mod = as.matrix(cbind(rep(1,num.obs.students),W,X))
                Y <- X.mod %*% ols.beta.original + eps
        }
        else{
                X.mod = as.matrix(cbind(rep(1,num.obs.students),W,X))
                Y <- X.mod %*% ols.beta.original +X$X1.sim*W + X$X2.sim*W + eps
        }       
        
        
        
        # Simple setting
        
        # Propensity score
        
        # Covariates
        # basic variables
 #       x0 <- rep(1,obs)
 #       x1_cafqt <- rnorm(obs,mean(basic_level_dta$cafqt),sd(basic_level_dta$cafqt)) 
 #       x2_msch  <- rnorm(obs,mean(basic_level_dta$mhgc),sd(basic_level_dta$mhgc))
 #       x3_numsibs <- rnorm(obs,mean(basic_level_dta$numsibs),sd(basic_level_dta$numsibs)) # It seems to follow and F-distribution
 #       x4_urban14 <- rbinom(obs,1,mean(basic_level_dta$urban14))
 #       prob.x_d <- as.vector(colMeans(basic_level_dta[,c(6:12)])) 
 #       x5_d  <- sample(c(57:63), obs, replace=TRUE, prob=prob.x_d)
 #       x6_exp <- rnorm(obs,mean(basic_level_dta$exp),sd(basic_level_dta$exp))
 #       x7_expsq <- sqrt(x6_exp) 
 #       # Create local ID to match individuals with locals
 #       obs_per_local <- sample(1:50, obs_local, replace=TRUE)
 #       prop.obs_local <- obs_per_local/sum(obs_per_local)
 #       local.id <- sample(1:obs_local,obs,replace=TRUE, prob= prop.obs_local)
 #       data.indv <- cbind.data.frame(local.id,x0,x1_cafqt,x2_msch,x3_numsibs,
                                      x4_urban14,x5_d,x6_exp,x7_expsq)
        
        # local variables
  #      avlocwage17 <- rnorm(obs_local,mean(local_level_dta$lavlocwage17),sd(local_level_dta$lavlocwage17))
  #      avurate <- rnorm(obs_local,mean(local_level_dta$avurate),sd(local_level_dta$avurate))
  #      lwage_14 <- rnorm(obs_local,mean(local_level_dta$lwage5),sd(local_level_dta$lwage5))
  #      lurate_14 <- rnorm(obs_local,mean(local_level_dta$lurate),sd(local_level_dta$lurate))
  #      lurate_17 <- rnorm(obs_local,mean(local_level_dta$lurate_17),sd(local_level_dta$lurate_17))
  #      tuit4c <- rnorm(obs_local,mean(local_level_dta$tuit4c),sd(local_level_dta$tuit4c))
  #      local.id <- sample(1:obs_local,obs_local,replace=FALSE)
        # Instruments
  #      pub4 <- rbinom(obs_local,1,mean(local_level_dta$pub4))
  #      lwage5_17 <- rnorm(obs_local,mean(local_level_dta$lwage5_17),sd(local_level_dta$lwage5_17))
        
  #      data.local <- cbind.data.frame(local.id,avlocwage17,avurate,lwage_14,
 #                                      lurate_14,lwage5_17,lurate_17,tuit4c,pub4)
        
        # Assign groups
  #      x.raw <- merge(data.indv,data.local,by="local.id")        
        
        # Treatment (vincular con los instrumentos)
  #      eps <- rnorm(obs,0,1)
  #      eps_x <- rnorm(obs,0,1)
  #      x8_omitted <- 1 + eps_x
  #      S <- x.raw$pub4+x.raw$lwage5_17+ eps_x
  #      S_dummy <- (S > mean(S))*1
        
        # Potential outcomes
  #      Y <- beta_0 + beta_1*
#}



Ash et al. (2018) and Chernozhukov et al. (2018a) use random forest and related methods to
select instruments for IV estimators of heavily over-identified models.







#################################################################################
### SIMULATION OF THE STUDY OF RETURNS TO COLLEGE
#################################################################################

setwd("C:/Users/user/Documents/GitHub/comp_stats_machine_learning_iv")
rm(list=ls())

# Load libraries
if(!require(haven)) {install.packages("haven")}
if(!require(MASS)) {install.packages("MASS")}
if(!require(ivreg)) install.packages("ivreg")
if(!require(randomForest)) install.packages("randomForest")
if(!require(gbm)) install.packages("gbm")
library(randomForest)
library(ivreg)
library(ivprobit)
library(MASS)
library(mvtnorm)
library(haven)
library(gbm)
library(caret)
install.packages("newboost", repos="http://R-Forge.R-project.org")
library(newboost)
set.seed(12345)


#=========================================================
# a. Generate the simulation of the data set
#=========================================================

local_level_dta <- read_dta('localvariables.dta')
basic_level_dta <- read_dta('basicvariables.dta')

# parameters for simulations
n_simulation <- 100
n_observations <- 1000
n_exogeneous_variables <- 5 # number of exogenous variables
mean_eps <- c(0,0)
sigma_eps <- matrix(c(1,0.6,0.6,1),2)
#obs_local <- obs/10

# DGP function 
#dgp <- function(individual_obs,local_obs,individual_data,local_data,continuous_iv=FALSE,weak=FALSE,long=FALSE){
#        if (continuous_iv== TRUE) {
                #x.2 <- matrix(rnorm(obs*n_exo,0,1),obs,n_exo)
dgp <- function(obs,mean_mu,sigma_mu){
                x.2 = rnorm(obs,0,1)
                u = mvrnorm(obs*2,mean_mu,sigma_mu)
                #z <- matrix(rnorm(obs*n_inst,0,1),obs,n_inst)
                z = rnorm(obs,0,1)
                #x.1 = (0.1 + 1.5*z - 0.8*x.2 + u[,2]  0)*1
                x.1 = 0.1 + 1.5*z - 0.8*x.2 + u[,2] 
                y = 0.1 + 1*x.1 - 0.5*x.2 + u[,1]
                data <- data.frame("y"=y, "x.1"=x.1,"x.2"=x.2,"z"=z)
                return(data)
        }

dgp_v2 <- function(obs,mean_mu,sigma_mu){
        err_c=rnorm(obs, 0, 1)
        theta=runif(obs, -1, 1)
        x.2 = 1+err_c 
        x.3 = rnorm(obs,0,1)
        u = mvrnorm(obs*2,mean_mu,sigma_mu)
        #z <- matrix(rnorm(obs*n_inst,0,1),obs,n_inst)
        z = rnorm(obs,0,1)
        x.1 = (0.1 + 1.5*z - 0.8*x.3 + err_c > 0.8)*1
        y = 0.1 + 1*x.1 +1*x.2 - 0.5*x.3 + theta
        data <- data.frame("y"=y, "x.1"=x.1,"x.2"=x.2,"z"=z)
        return(data)
}

train_data <- dgp(obs,mean_mu,sigma_mu)
test_data <- dgp(obs,mean_mu,sigma_mu)

naive_ols <- lm(y ~ x.1 + x.2,data=train_data)
train_MSE_OLS = mean(naive_ols$residuals^2)

test_mse<- function(test_data, model){
        value <- mean((test_data$y - predict(model))^2)
        return(value)
}

test_MSE_OLS = test_mse(test_data,lm(y ~ x.1 + x.2,data=test_data))

iv <- ivreg(y ~ x.1 + x.2| z + x.2,data=train_data)
summary(iv)
train_mse_iv = mean(iv$residuals^2)
test_mse_iv <- test_mse(test_data, ivreg(y ~ x.1 + x.2| z + x.2,data=test_data))

trees <- 100
rf <- randomForest(x.1 ~ z, importance=TRUE, data=train_data, ntree=trees)
x_hat <- predict(rf)
x_hat_test <- predict(rf, newdata=test_data)

train_mse_rf <- mean(rf$mse) 
test_mse_rf <- mean((x_hat_test - test_data$x.1)^2)

boost <- gbm(x.1 ~ z, data=train_data, distribution = "gaussian", 
             n.trees = trees, interaction.depth = 4)
x_hat_train <- predict(boost, n.trees = trees)
xhat_boost_test<-predict(boost, newdata= test_data, n.trees = trees)

train_mse_b <- mean(boost$train.error)
test_mse_b <- mean((xhat_boost_test - test_data$x.1)^2) 
teb <- L2BoostEffect(x=xs, y=ys, d=ds, post=FALSE)
tebp <- L2BoostEffect(x=xs, y=ys, d=ds)
tebo <- orthoL2BoostEffect(x=xs, y=ys, d=ds)


        
        
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







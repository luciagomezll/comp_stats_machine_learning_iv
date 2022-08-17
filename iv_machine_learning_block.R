
# Objective: This script implements ensemble methods to conduct causal inference. 
# The script defines functions and codes which allow to sequentially execute the 
# different ensemble methods and estimate instrumental variables regressions eventually

#################################################################################
### IMPLEMENTATION OF THE ENSEMBLE METHODS TO ESTIMATE INSTRUMENTAL VARIABLES
#################################################################################

# Remove all objects from the workspace
rm(list=ls())

# Load libraries
if(!require(haven)) {install.packages("haven")}
if(!require(mvtnorm)) {install.packages("mvtnorm")}
if(!require(MASS)) {install.packages("MASS")}
if(!require(ivreg)) install.packages("ivreg")
if(!require(caret)) install.packages("caret")
if(!require(randomForest)) install.packages("randomForest")
if(!require(grf)) install.packages("grf")
if(!require(gbm)) install.packages("gbm")
if(!require(xgboost)) install.packages("xgboost")
if(!require(newboost)) install.packages("newboost", repos="http://R-Forge.R-project.org")
library(haven)
library(mvtnorm)
library(MASS)
library(ivreg)
library(caret)
library(randomForest)
library(gbm)
library(newboost)
library(grf)
library(xgboost)

# Set seed
set.seed(12345)

################################################################################
# CONTINUOUS ENDOGENOUS EXPLANATORY VARIABLE
################################################################################

# parameter setting
n_simulations <- 500
n_observations <- 1000
mean_eps <- c(0,0)
sigma_eps <- matrix(c(1,0.6,0.6,1),2)
beta_1 <- 1
beta_2 <- 0.8
beta <- c(beta_1,beta_2)

# function setting
dgp <- function(n_obs,beta,mean_eps,sigma_eps,covariates=FALSE){  # other_specification="FALSE"

        if (covariates== FALSE) {
        u = mvrnorm(n_obs,mean_eps,sigma_eps)
        z.1 = rnorm(n_obs,0,1)
        x.1 = 1 + 0.8*z.1 + u[,1] 
        y = beta[1] + beta[2]*x.1 + u[,2]
        data <- data.frame("y"=y,"x.1"=x.1,"z.1"=z.1)
        return(data)
        }
        else{
        n_covariates <- 10
        beta.x.2 = runif(n_covariates,min=-1,max=1)
        x.2 = matrix(rnorm(n_obs*n_covariates,0,1), n_obs, n_covariates)
        u = mvrnorm(n_obs,mean_eps,sigma_eps)
        z.1 = rnorm(n_obs,0,1)
        x.1 = 1 + 0.8*z.1 + u[,1] 
        y = beta[1] + beta[2]*x.1 + x.2%*%beta.x.2 + u[,2]
        data <- data.frame("y"=y, "x.1"=x.1,x.2,"z.1"=z.1)
        return(data)
        }
}

#x.1 = (0.1 + 1.5*z - 0.8*x.2 + u[,2]  0)*1
#z <- matrix(rnorm(obs*n_inst,0,1),obs,n_inst)

# A.1 One instrument no covariates
#===============================================================================

names <- c("ols","rf","boost")
for (i in 1:length(names)){
        paste0("mse_","names[i]") <- matrix(NA,n_simulations,2)
}
ols
rf
boost


mse_ols   <- matrix(NA,n_simulations,2)
#mse_iv   <- matrix(NA,n_simulations,2)
mse_rf    <- matrix(NA,n_simulations,2)
mse_boost <- matrix(NA,n_simulations,2)

beta_1_ols   <- matrix(NA,n_simulations,4)
beta_1_iv    <- matrix(NA,n_simulations,4)
beta_1_tsls  <- matrix(NA,n_simulations,4)
beta_1_rf    <- matrix(NA,n_simulations,4)
beta_1_boost <- matrix(NA,n_simulations,4)

bias_beta_1_ols   <- rep(NA,n_simulations)
bias_beta_1_iv    <- rep(NA,n_simulations)
bias_beta_1_tsls  <- matrix(NA,n_simulations,2)
bias_beta_1_rf    <- matrix(NA,n_simulations,2)
bias_beta_1_boost <- matrix(NA,n_simulations,2)


for (i in 1:n_simulations){        

        data <- dgp(n_observations,beta,mean_eps,sigma_eps)
        train_index <- sample(seq_len(nrow(data)), size = nrow(data)*0.7)
        train_data  <- data[train_index, ]
        test_data   <- data[-train_index, ]
        #plot(train_data$x.1,train_data$y)
       
        # OLS model
        ols <- lm(y ~ x.1,data=data)
        beta_1_ols[i,1] = coefficients(ols)[2] # coefficient beta_1 for ols
        beta_1_ols[i,2] = sqrt(diag(vcov(ols)))[2] # sd of coefficient beta_1 for ols
        bias_beta_1_ols[i] = beta[2]- beta_1_ols[i,1]
        
        # IV model
        iv <- ivreg(y ~ x.1 | z.1,data=data)
        beta_1_iv[i,1] = coefficients(iv)[2] # coefficient beta_1 for iv
        beta_1_iv[i,2] = sqrt(diag(vcov(iv)))[2] # sd of coefficient beta_1 for iv
        bias_beta_1_iv[i] = beta[2]- beta_1_iv[i,1]
        
        # Estimation of the first step with OLS
        ols <- lm(x.1 ~ z.1,data=train_data)
        train_data$pred_ols_x.1 = predict(ols,data=train_data)
        test_data$pred_ols_x.1   = predict(ols,newdata=test_data)
        mse_ols[i,1] = mean((train_data$x.1 - train_data$pred_ols_x.1)^2) #mse of ols in train data
        mse_ols[i,2] = mean((test_data$x.1 - test_data$pred_ols_x.1)^2)   #mse of ols in test data
        tsls_train <- lm(y ~ pred_ols_x.1,data=train_data)
        tsls_test  <- lm(y ~ pred_ols_x.1,data=test_data)
        beta_1_tsls[i,1] = coefficients(tsls_train)[2]
        beta_1_tsls[i,2] = sqrt(diag(vcov(tsls_train)))[2]
        beta_1_tsls[i,3] = coefficients(tsls_test)[2]
        beta_1_tsls[i,4] = sqrt(diag(vcov(tsls_test)))[2]
        bias_beta_1_tsls[i,1] = beta[2]- beta_1_tsls[i,1]
        bias_beta_1_tsls[i,2] = beta[2]- beta_1_tsls[i,3]        
        
        # Estimation of the first step with random forest
        rf <- randomForest(x ~ z,data=train_data,ntrees=200)
        train_data$predict_rf_x <- predict(rf) # Predict on out-of-bag training samples
        test_data$predict_rf_x <- predict(rf,newdata = test_data)
        mse_rf[i,1] = mean((train_data$x.1 - train_data$predict_rf_x.1)^2)
        mse_rf[i,2] =mean((test_data$x.1 - test_data$predict_rf_x.1)^2)
        tsls_rf_train <- lm(y ~ predict_rf_x.1,data=train_data)
        tsls_rf_test  <- lm(y ~ predict_rf_x.1,data=test_data)
        beta_1_rf[i,1] = coefficients(tsls_rf_train)[2]
        beta_1_rf[i,2] = sqrt(diag(vcov(tsls_rf_train)))[2]
        beta_1_rf[i,3] = coefficients(tsls_rf_test)[2]
        beta_1_rf[i,4] = sqrt(diag(vcov(tsls_rf_test)))[2]
        bias_beta_1_rf[i,1] = beta[2]- beta_1_rf[i,1]
        bias_beta_1_rf[i,2] = beta[2]- beta_1_rf[i,3] 
        
        #data$x.2 <- 0
        #forest.iv = instrumental_forest(X=data$x.1,Y=data$y,W=data$x.1,Z=data$z.1, reduced.form.weight = 0, mtry = p)
        #preds.iv = predict(forest.iv, test_data, estimate.variance = TRUE)$predictions
        #tau.hat = preds.iv$predictions
        #var.hat = preds.iv$variance.estimates
        #output = data.frame(dummy, TAU=tau.hat, VAR=var.hat)
        
        # Estimation of the first step with boost
        
        boost <- gbm(x.1 ~ z.1,data=train_data,distribution="gaussian",
                     n.trees=trees,interaction.depth = 4)
        train_data$predict_boost_x.1 <- predict(boost) # Predict on out-of-bag training samples
        test_data$predict_boost_x.1 <- predict(boost,newdata = test_data,n.trees = trees)
        mse_boost[i,1] = mean((train_data$x.1 - train_data$predict_boost_x.1)^2)
        mse_boost[i,2] =mean((test_data$x.1 - test_data$predict_boost_x.1)^2)
        tsls_boost_train <- lm(y ~ predict_boost_x.1,data=train_data)
        tsls_boost_test  <- lm(y ~ predict_boost_x.1,data=test_data)
        beta_1_boost[i,1] = coefficients(tsls_boost_train)[2]
        beta_1_boost[i,2] = sqrt(diag(vcov(tsls_boost_train)))[2]
        beta_1_boost[i,3] = coefficients(tsls_boost_test)[2]
        beta_1_boost[i,4] = sqrt(diag(vcov(tsls_boost_test)))[2]
        bias_beta_1_boost[i,1] = beta[2]- beta_1_boost[i,1]
        bias_beta_1_boost[i,2] = beta[2]- beta_1_boost[i,3] 
}

mse <- matrix(c(mean(mse_ols[,1]),mean(mse_ols[,2]),
                mean(mse_rf[,1]),mean(mse_rf[,2]),
                mean(mse_boost[,1]),mean(mse_boost[,2])),
                3,2,byrow=TRUE)
colnames(mse) <- c("Train","Test")
rownames(mse) <- c("OLS","RF","Boost")
mse_1

beta_1 <- matrix(c(mean(beta_1_ols[,1]),mean(beta_1_ols[,2]),mean(bias_beta_1_ols),
                   mean(beta_1_iv[,1]),mean(beta_1_iv[,2]),mean(bias_beta_1_iv),
                   mean(beta_1_tsls[,1]),mean(beta_1_tsls[,2]),mean(bias_beta_1_tsls[,1]),
                   mean(beta_1_tsls[,3]),mean(beta_1_tsls[,4]),mean(bias_beta_1_tsls[,2]),
                   mean(beta_1_rf[,1]),mean(beta_1_rf[,2]),mean(bias_beta_1_rf[,1]),
                   mean(beta_1_rf[,3]),mean(beta_1_rf[,4]),mean(bias_beta_1_rf[,2]),
                   mean(beta_1_boost[,1]),mean(beta_1_boost[,2]),mean(bias_beta_1_boost[,1]),
                   mean(beta_1_boost[,3]),mean(beta_1_boost[,4]),mean(bias_beta_1_boost[,2])),
                8,3,byrow=TRUE)
colnames(beta_1) <- c("Coefficient","Std. dev.","Bias") 
rownames(beta_1) <- c("OLS","IV","TSLS train","TSLS test","RF train","RF test","Boost train","Boost test")
beta_1

# A.2 One instrument with covariates
#===============================================================================

for (i in 1:n_simulations){        
        
        data <- dgp(n_observations,beta,mean_eps,sigma_eps,covariates=TRUE)
        train_index <- sample(seq_len(nrow(data)), size = nrow(data)*0.7)
        train_data  <- data[train_index, ]
        test_data   <- data[-train_index, ]
        #plot(train_data$x.1,train_data$y)
        
        #data.mod <- data[,-which(colnames(data)%in%"z.1")]
        y <- data$y
        x.1 <- data$x.1
        x.2 <- as.matrix(data[,-which(colnames(data)%in% c("y","z.1"))])
        z.1 <- data$z.1
        
        # OLS model
        ols <- lm(y ~ x.1 + x.2)
        beta_1_ols[i,1] = coefficients(ols)[2] # coefficient beta_1 for ols
        beta_1_ols[i,2] = sqrt(diag(vcov(ols)))[2] # sd of coefficient beta_1 for ols
        bias_beta_1_ols[i] = beta[2]- beta_1_ols[i,1]
        
        # IV model
        iv <- ivreg(y ~ x.1 + x.2 | z.1 + x.2)
        beta_1_iv[i,1] = coefficients(iv)[2] # coefficient beta_1 for iv
        beta_1_iv[i,2] = sqrt(diag(vcov(iv)))[2] # sd of coefficient beta_1 for iv
        bias_beta_1_iv[i] = beta[2]- beta_1_iv[i,1]
        
        # Estimation of the first step with OLS
        ols <- lm(x.1 ~ z.1,data=train_data)
        train_data$pred_ols_x.1 = predict(ols,data=train_data)
        test_data$pred_ols_x.1   = predict(ols,newdata=test_data)
        mse_ols[i,1] = mean((train_data$x.1 - train_data$pred_ols_x.1)^2) #mse of ols in train data
        mse_ols[i,2] = mean((test_data$x.1 - test_data$pred_ols_x.1)^2)   #mse of ols in test data
        tsls_train <- lm(y ~ . - z.1,data=train_data)
        tsls_test  <- lm(y ~ . - z.1,data=test_data)
        beta_1_tsls[i,1] = coefficients(tsls_train)[2]
        beta_1_tsls[i,2] = sqrt(diag(vcov(tsls_train)))[2]
        beta_1_tsls[i,3] = coefficients(tsls_test)[2]
        beta_1_tsls[i,4] = sqrt(diag(vcov(tsls_test)))[2]
        bias_beta_1_tsls[i,1] = beta[2]- beta_1_tsls[i,1]
        bias_beta_1_tsls[i,2] = beta[2]- beta_1_tsls[i,3]        
        
        # Estimation of the first step with random forest
        trees <- 400
        rf <- randomForest(x.1 ~ z.1,data=train_data,ntrees=trees)
        train_data$predict_rf_x.1 <- predict(rf) # Predict on out-of-bag training samples
        test_data$predict_rf_x.1 <- predict(rf,newdata = test_data)
        mse_rf[i,1] = mean((train_data$x.1 - train_data$predict_rf_x.1)^2)
        mse_rf[i,2] =mean((test_data$x.1 - test_data$predict_rf_x.1)^2)
        tsls_rf_train <- lm(y ~ . - z.1,data=train_data)
        tsls_rf_test  <- lm(y ~ . - z.1,data=test_data)
        beta_1_rf[i,1] = coefficients(tsls_rf_train)[2]
        beta_1_rf[i,2] = sqrt(diag(vcov(tsls_rf_train)))[2]
        beta_1_rf[i,3] = coefficients(tsls_rf_test)[2]
        beta_1_rf[i,4] = sqrt(diag(vcov(tsls_rf_test)))[2]
        bias_beta_1_rf[i,1] = beta[2]- beta_1_rf[i,1]
        bias_beta_1_rf[i,2] = beta[2]- beta_1_rf[i,3] 
        
        #data$x.2 <- 0
        #forest.iv = instrumental_forest(X=data$x.1,Y=data$y,W=data$x.1,Z=data$z.1, reduced.form.weight = 0, mtry = p)
        #preds.iv = predict(forest.iv, test_data, estimate.variance = TRUE)$predictions
        #tau.hat = preds.iv$predictions
        #var.hat = preds.iv$variance.estimates
        #output = data.frame(dummy, TAU=tau.hat, VAR=var.hat)
        
        # Estimation of the first step with boost
        
        boost <- gbm(x.1 ~ z.1,data=train_data,distribution="gaussian",
                     n.trees=trees,interaction.depth = 4)
        train_data$predict_boost_x.1 <- predict(boost) # Predict on out-of-bag training samples
        test_data$predict_boost_x.1 <- predict(boost,newdata = test_data,n.trees = trees)
        mse_boost[i,1] = mean((train_data$x.1 - train_data$predict_boost_x.1)^2)
        mse_boost[i,2] =mean((test_data$x.1 - test_data$predict_boost_x.1)^2)
        tsls_boost_train <- lm(y ~ . - z.1,data=train_data)
        tsls_boost_test  <- lm(y ~ . - z.1,data=test_data)
        beta_1_boost[i,1] = coefficients(tsls_boost_train)[2]
        beta_1_boost[i,2] = sqrt(diag(vcov(tsls_boost_train)))[2]
        beta_1_boost[i,3] = coefficients(tsls_boost_test)[2]
        beta_1_boost[i,4] = sqrt(diag(vcov(tsls_boost_test)))[2]
        bias_beta_1_boost[i,1] = beta[2]- beta_1_boost[i,1]
        bias_beta_1_boost[i,2] = beta[2]- beta_1_boost[i,3] 
}


mse <- matrix(c(mean(mse_ols[,1]),mean(mse_ols[,2]),
                mean(mse_rf[,1]),mean(mse_rf[,2]),
                mean(mse_boost[,1]),mean(mse_boost[,2])),
              3,2,byrow=TRUE)
colnames(mse) <- c("Train","Test")
rownames(mse) <- c("OLS","RF","Boost")
mse

beta_1 <- matrix(c(mean(beta_1_ols[,1]),mean(beta_1_ols[,2]),mean(bias_beta_1_ols),
                   mean(beta_1_iv[,1]),mean(beta_1_iv[,2]),mean(bias_beta_1_iv),
                   mean(beta_1_tsls[,1]),mean(beta_1_tsls[,2]),mean(bias_beta_1_tsls[,1]),
                   mean(beta_1_tsls[,3]),mean(beta_1_tsls[,4]),mean(bias_beta_1_tsls[,2]),
                   mean(beta_1_rf[,1]),mean(beta_1_rf[,2]),mean(bias_beta_1_rf[,1]),
                   mean(beta_1_rf[,3]),mean(beta_1_rf[,4]),mean(bias_beta_1_rf[,2]),
                   mean(beta_1_boost[,1]),mean(beta_1_boost[,2]),mean(bias_beta_1_boost[,1]),
                   mean(beta_1_boost[,3]),mean(beta_1_boost[,4]),mean(bias_beta_1_boost[,2])),
                 8,3,byrow=TRUE)
colnames(beta_1) <- c("Coefficient","Std. dev.","Bias") 
rownames(beta_1) <- c("OLS","IV","TSLS train","TSLS test","RF train","RF test","Boost train","Boost test")
beta_1

# B.1 Many instruments with strong sparcity no covariates
#===============================================================================

# We first consider the case in which the first stage is sparse, that is, there 
# are many available instruments but there exist only a small set of strong
# instruments whose identity is unknown. 

# Instruments not correlated to each other

# parameter setting

n_simulations <- 500
sample_size <- 1000
#sample_size <- 1500
mean_eps <- c(0,0)
sigma_eps <- 0.3
#sigma_ev <- 0.5

n_strong_instruments <- 25
#n_strong_instruments <- 50
beta_1 <- 1
beta_2 <- 0.8
beta <- c(beta_1,beta_2)

# function setting
dgp <- function(n_obs,
                beta,
                mean_eps,
                sigma_eps,
                n_strong_instruments) 
{                
                # Parameter setting
                n_instruments <- 500
                alpha_1 <- 5/n_strong_instruments
                sigma_ev <- matrix(c(1,sigma_eps,sigma_eps,1),2)
                h = c(rep(0, n_instruments - n_strong_instruments), rep(1,n_strong_instruments))
                
                u = mvrnorm(n_obs,mean_eps,sigma_ev)
                z = matrix(rnorm(n_obs*n_instruments,0,1), n_obs, n_instruments)
                x = alpha_1*z%*%h + u[,1] 
                y = beta[1] + beta[2]*x + u[,2]
                data <- data.frame("y"=y,"x"=x,z)
                return(data)
}



# Instruments correlated to each other



# B.1 Many instruments with strong sparcity with covariates
#===============================================================================

# parameter setting

n_simulations <- 500
sample_size <- 1000
#sample_size <- 1500
mean_eps <- c(0,0)
sigma_eps <- 0.3
#sigma_ev <- 0.5

beta_0 <- -0-90
beta_1 <- 0.75
beta <- c(beta_1,beta_2)
alpha_0 <- 0.3
alpha_1 <- 0.3
# alpha_1 <- 0.5
alpha <- c(alpha_0,alpha_1)

# function setting
dgp <- function(n_obs,beta,alpha,mean_eps,sigma_eps){  # other_specification="FALSE"
        
        n_instruments <- 500
        sigma_ev <- matrix(c(1,sigma_eps,sigma_eps,1),2)
        h = c(rep(1,n_instruments))
        
        u = mvrnorm(n_obs,mean_eps,sigma_ev)
        z = matrix(rnorm(n_obs*n_instruments,0,1), n_obs, n_instruments)
        x = alpha[1]+alpha[2]*z%*%h + u[,1] 
        y = beta[1] + beta[2]*x + u[,2]
        data <- data.frame("y"=y,"x"=x,z)
        return(data)
}


# Instruments not correlated to each other

# Instruments correlated to each other

# C. Many week with strong sparcity no covariates
#===============================================================================

#We next consider the case of a linear causal model with many instruments, all of which are week 
#(sparcity assumption breaks down). The issue with many weak instruments is
#that, when the sparsity assumption breaks down, variable selection methods like lasso or
#boosting tend not to select any variable or select all variables, which leads to poor asymptotics.

# Instruments not correlated to each other

# Many weak technical instruments. Constructing multiple polynomial functions 
# from the available exogeneous data
# parameter setting

n_simulations <- 500
n_observations <- 1000
#sample_size <- 1500
mean_eps <- c(0,0)
sigma_eps <- 0.3
#sigma_ev <- 0.5

beta_0 <- -0-90
beta_1 <- 0.75
beta <- c(beta_0,beta_1)
alpha_0 <- 0.3
alpha_1 <- 1
# alpha_1 <- 0.5
alpha <- c(alpha_0,alpha_1)

# function setting
dgp <- function(n_obs,beta,alpha,mean_eps,sigma_eps){  # other_specification="FALSE"

        #n_instruments <- 5
        sigma_ev <- matrix(c(1,sigma_eps,sigma_eps,1),2)

        u = mvrnorm(n_obs,mean_eps,sigma_ev)
        #z <-  matrix(rnorm(n_obs*n_instruments,0,1), n_obs, n_instruments)
        z = rnorm(n_obs,0,1)
        x = (alpha[1]+ alpha[2]*z + u[,1] > 0)*1 
        y = beta[1] + beta[2]*x + u[,2]
        data <- data.frame("y"=y,"x"=x,"z"=z)
        return(data)
}

# A.1 One instrument no covariates
#===============================================================================

mse_ols   <- matrix(NA,n_simulations,2)
#mse_iv   <- matrix(NA,n_simulations,2)
mse_logit <- matrix(NA,n_simulations,2)
mse_rf    <- matrix(NA,n_simulations,2)
mse_boost <- matrix(NA,n_simulations,2)
mse_xgboost <- matrix(NA,n_simulations,2)

beta_1_ols   <- matrix(NA,n_simulations,4)
#beta_1_iv    <- matrix(NA,n_simulations,4)
beta_1_tsls  <- matrix(NA,n_simulations,4)
beta_1_rf    <- matrix(NA,n_simulations,4)
beta_1_boost <- matrix(NA,n_simulations,4)
beta_1_xgboost <- matrix(NA,n_simulations,4)

bias_beta_1_ols   <- rep(NA,n_simulations)
#bias_beta_1_iv    <- rep(NA,n_simulations)
bias_beta_1_tsls  <- matrix(NA,n_simulations,2)
bias_beta_1_rf    <- matrix(NA,n_simulations,2)
bias_beta_1_boost <- matrix(NA,n_simulations,2)
bias_beta_1_xgboost <- matrix(NA,n_simulations,2)

for (i in 1:n_simulations){        
        
        data <- dgp(n_observations,beta,alpha,mean_eps,sigma_eps)
        train_index <- sample(seq_len(nrow(data)), size = nrow(data)*0.7)
        train_data  <- data[train_index, ]
        test_data   <- data[-train_index, ]
        #plot(train_data$x.1,train_data$y)
        
        # OLS model
        ols <- lm(y ~ x,data=data)
        beta_1_ols[i,1] = coefficients(ols)[2] # coefficient beta_1 for ols
        beta_1_ols[i,2] = sqrt(diag(vcov(ols)))[2] # sd of coefficient beta_1 for ols
        bias_beta_1_ols[i] = beta[2]- beta_1_ols[i,1]
        
        # Estimation of the first step with OLS
        logit <- glm(x ~ z, family = "binomial", data = train_data)
        probs.tr <- predict(logit,type ="response")
        x_hat.tr <- ifelse(probs.tr > 0.5,1,0)
        probs.ts <- predict(logit,newdata=test_data,type ="response")
        x_hat.ts <- ifelse(probs.ts > 0.5,1,0)
        mse_logit[i,1] = mean(train_data$x != x_hat.tr)  #mse of logit in train data
        mse_logit[i,1] = mean(test_data$x != x_hat.ts) #mse of logit in test data
        
        tsls_train <- lm(y ~ x_hat.tr,data=train_data)
        tsls_test  <- lm(y ~ x_hat.ts,data=test_data)
        beta_1_tsls[i,1] = coefficients(tsls_train)[2]
        beta_1_tsls[i,2] = sqrt(diag(vcov(tsls_train)))[2]
        beta_1_tsls[i,3] = coefficients(tsls_test)[2]
        beta_1_tsls[i,4] = sqrt(diag(vcov(tsls_test)))[2]
        bias_beta_1_tsls[i,1] = beta[2]- beta_1_tsls[i,1]
        bias_beta_1_tsls[i,2] = beta[2]- beta_1_tsls[i,3]   
        
        # Estimation of the first step with random forest
        trees <- 100
        rf <- randomForest(x ~ z,data=train_data,ntrees=trees)
        train_data$predict_rf_x <- predict(rf,type="class") # Predict on out-of-bag training samples
        x_hat.tr <- ifelse(train_data$predict_rf_x > 0.5,1,0)
        test_data$predict_rf_x <- predict(rf,type="class",newdata = test_data)
        x_hat.ts <- ifelse(test_data$predict_rf_x > 0.5,1,0)
        mse_rf[i,1] = mean(train_data$x != x_hat.tr)  #mse of logit in train data
        mse_rf[i,1] = mean(test_data$x != x_hat.ts) #mse of logit in test data
        
        tsls_rf_train <- lm(y ~ predict_rf_x,data=train_data)
        tsls_rf_test  <- lm(y ~ predict_rf_x,data=test_data)
        beta_1_rf[i,1] = coefficients(tsls_rf_train)[2]
        beta_1_rf[i,2] = sqrt(diag(vcov(tsls_rf_train)))[2]
        beta_1_rf[i,3] = coefficients(tsls_rf_test)[2]
        beta_1_rf[i,4] = sqrt(diag(vcov(tsls_rf_test)))[2]
        bias_beta_1_rf[i,1] = beta[2]- beta_1_rf[i,1]
        bias_beta_1_rf[i,2] = beta[2]- beta_1_rf[i,3] 
        
        #data$x.2 <- 0
        #forest.iv = instrumental_forest(X=data$x.1,Y=data$y,W=data$x.1,Z=data$z.1, reduced.form.weight = 0, mtry = p)
        #preds.iv = predict(forest.iv, test_data, estimate.variance = TRUE)$predictions
        #tau.hat = preds.iv$predictions
        #var.hat = preds.iv$variance.estimates
        #output = data.frame(dummy, TAU=tau.hat, VAR=var.hat)
        
        # Estimation of the first step with boost
        
        boost <- gbm(x ~ z,data=train_data,distribution="bernoulli",
                     n.trees=trees)
        train_data$predict_boost_x <- predict(boost) # Predict on out-of-bag training samples
        x_hat.tr <- ifelse(train_data$predict_boost_x > 0.5,1,0)
        test_data$predict_boost_x <- predict(boost,newdata = test_data,n.trees = trees)
        x_hat.ts <- ifelse(test_data$predict_boost_x > 0.5,1,0)
        mse_boost[i,1] = mean(train_data$x != x_hat.tr)  #mse of logit in train data
        mse_boost[i,1] = mean(test_data$x != x_hat.ts) #mse of logit in test data
        
        tsls_boost_train <- lm(y ~ predict_boost_x,data=train_data)
        tsls_boost_test  <- lm(y ~ predict_boost_x,data=test_data)
        beta_1_boost[i,1] = coefficients(tsls_boost_train)[2]
        beta_1_boost[i,2] = sqrt(diag(vcov(tsls_boost_train)))[2]
        beta_1_boost[i,3] = coefficients(tsls_boost_test)[2]
        beta_1_boost[i,4] = sqrt(diag(vcov(tsls_boost_test)))[2]
        bias_beta_1_boost[i,1] = beta[2]- beta_1_boost[i,1]
        bias_beta_1_boost[i,2] = beta[2]- beta_1_boost[i,3] 

        # Estimation of the first step with xgboost
        # Define predictor and response variables in training and test sets
        train_data_mod <- as.matrix(train_data)
        test_data_mod <- as.matrix(test_data)
        train_x <- as.vector(train_data_mod[,2])
        train_z <- cbind(rep(1,nrow(train_data_mod)),train_data_mod[,3])
        test_x <- as.vector(test_data_mod[2])
        test_z <- cbind(rep(1,nrow(test_data_mod)),test_data_mod[,3])
        # Define final training and testing sets
        #xgb_train <- xgb.DMatrix(data = train_z,label = train_x)
        #xgb_test <- xgb.DMatrix(data = test_z, label = test_x)
        #model <- xgb.train(data=xgb_train,max.depth=3,nrounds=70)
        #boost <- xgboost(data=xgb_train,max.depth=3,nrounds=70,distribution="bernoulli",
        #                 n.trees=trees)
    
        xg_boost <- xgboost(data=train_z,label=train_x,
                         max.depth=3,nrounds=70,objective = "binary:logistic",
                         verbose=FALSE)
        train_data$predict_xg_boost_x <- predict(xg_boost,data.matrix(train_z)) # Predict on out-of-bag training samples
        x_hat.tr <- ifelse(train_data$predict_xg_boost_x > 0.5,1,0)
        test_data$predict_xg_boost_x <- predict(xg_boost,data.matrix(test_z))
        x_hat.tr <- ifelse(train_data$predict_xg_boost_x > 0.5,1,0)
        mse_xgboost[i,1] = mean(train_data$x != x_hat.tr)  #mse of xgboost in train data
        mse_xgboost[i,1] = mean(test_data$x != x_hat.ts) #mse of xgboost in test data
        
        tsls_xgboost_train <- lm(y ~ predict_xg_boost_x,data=train_data)
        tsls_xgboost_test  <- lm(y ~ predict_xg_boost_x,data=test_data)
        beta_1_xgboost[i,1] = coefficients(tsls_xgboost_train)[2]
        beta_1_xgboost[i,2] = sqrt(diag(vcov(tsls_xgboost_train)))[2]
        beta_1_xgboost[i,3] = coefficients(tsls_xgboost_test)[2]
        beta_1_xgboost[i,4] = sqrt(diag(vcov(tsls_xgboost_test)))[2]
        bias_beta_1_xgboost[i,1] = beta[2]- beta_1_xgboost[i,1]
        bias_beta_1_xgboost[i,2] = beta[2]- beta_1_xgboost[i,3] 
}

mse <- matrix(c(mean(mse_logit[,1]),mean(mse_logit[,2]),
                mean(mse_rf[,1]),mean(mse_rf[,2]),
                mean(mse_boost[,1]),mean(mse_boost[,2]),
                mean(mse_xgboost[,1]),mean(mse_xgboost[,2])),
              4,2,byrow=TRUE)
colnames(mse) <- c("Train","Test")
rownames(mse) <- c("Logit","RF","Boost","Xgboost")
mse

beta_1 <- matrix(c(mean(beta_1_ols[,1]),mean(beta_1_ols[,2]),mean(bias_beta_1_ols),
                   mean(beta_1_tsls[,1]),mean(beta_1_tsls[,2]),mean(bias_beta_1_tsls[,1]),
                   mean(beta_1_tsls[,3]),mean(beta_1_tsls[,4]),mean(bias_beta_1_tsls[,2]),
                   mean(beta_1_rf[,1]),mean(beta_1_rf[,2]),mean(bias_beta_1_rf[,1]),
                   mean(beta_1_rf[,3]),mean(beta_1_rf[,4]),mean(bias_beta_1_rf[,2]),
                   mean(beta_1_boost[,1]),mean(beta_1_boost[,2]),mean(bias_beta_1_boost[,1]),
                   mean(beta_1_boost[,3]),mean(beta_1_boost[,4]),mean(bias_beta_1_boost[,2]),
                   mean(beta_1_xgboost[,1]),mean(beta_1_xgboost[,2]),mean(bias_beta_1_xgboost[,1]),
                   mean(beta_1_xgboost[,3]),mean(beta_1_xgboost[,4]),mean(bias_beta_1_xgboost[,2])),
                 9,3,byrow=TRUE)
colnames(beta_1) <- c("Coefficient","Std. dev.","Bias") 
rownames(beta_1) <- c("OLS","TSLS train","TSLS test","RF train","RF test",
                      "Boost train","Boost test","Xgboost train","Xgboost test")
beta_1




# /!\ Warning 
# Change tunning parameters \lambda por boosting regression

################################################################################
# BINARY ENDOGENEOUS EXPLANATORY VARIABLE
################################################################################

# A. One instrument no covariates
#===============================================================================

# B. Many instruments with strong sparcity no covariates
#===============================================================================

# Instruments not correlated to each other

# Instruments correlated to each other

# C. Many week with strong sparcity no covariates
#===============================================================================

# Instruments not correlated to each other

if (other_specification== FALSE) {
        u = mvrnorm(n_obs,mean_eps,sigma_eps)
        z.1 = rnorm(n_obs,0,1)
        x.1 = 1 + 0.5*z.1 + u[,1] 
        y = 1 + 0.3*x.1 + u[,2]
        data <- data.frame("y"=y, "x.1"=x.1,"z.1"=z.1)
        return(data)
}
else{
        u = mvrnorm(n_obs,mean_eps,sigma_eps)
        z.1 = rnorm(n_obs,0,1)
        x.1 = 1 + 0.5*z.1 + u[,1] 
        y = 1 + 0.3*x.1 + u[,2]
        data <- data.frame("y"=y, "x.1"=x.1,"z.1"=z.1)
        return(data)
}

#pred_ols_train_y = predict(ols,data=train_data)
#pred_ols_test_y = predict(ols,newdata=test_data)
#mse_ols[i,1] = mean((train_data$y - pred_ols_train_y)^2) #mse of ols in train data
#mse_ols[i,2] = mean((test_data$y - pred_ols_test_y)^2) # mse of ols in test data

#pred_iv_train_y = predict(iv,data=train_data)
#pred_iv_test_y  = predict(iv,newdata=test_data)
#mse_iv[i,1] = mean((train_data$y - pred_iv_train_y)^2)
#mse_iv[i,2] = mean((test_data$y - pred_iv_test_y)^2)

function(data,covariates=FALSE){
        if (covariates== FALSE) ols <- lm(y ~ x.1,data=data)
        else{
                X.mod = as.matrix(cbind(rep(1,num.obs.students),W,X))
                Y <- X.mod %*% ols.beta.original +X$X1.sim*W + X$X2.sim*W + eps
        }        
}

logit.coef <- logit$coefficients
z <- as.matrix(cbind(rep(1,nrow(train_data)),train_data$z))
pi_x.tr <- exp(z %*% as.vector(logit.coef))/(1 + exp(z %*% as.vector(logit.coef)))
pi_x.te <- exp(test_data$x %*% logit.coef)/(1 + exp(test_data$x %*% logit.coef))
x.tr_hat <- ifelse(pi_x.tr > 0.5,1,0)
Y.te_hat <- ifelse(pi_x.te > 0.5,1,0)
c.mse[i] <- mean(train_data$x!=x.tr_hat)
c.pe[i] <- mean(Y.te!=Y.te_hat)
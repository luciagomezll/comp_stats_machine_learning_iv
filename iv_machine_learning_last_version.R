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
#if(!require(grf)) install.packages("grf")
if(!require(gbm)) install.packages("gbm")
if(!require(xgboost)) install.packages("xgboost")
#if(!require(newboost)) install.packages("newboost", repos="http://R-Forge.R-project.org")
library(haven)
library(mvtnorm)
library(MASS)
library(ivreg)
library(caret)
library(randomForest)
#library(grf)
library(gbm)
library(xgboost)
#library(newboost)

# Setup-------------------------------------------------------------------------
# Set working directory using main project path
main_dir <- getwd()
code_dir <- file.path(main_dir,"Code")
data_dir <- file.path(main_dir,"Data")
output_dir <- file.path(main_dir,"Output")
# Set seed
set.seed(12345)

# Functions --------------------------------------------------------------------

# Data generating process with one instrumental variable

dgp_one_cont_iv <-  function(n_obs,
                             beta,
                             alpha,
                             mean_eps,
                             sigma_eps,
                             type_x)
{ 
        sigma_ev <- matrix(c(1,sigma_eps,sigma_eps,1),2)
        
        u = mvrnorm(n_obs,mean_eps,sigma_ev)
        z = rnorm(n_obs,0,1)
        if (type_x == "continuous"){
        x = alpha[1]+ alpha[2]*z + u[,1]
        } else if (type_x == "binary"){
        x = (alpha[1]+ alpha[2]*z + u[,1] > 0)*1     
        }
        y = beta[1] + beta[2]*x + u[,2]
        
        data <- data.frame("y"=y,"x"=x,"z"=z)
        return(data)
}

# Data generating process with many instrumental variables

dgp_many_cont_strong_iv <- function(n_obs,
                                    beta,
                                    alpha,
                                    mean_eps,
                                    sigma_eps,
                                    type_x,
                                    n_strong_instruments) {                
        # Parameter setting
        n_instruments <- 300
        sigma_ev <- matrix(c(1,sigma_eps,sigma_eps,1),2)
        if (n_strong_instruments == 0) {
                h <- c(rep(1,n_instruments))
        } else if (n_strong_instruments > 0) {
                h <- c(rep(0,n_instruments-n_strong_instruments), rep(1,n_strong_instruments)) 
        }
        
        # Defining the main variables
        u = mvrnorm(n_obs,mean_eps,sigma_ev)
        z = matrix(rnorm(n_obs*n_instruments,0,1),n_obs,n_instruments)
        if (type_x == "continuous"){
        x = alpha[1]+alpha[2]*z%*%h + u[,1]
        } else if (type_x == "binary"){
        x = (alpha[1]+alpha[2]*z%*%h + u[,1] > 0)*1
        }
        #x = alpha*z%*%h + u[,1] 
        y = beta[1] + beta[2]*x + u[,2]
        
        data <- data.frame("y"=y,"x"=x,z)
        return(data)
}


################# CONTINUOUS ENDOGENOUS EXPLANATORY VARIABLE ###################

# Case 1: One instrument no covariates -----------------------------------------

n_simulations <- 100
n_observations <- 1000
#sample_size <- 1500
mean_eps <- c(0,0)
sigma_eps <- 0.3
#sigma_eps <- 0.5
beta_0 <- -0-90
beta_1 <- 0.75
beta <- c(beta_0,beta_1)
alpha_0 <- 0.3
alpha_1 <- 0.3
# alpha_1 <- 0.5
alpha <- c(alpha_0,alpha_1)

mse_ols_cc1   <- matrix(NA,n_simulations,2)
mse_rf_cc1    <- matrix(NA,n_simulations,2)
mse_boost_cc1 <- matrix(NA,n_simulations,2)

beta_1_ols_cc1   <- matrix(NA,n_simulations,4)
beta_1_iv_cc1    <- matrix(NA,n_simulations,4)
beta_1_tsls_cc1  <- matrix(NA,n_simulations,4)
beta_1_rf_cc1    <- matrix(NA,n_simulations,4)
beta_1_boost_cc1 <- matrix(NA,n_simulations,4)

bias_beta_1_ols_cc1   <- rep(NA,n_simulations)
bias_beta_1_iv_cc1    <- rep(NA,n_simulations)
bias_beta_1_tsls_cc1  <- matrix(NA,n_simulations,2)
bias_beta_1_rf_cc1    <- matrix(NA,n_simulations,2)
bias_beta_1_boost_cc1 <- matrix(NA,n_simulations,2)

for (i in 1:n_simulations){        
        
        train_data <- dgp_one_cont_iv(n_observations,beta,alpha,mean_eps,sigma_eps,type_x="continuous")
        test_data <- dgp_one_cont_iv(n_observations,beta,alpha,mean_eps,sigma_eps,type_x="continuous")
#       train_index <- sample(seq_len(nrow(data)), size = nrow(data)*0.7)
#       train_data  <- data[train_index, ]
#       test_data   <- data[-train_index, ]

        # OLS model
        ols <- lm(y ~ x,data=test_data)
        beta_1_ols_cc1[i,1] = coefficients(ols)[2] 
        beta_1_ols_cc1[i,2] = sqrt(diag(vcov(ols)))[2] 
        bias_beta_1_ols_cc1[i] = beta[2]- beta_1_ols_cc1[i,1]
        
        # IV model
        iv <- ivreg(y ~ x|z,data=test_data)
        beta_1_iv_cc1[i,1] = coefficients(iv)[2]
        beta_1_iv_cc1[i,2] = sqrt(diag(vcov(iv)))[2] 
        bias_beta_1_iv_cc1[i] = beta[2]- beta_1_iv_cc1[i,1]
        
        # Estimation of the first stage with OLS
        ols <- lm(x ~ z,data=train_data)
        train_data$pred_ols_x = predict(ols,data=train_data)
        test_data$pred_ols_x  = predict(ols,newdata=test_data)
        mse_ols_cc1[i,1] = mean((train_data$x - train_data$pred_ols_x)^2) 
        mse_ols_cc1[i,2] = mean((test_data$x - test_data$pred_ols_x)^2)   
        # Estimation of the second stage with OLS
        tsls_train <- lm(y ~ pred_ols_x,data=train_data)
        tsls_test  <- lm(y ~ pred_ols_x,data=test_data)
        beta_1_tsls_cc1[i,1] = coefficients(tsls_train)[2]
        beta_1_tsls_cc1[i,2] = sqrt(diag(vcov(tsls_train)))[2]
        beta_1_tsls_cc1[i,3] = coefficients(tsls_test)[2]
        beta_1_tsls_cc1[i,4] = sqrt(diag(vcov(tsls_test)))[2]
        bias_beta_1_tsls_cc1[i,1] = beta[2] - beta_1_tsls_cc1[i,1]
        bias_beta_1_tsls_cc1[i,2] = beta[2] - beta_1_tsls_cc1[i,3]
        
        # Estimation of the first stage with random forest
        rf <- randomForest(x ~ z,data=train_data,ntrees=200)
        train_data$predict_rf_x <- predict(rf) # Predict on out-of-bag training samples
        test_data$predict_rf_x  <- predict(rf,newdata = test_data)
        mse_rf_cc1[i,1] = mean((train_data$x - train_data$predict_rf_x)^2)
        mse_rf_cc1[i,2] = mean((test_data$x - test_data$predict_rf_x)^2)
        # Estimation of the second stage with OLS
        tsls_rf_train <- lm(y ~ predict_rf_x,data=train_data)
        tsls_rf_test  <- lm(y ~ predict_rf_x,data=test_data)
        beta_1_rf_cc1[i,1] = coefficients(tsls_rf_train)[2]
        beta_1_rf_cc1[i,2] = sqrt(diag(vcov(tsls_rf_train)))[2]
        beta_1_rf_cc1[i,3] = coefficients(tsls_rf_test)[2]
        beta_1_rf_cc1[i,4] = sqrt(diag(vcov(tsls_rf_test)))[2]
        bias_beta_1_rf_cc1[i,1] = beta[2]- beta_1_rf_cc1[i,1]
        bias_beta_1_rf_cc1[i,2] = beta[2]- beta_1_rf_cc1[i,3]
        
        # Estimation of the first stage with boost
        boost <- gbm(x ~ z,data=train_data,distribution="gaussian",
                     n.trees=200,interaction.depth = 4)
        train_data$predict_boost_x <- predict(boost) # Predict on out-of-bag training samples
        test_data$predict_boost_x  <- predict(boost,newdata = test_data)
        mse_boost_cc1[i,1] = mean((train_data$x - train_data$predict_boost_x)^2)
        mse_boost_cc1[i,2] = mean((test_data$x - test_data$predict_boost_x)^2)
        # Estimation of the second stage with OLS
        tsls_boost_train <- lm(y ~ predict_boost_x,data=train_data)
        tsls_boost_test  <- lm(y ~ predict_boost_x,data=test_data)
        beta_1_boost_cc1[i,1] = coefficients(tsls_boost_train)[2]
        beta_1_boost_cc1[i,2] = sqrt(diag(vcov(tsls_boost_train)))[2]
        beta_1_boost_cc1[i,3] = coefficients(tsls_boost_test)[2]
        beta_1_boost_cc1[i,4] = sqrt(diag(vcov(tsls_boost_test)))[2]
        bias_beta_1_boost_cc1[i,1] = beta[2]- beta_1_boost_cc1[i,1]
        bias_beta_1_boost_cc1[i,2] = beta[2]- beta_1_boost_cc1[i,3] 
}

mse_cc1 <- matrix(c(mean(mse_ols_cc1[,1]),mean(mse_ols_cc1[,2]),
                 mean(mse_rf_cc1[,1]),mean(mse_rf_cc1[,2]),
                 mean(mse_boost_cc1[,1]),mean(mse_boost_cc1[,2])),
                 3,2,byrow=TRUE)
colnames(mse_cc1) <- c("Train","Test")
rownames(mse_cc1) <- c("OLS","RF","Boost")
mse_cc1

beta_1_cc1 <- matrix(c(mean(beta_1_ols_cc1[,1]),mean(beta_1_ols_cc1[,2]),mean(bias_beta_1_ols_cc1),
                     mean(beta_1_iv_cc1[,1]),mean(beta_1_iv_cc1[,2]),mean(bias_beta_1_iv_cc1),
                     mean(beta_1_tsls_cc1[,1]),mean(beta_1_tsls_cc1[,2]),mean(bias_beta_1_tsls_cc1[,1]),
                     mean(beta_1_tsls_cc1[,3]),mean(beta_1_tsls_cc1[,4]),mean(bias_beta_1_tsls_cc1[,2]),
                     mean(beta_1_rf_cc1[,1]),mean(beta_1_rf_cc1[,2]),mean(bias_beta_1_rf_cc1[,1]),
                     mean(beta_1_rf_cc1[,3]),mean(beta_1_rf_cc1[,4]),mean(bias_beta_1_rf_cc1[,2]),
                     mean(beta_1_boost_cc1[,1]),mean(beta_1_boost_cc1[,2]),mean(bias_beta_1_boost_cc1[,1]),
                     mean(beta_1_boost_cc1[,3]),mean(beta_1_boost_cc1[,4]),mean(bias_beta_1_boost_cc1[,2])),
                     8,3,byrow=TRUE)
colnames(beta_1_cc1) <- c("Coefficient","Std. dev.","Bias") 
rownames(beta_1_cc1) <- c("OLS","IV","TSLS train","TSLS test","TSLS-RF train","TSLS-RF test",
                          "TSLS-Boost train","TSLS-Boost test")
beta_1_cc1


# Case 2: Many instruments with strong sparcity no covariates ------------------

n_simulations <- 100
n_observations <- 1000
#sample_size <- 1500
n_strong_instruments <- 25
#n_strong_instruments <- 50
mean_eps <- c(0,0)
sigma_eps <- 0.3
#sigma_eps <- 0.5
beta_0 <- -0-90
beta_1 <- 0.75
beta <- c(beta_0,beta_1)
alpha_0 <- 0.3
alpha_1 <- 3/n_strong_instruments
# alpha_1 <- 6/n_strong_instruments
alpha <- c(alpha_0,alpha_1)

mse_ols_cc2   <- matrix(NA,n_simulations,2)
mse_rf_cc2    <- matrix(NA,n_simulations,2)
mse_boost_cc2 <- matrix(NA,n_simulations,2)

beta_1_ols_cc2   <- matrix(NA,n_simulations,4)
beta_1_iv_cc2    <- matrix(NA,n_simulations,4)
beta_1_tsls_cc2  <- matrix(NA,n_simulations,4)
beta_1_rf_cc2    <- matrix(NA,n_simulations,4)
beta_1_boost_cc2 <- matrix(NA,n_simulations,4)

bias_beta_1_ols_cc2   <- rep(NA,n_simulations)
bias_beta_1_iv_cc2    <- rep(NA,n_simulations)
bias_beta_1_tsls_cc2  <- matrix(NA,n_simulations,2)
bias_beta_1_rf_cc2    <- matrix(NA,n_simulations,2)
bias_beta_1_boost_cc2 <- matrix(NA,n_simulations,2)

for (i in 1:n_simulations){        

        train_data <- dgp_many_cont_strong_iv(n_observations,beta,alpha,mean_eps,sigma_eps,type_x="continuous",
                                        n_strong_instruments)  
        test_data <- dgp_many_cont_strong_iv(n_observations,beta,alpha,mean_eps,sigma_eps,type_x="continuous",
                                              n_strong_instruments)  
#        train_index <- sample(seq_len(nrow(data)), size = nrow(data)*0.7)
#        train_data  <- data[train_index, ]
#        test_data   <- data[-train_index, ]

        # OLS model
        ols <- lm(y ~ x, data = data)
        beta_1_ols_cc2[i,1] = coefficients(ols)[2] 
        beta_1_ols_cc2[i,2] = sqrt(diag(vcov(ols)))[2] 
        bias_beta_1_ols_cc2[i] = beta[2]- beta_1_ols_cc2[i,1]
        
        # IV model
        y <- data$y
        x <- data$x
        z <- as.matrix(data[,-which(colnames(data)%in% c("y","x"))])
        iv <- ivreg(y ~ x| z)
        beta_1_iv_cc2[i,1] = coefficients(iv)[2] 
        beta_1_iv_cc2[i,2] = sqrt(diag(vcov(iv)))[2] 
        bias_beta_1_iv_cc2[i] = beta[2]- beta_1_iv_cc2[i,1]
        
        # Estimation of the first stage with OLS
        ols <- lm(x ~ .-y,data=train_data)
        train_data$pred_ols_x = predict(ols)
        test_data$pred_ols_x   = predict(ols,newdata=test_data)
        mse_ols_cc2[i,1] = mean((train_data$x - train_data$pred_ols_x)^2) 
        mse_ols_cc2[i,2] = mean((test_data$x - test_data$pred_ols_x)^2)   
        # Estimation of the second stage with OLS
        tsls_train <- lm(y ~ pred_ols_x,data=train_data)
        tsls_test  <- lm(y ~ pred_ols_x,data=test_data)
        beta_1_tsls_cc2[i,1] = coefficients(tsls_train)[2]
        beta_1_tsls_cc2[i,2] = sqrt(diag(vcov(tsls_train)))[2]
        beta_1_tsls_cc2[i,3] = coefficients(tsls_test)[2]
        beta_1_tsls_cc2[i,4] = sqrt(diag(vcov(tsls_test)))[2]
        bias_beta_1_tsls_cc2[i,1] = beta[2]- beta_1_tsls_cc2[i,1]
        bias_beta_1_tsls_cc2[i,2] = beta[2]- beta_1_tsls_cc2[i,3]        
        
        # Estimation of the first stage with random forest
        rf <- randomForest(x ~ .-y,data=train_data,ntrees=200)
        train_data$predict_rf_x <- predict(rf) # Predict on out-of-bag training samples
        test_data$predict_rf_x  <- predict(rf,newdata = test_data)
        mse_rf_cc2[i,1] = mean((train_data$x - train_data$predict_rf_x)^2)
        mse_rf_cc2[i,2] = mean((test_data$x - test_data$predict_rf_x)^2)
        # Estimation of the second stage with OLS
        tsls_rf_train <- lm(y ~ predict_rf_x,data=train_data)
        tsls_rf_test  <- lm(y ~ predict_rf_x,data=test_data)
        beta_1_rf_cc2[i,1] = coefficients(tsls_rf_train)[2]
        beta_1_rf_cc2[i,2] = sqrt(diag(vcov(tsls_rf_train)))[2]
        beta_1_rf_cc2[i,3] = coefficients(tsls_rf_test)[2]
        beta_1_rf_cc2[i,4] = sqrt(diag(vcov(tsls_rf_test)))[2]
        bias_beta_1_rf_cc2[i,1] = beta[2]- beta_1_rf_cc2[i,1]
        bias_beta_1_rf_cc2[i,2] = beta[2]- beta_1_rf_cc2[i,3]
        
        # Estimation of the first stage with boost
        boost <- gbm(x ~ .-y,data=train_data,distribution="gaussian",
                     n.trees=200,interaction.depth = 4)
        train_data$predict_boost_x <- predict(boost) # Predict on out-of-bag training samples
        test_data$predict_boost_x  <- predict(boost,newdata = test_data)
        mse_boost_cc2[i,1] = mean((train_data$x - train_data$predict_boost_x)^2)
        mse_boost_cc2[i,2] = mean((test_data$x - test_data$predict_boost_x)^2)
        # Estimation of the second stage with OLS
        tsls_boost_train <- lm(y ~ predict_boost_x,data=train_data)
        tsls_boost_test  <- lm(y ~ predict_boost_x,data=test_data)
        beta_1_boost_cc2[i,1] = coefficients(tsls_boost_train)[2]
        beta_1_boost_cc2[i,2] = sqrt(diag(vcov(tsls_boost_train)))[2]
        beta_1_boost_cc2[i,3] = coefficients(tsls_boost_test)[2]
        beta_1_boost_cc2[i,4] = sqrt(diag(vcov(tsls_boost_test)))[2]
        bias_beta_1_boost_cc2[i,1] = beta[2]- beta_1_boost_cc2[i,1]
        bias_beta_1_boost_cc2[i,2] = beta[2]- beta_1_boost_cc2[i,3] 
}

mse_cc2 <- matrix(c(mean(mse_ols_cc2[,1]),mean(mse_ols_cc2[,2]),
                    mean(mse_rf_cc2[,1]),mean(mse_rf_cc2[,2]),
                    mean(mse_boost_cc2[,1]),mean(mse_boost_cc2[,2])),
                  3,2,byrow=TRUE)
colnames(mse_cc2) <- c("Train","Test")
rownames(mse_cc2) <- c("OLS","RF","Boost")
mse_cc2

beta_1_cc2 <- matrix(c(mean(beta_1_ols_cc2[,1]),mean(beta_1_ols_cc2[,2]),mean(bias_beta_1_ols_cc2),
                       mean(beta_1_iv_cc2[,1]),mean(beta_1_iv_cc2[,2]),mean(bias_beta_1_iv_cc2),
                       mean(beta_1_tsls_cc2[,1]),mean(beta_1_tsls_cc2[,2]),mean(bias_beta_1_tsls_cc2[,1]),
                       mean(beta_1_tsls_cc2[,3]),mean(beta_1_tsls_cc2[,4]),mean(bias_beta_1_tsls_cc2[,2]),
                       mean(beta_1_rf_cc2[,1]),mean(beta_1_rf_cc2[,2]),mean(bias_beta_1_rf_cc2[,1]),
                       mean(beta_1_rf_cc2[,3]),mean(beta_1_rf_cc2[,4]),mean(bias_beta_1_rf_cc2[,2]),
                       mean(beta_1_boost_cc2[,1]),mean(beta_1_boost_cc2[,2]),mean(bias_beta_1_boost_cc2[,1]),
                       mean(beta_1_boost_cc2[,3]),mean(beta_1_boost_cc2[,4]),mean(bias_beta_1_boost_cc2[,2])),
                     8,3,byrow=TRUE)
colnames(beta_1_cc2) <- c("Coefficient","Std. dev.","Bias") 
rownames(beta_1_cc2) <- c("OLS","IV","TSLS train","TSLS test","TSLS-RF train","TSLS-RF test",
                          "TSLS-Boost train","TSLS-Boost test")
beta_1_cc2

# Case 3: Many week instruments no covariates ----------------------------------

n_simulations <- 100
n_observations <- 1000
#sample_size <- 1500
n_strong_instruments <- 0
mean_eps <- c(0,0)
sigma_eps <- 0.3
#sigma_eps <- 0.5
beta_0 <- -0-90
beta_1 <- 0.75
beta <- c(beta_0,beta_1)
alpha_0 <- 0.3
alpha_1 <- 0.3
# alpha_1 <- 0.5
alpha <- c(alpha_0,alpha_1)

mse_ols_cc3   <- matrix(NA,n_simulations,2)
mse_rf_cc3    <- matrix(NA,n_simulations,2)
mse_boost_cc3 <- matrix(NA,n_simulations,2)

beta_1_ols_cc3   <- matrix(NA,n_simulations,4)
beta_1_iv_cc3    <- matrix(NA,n_simulations,4)
beta_1_tsls_cc3  <- matrix(NA,n_simulations,4)
beta_1_rf_cc3    <- matrix(NA,n_simulations,4)
beta_1_boost_cc3 <- matrix(NA,n_simulations,4)

bias_beta_1_ols_cc3   <- rep(NA,n_simulations)
bias_beta_1_iv_cc3    <- rep(NA,n_simulations)
bias_beta_1_tsls_cc3  <- matrix(NA,n_simulations,2)
bias_beta_1_rf_cc3    <- matrix(NA,n_simulations,2)
bias_beta_1_boost_cc3 <- matrix(NA,n_simulations,2)

for (i in 1:n_simulations){        

        train_data <- dgp_many_cont_strong_iv(n_observations,beta,alpha,
                                              mean_eps,sigma_eps,type_x="continuous",
                                              n_strong_instruments=0)  
        test_data <- dgp_many_cont_strong_iv(n_observations,beta,alpha,
                                              mean_eps,sigma_eps,type_x="continuous",
                                              n_strong_instruments=0)  
#       train_index <- sample(seq_len(nrow(data)), size = nrow(data)*0.7)
#       train_data  <- data[train_index, ]
#       test_data   <- data[-train_index, ]
        
        # OLS model
        ols <- lm(y ~ x, data = test_data)
        beta_1_ols_cc3[i,1] = coefficients(ols)[2] 
        beta_1_ols_cc3[i,2] = sqrt(diag(vcov(ols)))[2] 
        bias_beta_1_ols_cc3[i] = beta[2]- beta_1_ols_cc3[i,1]
        
        # IV model
        y <- test_data$y
        x <- test_data$x
        z <- as.matrix(test_data[,-which(colnames(test_data)%in% c("y","x"))])
        iv <- ivreg(y ~ x| z)
        beta_1_iv_cc3[i,1] = coefficients(iv)[2] 
        beta_1_iv_cc3[i,2] = sqrt(diag(vcov(iv)))[2] 
        bias_beta_1_iv_cc3[i] = beta[2]- beta_1_iv_cc3[i,1]
        
        # Estimation of the first stage with OLS
        ols <- lm(x ~ .-y,data=train_data)
        train_data$pred_ols_x = predict(ols)
        test_data$pred_ols_x   = predict(ols,newdata=test_data)
        mse_ols_cc3[i,1] = mean((train_data$x - train_data$pred_ols_x)^2) 
        mse_ols_cc3[i,2] = mean((test_data$x - test_data$pred_ols_x)^2)  
        # Estimation of the second stage with OLS
        tsls_train <- lm(y ~ pred_ols_x,data=train_data)
        tsls_test  <- lm(y ~ pred_ols_x,data=test_data)
        beta_1_tsls_cc3[i,1] = coefficients(tsls_train)[2]
        beta_1_tsls_cc3[i,2] = sqrt(diag(vcov(tsls_train)))[2]
        beta_1_tsls_cc3[i,3] = coefficients(tsls_test)[2]
        beta_1_tsls_cc3[i,4] = sqrt(diag(vcov(tsls_test)))[2]
        bias_beta_1_tsls_cc3[i,1] = beta[2]- beta_1_tsls_cc3[i,1]
        bias_beta_1_tsls_cc3[i,2] = beta[2]- beta_1_tsls_cc3[i,3]        
        
        # Estimation of the first stage with random forest
        rf <- randomForest(x ~ .-y,data=train_data,ntrees=200)
        train_data$predict_rf_x <- predict(rf) # Predict on out-of-bag training samples
        test_data$predict_rf_x  <- predict(rf,newdata = test_data)
        mse_rf_cc3[i,1] = mean((train_data$x - train_data$predict_rf_x)^2)
        mse_rf_cc3[i,2] = mean((test_data$x - test_data$predict_rf_x)^2)
        # Estimation of the second stage with OLS
        tsls_rf_train <- lm(y ~ predict_rf_x,data=train_data)
        tsls_rf_test  <- lm(y ~ predict_rf_x,data=test_data)
        beta_1_rf_cc3[i,1] = coefficients(tsls_rf_train)[2]
        beta_1_rf_cc3[i,2] = sqrt(diag(vcov(tsls_rf_train)))[2]
        beta_1_rf_cc3[i,3] = coefficients(tsls_rf_test)[2]
        beta_1_rf_cc3[i,4] = sqrt(diag(vcov(tsls_rf_test)))[2]
        bias_beta_1_rf_cc3[i,1] = beta[2]- beta_1_rf_cc3[i,1]
        bias_beta_1_rf_cc3[i,2] = beta[2]- beta_1_rf_cc3[i,3]
        
        # Estimation of the first stage with boost
        boost <- gbm(x ~ .-y,data=train_data,distribution="gaussian",
                     n.trees=200,interaction.depth = 4)
        train_data$predict_boost_x <- predict(boost) # Predict on out-of-bag training samples
        test_data$predict_boost_x  <- predict(boost,newdata = test_data)
        mse_boost_cc3[i,1] = mean((train_data$x - train_data$predict_boost_x)^2)
        mse_boost_cc3[i,2] = mean((test_data$x - test_data$predict_boost_x)^2)
        # Estimation of the second stage with OLS
        tsls_boost_train <- lm(y ~ predict_boost_x,data=train_data)
        tsls_boost_test  <- lm(y ~ predict_boost_x,data=test_data)
        beta_1_boost_cc3[i,1] = coefficients(tsls_boost_train)[2]
        beta_1_boost_cc3[i,2] = sqrt(diag(vcov(tsls_boost_train)))[2]
        beta_1_boost_cc3[i,3] = coefficients(tsls_boost_test)[2]
        beta_1_boost_cc3[i,4] = sqrt(diag(vcov(tsls_boost_test)))[2]
        bias_beta_1_boost_cc3[i,1] = beta[2]- beta_1_boost_cc3[i,1]
        bias_beta_1_boost_cc3[i,2] = beta[2]- beta_1_boost_cc3[i,3] 
}

mse_cc3 <- matrix(c(mean(mse_ols_cc3[,1]),mean(mse_ols_cc3[,2]),
                    mean(mse_rf_cc3[,1]),mean(mse_rf_cc3[,2]),
                    mean(mse_boost_cc3[,1]),mean(mse_boost_cc3[,2])),
                  3,2,byrow=TRUE)
colnames(mse_cc3) <- c("Train","Test")
rownames(mse_cc3) <- c("OLS","RF","Boost")
mse_cc3

beta_1_cc3 <- matrix(c(mean(beta_1_ols_cc3[,1]),mean(beta_1_ols_cc3[,2]),mean(bias_beta_1_ols_cc3),
                       mean(beta_1_iv_cc3[,1]),mean(beta_1_iv_cc3[,2]),mean(bias_beta_1_iv_cc3),
                       mean(beta_1_tsls_cc3[,1]),mean(beta_1_tsls_cc3[,2]),mean(bias_beta_1_tsls_cc3[,1]),
                       mean(beta_1_tsls_cc3[,3]),mean(beta_1_tsls_cc3[,4]),mean(bias_beta_1_tsls_cc3[,2]),
                       mean(beta_1_rf_cc3[,1]),mean(beta_1_rf_cc3[,2]),mean(bias_beta_1_rf_cc3[,1]),
                       mean(beta_1_rf_cc3[,3]),mean(beta_1_rf_cc3[,4]),mean(bias_beta_1_rf_cc3[,2]),
                       mean(beta_1_boost_cc3[,1]),mean(beta_1_boost_cc3[,2]),mean(bias_beta_1_boost_cc3[,1]),
                       mean(beta_1_boost_cc3[,3]),mean(beta_1_boost_cc3[,4]),mean(bias_beta_1_boost_cc3[,2])),
                     8,3,byrow=TRUE)
colnames(beta_1_cc3) <- c("Coefficient","Std. dev.","Bias") 
rownames(beta_1_cc3) <- c("OLS","IV","TSLS train","TSLS test","TSLS-RF train","TSLS-RF test",
                          "TSLS-Boost train","TSLS-Boost test")
beta_1_cc3


################### BINARY ENDOGENOUS EXPLANATORY VARIABLE #####################

# Case 1: One instrument no covariates -----------------------------------------

n_simulations <- 100
n_observations <- 1000
#sample_size <- 1500
mean_eps <- c(0,0)
sigma_eps <- 0.3
#sigma_eps <- 0.5
beta_0 <- -0-90
beta_1 <- 0.75
beta <- c(beta_0,beta_1)
alpha_0 <- 0.3
alpha_1 <- 0.3
# alpha_1 <- 0.5
alpha <- c(alpha_0,alpha_1)

mse_logit_bc1 <- matrix(NA,n_simulations,2)
mse_rf_bc1    <- matrix(NA,n_simulations,2)
mse_boost_bc1 <- matrix(NA,n_simulations,2)
#mse_xgboost_bc1 <- matrix(NA,n_simulations,2)

beta_1_ols_bc1   <- matrix(NA,n_simulations,4)
beta_1_tsls_bc1  <- matrix(NA,n_simulations,4)
beta_1_rf_bc1    <- matrix(NA,n_simulations,4)
beta_1_boost_bc1 <- matrix(NA,n_simulations,4)
#beta_1_xgboost_bc1 <- matrix(NA,n_simulations,4)

bias_beta_1_ols_bc1   <- rep(NA,n_simulations)
bias_beta_1_tsls_bc1  <- matrix(NA,n_simulations,2)
bias_beta_1_rf_bc1    <- matrix(NA,n_simulations,2)
bias_beta_1_boost_bc1 <- matrix(NA,n_simulations,2)
#bias_beta_1_xgboost_bc1 <- matrix(NA,n_simulations,2)

for (i in 1:n_simulations){        
        
        train_data <- dgp_one_cont_iv(n_observations,beta,alpha,mean_eps,sigma_eps,type_x="binary")
        test_data <- dgp_one_cont_iv(n_observations,beta,alpha,mean_eps,sigma_eps,type_x="binary")
#       train_index <- sample(seq_len(nrow(data)), size = nrow(data)*0.7)
#       train_data  <- data[train_index, ]
#       test_data   <- data[-train_index, ]

        # OLS model
        ols <- lm(y ~ x,data=test_data)
        beta_1_ols_bc1[i,1] = coefficients(ols)[2] 
        beta_1_ols_bc1[i,2] = sqrt(diag(vcov(ols)))[2] 
        bias_beta_1_ols_bc1[i] = beta[2]- beta_1_ols_bc1[i,1]

        # Estimation of the first stage with logit
        logit <- glm(x ~ z, family = "binomial", data = train_data)
        train_data$pred_logit_x =  predict(logit,type ="response")
        test_data$pred_logit_x   = predict(logit,newdata=test_data,type ="response")
        x_hat_logit_tr <- ifelse(train_data$pred_logit_x  > 0.5,1,0)
        x_hat_logit_ts <- ifelse(test_data$pred_logit_x > 0.5,1,0)
        mse_logit_bc1[i,1] = mean(train_data$x != x_hat_logit_tr)  
        mse_logit_bc1[i,2] = mean(test_data$x != x_hat_logit_ts) 
        # Estimation of the second stage with OLS        
        tsls_train <- lm(y ~ pred_logit_x,data=train_data)
        tsls_test  <- lm(y ~ pred_logit_x,data=test_data)
        beta_1_tsls_bc1[i,1] = coefficients(tsls_train)[2]
        beta_1_tsls_bc1[i,2] = sqrt(diag(vcov(tsls_train)))[2]
        beta_1_tsls_bc1[i,3] = coefficients(tsls_test)[2]
        beta_1_tsls_bc1[i,4] = sqrt(diag(vcov(tsls_test)))[2]
        bias_beta_1_tsls_bc1[i,1] = beta[2]- beta_1_tsls_bc1[i,1]
        bias_beta_1_tsls_bc1[i,2] = beta[2]- beta_1_tsls_bc1[i,3]   
        
        # Estimation of the first stage with random forest
        rf <- randomForest(x ~ z,data=train_data,ntrees=200)
        train_data$predict_rf_x <- predict(rf,type="class") # Predict on out-of-bag training samples
        x_hat_rf_tr <- ifelse(train_data$predict_rf_x > 0.5,1,0)
        test_data$predict_rf_x <- predict(rf,type="class",newdata = test_data)
        x_hat_rf_ts <- ifelse(test_data$predict_rf_x > 0.5,1,0)
        mse_rf_bc1[i,1] = mean(train_data$x != x_hat_rf_tr) 
        mse_rf_bc1[i,1] = mean(test_data$x != x_hat_rf_ts) 
        # Estimation of the second stage with OLS                
        tsls_rf_train <- lm(y ~ predict_rf_x,data=train_data)
        tsls_rf_test  <- lm(y ~ predict_rf_x,data=test_data)
        beta_1_rf_bc1[i,1] = coefficients(tsls_rf_train)[2]
        beta_1_rf_bc1[i,2] = sqrt(diag(vcov(tsls_rf_train)))[2]
        beta_1_rf_bc1[i,3] = coefficients(tsls_rf_test)[2]
        beta_1_rf_bc1[i,4] = sqrt(diag(vcov(tsls_rf_test)))[2]
        bias_beta_1_rf_bc1[i,1] = beta[2]- beta_1_rf_bc1[i,1]
        bias_beta_1_rf_bc1[i,2] = beta[2]- beta_1_rf_bc1[i,3] 
        
        #data$x.2 <- 0
        #forest.iv = instrumental_forest(X=data$x.1,Y=data$y,W=data$x.1,Z=data$z.1, reduced.form.weight = 0, mtry = p)
        #preds.iv = predict(forest.iv, test_data, estimate.variance = TRUE)$predictions
        #tau.hat = preds.iv$predictions
        #var.hat = preds.iv$variance.estimates
        #output = data.frame(dummy, TAU=tau.hat, VAR=var.hat)
        
        # Estimation of the first stage with boost
        boost <- gbm(x ~ z,data=train_data,distribution="bernoulli",
                     n.trees=trees)
        train_data$predict_boost_x <- predict(boost) # Predict on out-of-bag training samples
        x_hat_boost_tr <- ifelse(train_data$predict_boost_x > 0.5,1,0)
        test_data$predict_boost_x <- predict(boost,newdata = test_data,n.trees = trees)
        x_hat_boost_ts <- ifelse(test_data$predict_boost_x > 0.5,1,0)
        mse_boost_bc1[i,1] = mean(train_data$x != x_hat_boost_tr) 
        mse_boost_bc1[i,1] = mean(test_data$x != x_hat_boost_ts) 
        
        # Estimation of the second stage with OLS                
        tsls_boost_train <- lm(y ~ predict_boost_x,data=train_data)
        tsls_boost_test  <- lm(y ~ predict_boost_x,data=test_data)
        beta_1_boost_bc1[i,1] = coefficients(tsls_boost_train)[2]
        beta_1_boost_bc1[i,2] = sqrt(diag(vcov(tsls_boost_train)))[2]
        beta_1_boost_bc1[i,3] = coefficients(tsls_boost_test)[2]
        beta_1_boost_bc1[i,4] = sqrt(diag(vcov(tsls_boost_test)))[2]
        bias_beta_1_boost_bc1[i,1] = beta[2]- beta_1_boost_bc1[i,1]
        bias_beta_1_boost_bc1[i,2] = beta[2]- beta_1_boost_bc1[i,3] 
}

mse_bc1 <- matrix(c(mean(mse_logit_bc1[,1]),mean(mse_logit_bc1[,2]),
                mean(mse_rf_bc1[,1]),mean(mse_rf_bc1[,2]),
                mean(mse_boost_bc1[,1]),mean(mse_boost_bc1[,2])),
                3,2,byrow=TRUE)
colnames(mse_bc1) <- c("Train","Test")
rownames(mse_bc1) <- c("Logit","RF","Boost")
mse_bc1

beta_1_bc1 <- matrix(c(mean(beta_1_ols_bc1[,1]),mean(beta_1_ols_bc1[,2]),mean(bias_beta_1_ols_bc1),
                   mean(beta_1_tsls_bc1[,1]),mean(beta_1_tsls_bc1[,2]),mean(bias_beta_1_tsls_bc1[,1]),
                   mean(beta_1_tsls_bc1[,3]),mean(beta_1_tsls_bc1[,4]),mean(bias_beta_1_tsls_bc1[,2]),
                   mean(beta_1_rf_bc1[,1]),mean(beta_1_rf_bc1[,2]),mean(bias_beta_1_rf_bc1[,1]),
                   mean(beta_1_rf_bc1[,3]),mean(beta_1_rf_bc1[,4]),mean(bias_beta_1_rf_bc1[,2]),
                   mean(beta_1_boost_bc1[,1]),mean(beta_1_boost_bc1[,2]),mean(bias_beta_1_boost_bc1[,1]),
                   mean(beta_1_boost_bc1[,3]),mean(beta_1_boost_bc1[,4]),mean(bias_beta_1_boost_bc1[,2])),
                   7,3,byrow=TRUE)
colnames(beta_1_bc1) <- c("Coefficient","Std. dev.","Bias") 
rownames(beta_1_bc1) <- c("OLS","TSLS train","TSLS test","TSLS-RF train","TSLS-RF test",
                          "TSLS-Boost train","TSLS-Boost test")

# Case 2: Many instruments with strong sparcity no covariates ------------------

n_simulations <- 100
n_observations <- 1000
#sample_size <- 1500
n_strong_instruments <- 25
#n_strong_instruments <- 50
mean_eps <- c(0,0)
sigma_eps <- 0.3
#sigma_eps <- 0.5
beta_0 <- -0-90
beta_1 <- 0.75
beta <- c(beta_0,beta_1)
alpha_0 <- 0.3
alpha_1 <- 3/n_strong_instruments
# alpha_1 <- 6/n_strong_instruments
alpha <- c(alpha_0,alpha_1)

mse_logit_bc1 <- matrix(NA,n_simulations,2)
mse_rf_bc1    <- matrix(NA,n_simulations,2)
mse_boost_bc1 <- matrix(NA,n_simulations,2)
#mse_xgboost_bc1 <- matrix(NA,n_simulations,2)

beta_1_ols_bc1   <- matrix(NA,n_simulations,4)
beta_1_tsls_bc1  <- matrix(NA,n_simulations,4)
beta_1_rf_bc1    <- matrix(NA,n_simulations,4)
beta_1_boost_bc1 <- matrix(NA,n_simulations,4)
#beta_1_xgboost_bc1 <- matrix(NA,n_simulations,4)

bias_beta_1_ols_bc1   <- rep(NA,n_simulations)
bias_beta_1_tsls_bc1  <- matrix(NA,n_simulations,2)
bias_beta_1_rf_bc1    <- matrix(NA,n_simulations,2)
bias_beta_1_boost_bc1 <- matrix(NA,n_simulations,2)
#bias_beta_1_xgboost_bc1 <- matrix(NA,n_simulations,2)

for (i in 1:n_simulations){        
        
        train_data <- dgp_many_cont_strong_iv(n_observations,beta,alpha,
                                              mean_eps,sigma_eps,type_x="binary",
                                              n_strong_instruments)  
        test_data <- dgp_many_cont_strong_iv(n_observations,beta,alpha,
                                             mean_eps,sigma_eps,type_x="binary",
                                             n_strong_instruments) 
        #       train_index <- sample(seq_len(nrow(data)), size = nrow(data)*0.7)
        #       train_data  <- data[train_index, ]
        #       test_data   <- data[-train_index, ]
        
        # OLS model
        ols <- lm(y ~ x,data=test_data)
        beta_1_ols_bc2[i,1] = coefficients(ols)[2] 
        beta_1_ols_bc2[i,2] = sqrt(diag(vcov(ols)))[2] 
        bias_beta_1_ols_bc2[i] = beta[2]- beta_1_ols_bc2[i,1]
        
        # Estimation of the first stage with logit
        logit <- glm(x ~ . -y, family = "binomial", data = train_data)
        train_data$pred_logit_x =  predict(logit,type ="response")
        test_data$pred_logit_x   = predict(logit,newdata=test_data,type ="response")
        x_hat_logit_tr <- ifelse(train_data$pred_logit_x  > 0.5,1,0)
        x_hat_logit_ts <- ifelse(test_data$pred_logit_x > 0.5,1,0)
        mse_logit_bc2[i,1] = mean(train_data$x != x_hat_logit_tr)  
        mse_logit_bc2[i,2] = mean(test_data$x != x_hat_logit_ts) 
        # Estimation of the second stage with OLS        
        tsls_train <- lm(y ~ pred_logit_x,data=train_data)
        tsls_test  <- lm(y ~ pred_logit_x,data=test_data)
        beta_1_tsls_bc2[i,1] = coefficients(tsls_train)[2]
        beta_1_tsls_bc2[i,2] = sqrt(diag(vcov(tsls_train)))[2]
        beta_1_tsls_bc2[i,3] = coefficients(tsls_test)[2]
        beta_1_tsls_bc2[i,4] = sqrt(diag(vcov(tsls_test)))[2]
        bias_beta_1_tsls_bc2[i,1] = beta[2]- beta_1_tsls_bc2[i,1]
        bias_beta_1_tsls_bc2[i,2] = beta[2]- beta_1_tsls_bc2[i,3]   
        
        # Estimation of the first stage with random forest
        rf <- randomForest(x ~ .-y,data=train_data,ntrees=200)
        train_data$predict_rf_x <- predict(rf,type="class") # Predict on out-of-bag training samples
        x_hat_rf_tr <- ifelse(train_data$predict_rf_x > 0.5,1,0)
        test_data$predict_rf_x <- predict(rf,type="class",newdata = test_data)
        x_hat_rf_ts <- ifelse(test_data$predict_rf_x > 0.5,1,0)
        mse_rf_bc2[i,1] = mean(train_data$x != x_hat_rf_tr) 
        mse_rf_bc2[i,1] = mean(test_data$x != x_hat_rf_ts) 
        # Estimation of the second stage with OLS                
        tsls_rf_train <- lm(y ~ predict_rf_x,data=train_data)
        tsls_rf_test  <- lm(y ~ predict_rf_x,data=test_data)
        beta_1_rf_bc2[i,1] = coefficients(tsls_rf_train)[2]
        beta_1_rf_bc2[i,2] = sqrt(diag(vcov(tsls_rf_train)))[2]
        beta_1_rf_bc2[i,3] = coefficients(tsls_rf_test)[2]
        beta_1_rf_bc2[i,4] = sqrt(diag(vcov(tsls_rf_test)))[2]
        bias_beta_1_rf_bc2[i,1] = beta[2]- beta_1_rf_bc2[i,1]
        bias_beta_1_rf_bc2[i,2] = beta[2]- beta_1_rf_bc2[i,3] 
        
        #data$x.2 <- 0
        #forest.iv = instrumental_forest(X=data$x.1,Y=data$y,W=data$x.1,Z=data$z.1, reduced.form.weight = 0, mtry = p)
        #preds.iv = predict(forest.iv, test_data, estimate.variance = TRUE)$predictions
        #tau.hat = preds.iv$predictions
        #var.hat = preds.iv$variance.estimates
        #output = data.frame(dummy, TAU=tau.hat, VAR=var.hat)
        
        # Estimation of the first stage with boost
        boost <- gbm(x ~ .-y,data=train_data,distribution="bernoulli",
                     n.trees=trees)
        train_data$predict_boost_x <- predict(boost) # Predict on out-of-bag training samples
        x_hat_boost_tr <- ifelse(train_data$predict_boost_x > 0.5,1,0)
        test_data$predict_boost_x <- predict(boost,newdata = test_data,n.trees = trees)
        x_hat_boost_ts <- ifelse(test_data$predict_boost_x > 0.5,1,0)
        mse_boost_bc2[i,1] = mean(train_data$x != x_hat_boost_tr) 
        mse_boost_bc2[i,1] = mean(test_data$x != x_hat_boost_ts) 
        # Estimation of the second stage with OLS                
        tsls_boost_train <- lm(y ~ predict_boost_x,data=train_data)
        tsls_boost_test  <- lm(y ~ predict_boost_x,data=test_data)
        beta_1_boost_bc2[i,1] = coefficients(tsls_boost_train)[2]
        beta_1_boost_bc2[i,2] = sqrt(diag(vcov(tsls_boost_train)))[2]
        beta_1_boost_bc2[i,3] = coefficients(tsls_boost_test)[2]
        beta_1_boost_bc2[i,4] = sqrt(diag(vcov(tsls_boost_test)))[2]
        bias_beta_1_boost_bc2[i,1] = beta[2]- beta_1_boost_bc2[i,1]
        bias_beta_1_boost_bc2[i,2] = beta[2]- beta_1_boost_bc2[i,3] 
}

mse_bc2 <- matrix(c(mean(mse_logit_bc2[,1]),mean(mse_logit_bc2[,2]),
                    mean(mse_rf_bc2[,1]),mean(mse_rf_bc2[,2]),
                    mean(mse_boost_bc2[,1]),mean(mse_boost_bc2[,2])),
                  3,2,byrow=TRUE)
colnames(mse_bc2) <- c("Train","Test")
rownames(mse_bc2) <- c("Logit","RF","Boost")
mse_bc2

beta_1_bc2 <- matrix(c(mean(beta_1_ols_bc2[,1]),mean(beta_1_ols_bc2[,2]),mean(bias_beta_1_ols_bc2),
                       mean(beta_1_tsls_bc2[,1]),mean(beta_1_tsls_bc2[,2]),mean(bias_beta_1_tsls_bc2[,1]),
                       mean(beta_1_tsls_bc2[,3]),mean(beta_1_tsls_bc2[,4]),mean(bias_beta_1_tsls_bc2[,2]),
                       mean(beta_1_rf_bc2[,1]),mean(beta_1_rf_bc2[,2]),mean(bias_beta_1_rf_bc2[,1]),
                       mean(beta_1_rf_bc2[,3]),mean(beta_1_rf_bc2[,4]),mean(bias_beta_1_rf_bc2[,2]),
                       mean(beta_1_boost_bc2[,1]),mean(beta_1_boost_bc2[,2]),mean(bias_beta_1_boost_bc2[,1]),
                       mean(beta_1_boost_bc2[,3]),mean(beta_1_boost_bc2[,4]),mean(bias_beta_1_boost_bc2[,2])),
                     7,3,byrow=TRUE)
colnames(beta_1_bc2) <- c("Coefficient","Std. dev.","Bias") 
rownames(beta_1_bc2) <- c("OLS","TSLS train","TSLS test","RF train","RF test",
                          "Boost train","Boost test")

# Case 3: Many week instruments no covariates ----------------------------------

n_simulations <- 100
n_observations <- 1000
n_strong_instruments <- 0
#sample_size <- 1500
mean_eps <- c(0,0)
sigma_eps <- 0.3
#sigma_eps <- 0.5
beta_0 <- -0-90
beta_1 <- 0.75
beta <- c(beta_0,beta_1)
alpha_0 <- 0.3
alpha_1 <- 0.3
# alpha_1 <- 0.5
alpha <- c(alpha_0,alpha_1)





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

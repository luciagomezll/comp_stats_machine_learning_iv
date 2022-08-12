#################################################################################
### SIMULATION OF THE STUDY OF RETURNS TO COLLEGE
#################################################################################

rm(list=ls())

set.seed(12345)

# Load libraries
if(!require(haven)) {install.packages("haven")}
if(!require(mvtnorm)) {install.packages("mvtnorm")}
if(!require(MASS)) {install.packages("MASS")}
if(!require(ivreg)) install.packages("ivreg")
if(!require(caret)) install.packages("caret")
if(!require(randomForest)) install.packages("randomForest")
if(!require(grf)) install.packages("grf")
if(!require(gbm)) install.packages("gbm")
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
        trees <- 400
        rf <- randomForest(x.1 ~ z.1,data=train_data,ntrees=trees)
        train_data$predict_rf_x.1 <- predict(rf) # Predict on out-of-bag training samples
        test_data$predict_rf_x.1 <- predict(rf,newdata = test_data)
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

# Instruments not correlated to each other

# Instruments correlated to each other

# B.1 Many instruments with strong sparcity with covariates
#===============================================================================

# Instruments not correlated to each other

# Instruments correlated to each other

# C. Many week with strong sparcity no covariates
#===============================================================================

# Instruments not correlated to each other





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
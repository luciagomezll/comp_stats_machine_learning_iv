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

summary(local_level_dta) <- read_dta('localvariables.dta')
basic_level_dta <- read_dta('basicvariables.dta')

obs <- 1000
obs_local <- obs/10

# parameters for simulations
dgp <- function(obs,){
        # Create sample
        
        # Instruments
        
        # Propensity score
        
        # Group shares
        
        # Assign groups
        
        # Treatment (vincular con los instrumentos)
        S <- rbinom(obs,1,mean(basic_level_dta$state))
        
        # Covariates
        # basic variables
        x0 <- rep(1,obs)
        x1_cafqt <- rnorm(obs,mean(basic_level_dta$cafqt),sd(basic_level_dta$cafqt)) 
        x2_msch  <- rnorm(obs,mean(basic_level_dta$mhgc),sd(basic_level_dta$mhgc))
        x3_numsibs <- rnorm(obs,mean(basic_level_dta$numsibs),sd(basic_level_dta$numsibs)) # It seems to follow and F-distribution
        x4_urban14 <- rbinom(obs,1,mean(basic_level_dta$urban14))
        prob.x_d <- as.vector(colMeans(basic_level_dta[,c(6:12)])) 
        x5_d  <- sample(c(57:63), obs, replace=TRUE, prob=prob.x_d)
        x6_exp <- rnorm(obs,mean(basic_level_dta$exp),sd(basic_level_dta$exp))
        x7_expsq <- sq(x6_exp) 
        eps <- rnorm(obs,0,1)
        
        # Create school ID to match students with schools
        obs_per_local <- sample(1:50, obs_local, replace=TRUE)
        prop.obs_local <- obs_per_local/sum(obs_per_local)
        local.id <- sample(1:obs_local,obs,replace=TRUE, prob= prop.obs_local)
        
        # Local level variables
        data.x <- summary(local_level_dta) %>%                                        # Create ID by group
                 group_by(local.id) %>%
                 mutate(rnorm(obs,))
        
        
        # local variables
        avlocwage17 <- rnorm(obs_local,mean(local_level_dta$lavlocwage17),sd(local_level_dta$lavlocwage17))
        avurate <- rnorm(obs_local,mean(local_level_dta$avurate),sd(local_level_dta$avurate))
        lwage_14 <- rnorm(obs_local,mean(local_level_dta$lwage5),sd(local_level_dta$lwage5))
        lurate_14 <- rnorm(obs_local,mean(local_level_dta$lurate),sd(local_level_dta$lurate))
        lwage5_17 <- rnorm(obs_local,mean(local_level_dta$lwage5_17),sd(local_level_dta$lwage5_17))
        lurate_17 <- rnorm(obs_local,mean(local_level_dta$lurate_17),sd(local_level_dta$lurate_17))
        tuit4c <- rnorm(obs_local,mean(local_level_dta$tuit4c),sd(local_level_dta$tuit4c))
        pub4 <- rbinom(obs,1,mean(basic_level_dta$pub4))
        
        # Potential outcomes
        
        }









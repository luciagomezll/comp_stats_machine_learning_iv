#################################################################################
### SIMULATION OF THE STUDY OF RETURNS TO COLLEGE
#################################################################################

rm(list=ls())

# Load libraries
if (!require("haven")) {install.packages("haven")}
library(haven)

# Set up the working directory
# setwd("F:/MASTER_BONN/SECOND SEMESTER/comp_stats/project/135041-V1/data/") # Change the user
# shapefiles <- paste0(main, "clean/")

# a. Generate the simulation of the data set

dgp <- function(){
        # Create sample
        # Instruments
        # Propensity score
        # Group shares
        # Assign groups
        # Treatment
        # Covariates
        # Potential outcomes
        
#        local_level_data <- read_dta(simulated_data.dta)
        basic_level_data <- read_dta(localvariables.dta)
        local_level_data <- read_dta(basicvariables.dta)
}










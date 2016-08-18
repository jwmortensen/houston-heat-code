# Runs initial code to prepare the model to be fit ============================
source("../code/setup.R")

# Initialize vars for MCMC ====================================================
mod_name <- "fit11"
n_draws <- 10000
thin_factor <- 10
burn <- 25000
X <- model.matrix(~1 + (HI_MAX + over65PCT + under5PCT + BlackPCT + HispanicPCT)^2,
 intercept_df)

# Fits model using MCMC and saves draws in ./data/mod_name.Rdata ==============
source("../code/fit_model.R")
library(Rcpp)
library(RcppArmadillo)
library(FNN)
library(fields)
library(gtools)

load("../data/Spatial911PtPtrn.RData") # calls, hrdlas.grid, pp.grid, houston, kp.gp
load("../data/AllTempData.RData") # temp.data.nomiss
load("../data/RevisedInterceptData.RData") # intercept.df
sourceCpp("../code/estimateParametersLowerMu.cpp")
death <- read.csv("../data/cleanDeathData.csv")

# Functions
fac_to_num <- function(x) {
  as.numeric(as.character(x))
}

# Define Global Vars
n_predlocs <- nrow(pp.grid)
n_blocks <- 21 # Factors evenly into n_predlocs
n_per_block <- n_predlocs / n_blocks
block_pts <- pp.grid[seq.int(from = 1, to = n_predlocs, length = n_blocks), ]

close_pts <- vector("list", length = n_blocks)
cp_ind <- vector("list", length = n_blocks)
pp_copy <- pp.grid

for (i in 1:n_blocks) {
  nn <- get.knnx(pp_copy, block_pts[i, ], n_per_block)
  close_pts[[i]] <- pp_copy[nn$nn.index, ]
  cp_ind[[i]] <- which(paste(pp.grid$Latitude, pp.grid$Longitude)
                       %in% paste(close_pts[[i]]$Latitude, close_pts[[i]]$Longitude))
  pp_copy <- pp_copy[-nn$nn.index, ]
}

# Get nearest prediction locations for each observed call
call_locs <- calls[, 1:2]
call_nn <- get.knnx(pp.grid, call_locs, 1)
call_ind <- call_nn$nn.index
Nk_911 <- numeric(n_predlocs)
for (i in call_ind) Nk_911[i] <- Nk_911[i] + 1

# Get nearest prediction locations for each death
death_locs <- data.frame(lat = death$D_LAT, lon = death$D_LONG)
death_nn <- get.knnx(pp.grid, death_locs, 1)
death_ind <- death_nn$nn.index
Nk_death <- numeric(n_predlocs)
for (i in death_ind) Nk_death[i] <- Nk_death[i] + 1

# Total number of events in each grid cell
Nk_mu <- Nk_911 + Nk_death

# Calculate expected number of events if only population matters
E_911 <- sum(Nk_911) * intercept_df$TotalPopulation / sum(intercept_df$TotalPopulation)
E_death <- sum(Nk_death) * intercept_df$TotalPopulation / sum(intercept_df$TotalPopulation)

# Calculate distance between elements for use in GP
spat_dist <- rdist(pp.grid)

# Fix nu
nu <- 3.5

# Fix phi based on selection of previous model runs
lambda_phi <- 500
lambda_matern <- Matern(spat_dist, alpha = lambda_phi, nu = nu)
chol_matern <- chol(lambda_matern)
lambda_invmat <- solve(lambda_matern)

# Initialize AMCMC
amcmc_911 <- vector("list", length = n_blocks)
amcmc_death <- vector("list", length = n_blocks)
amcmc_mu <- vector("list", length = n_blocks)

for (i in 1:n_blocks) {
  cpl <- length(cp_ind[[i]])
  amcmc_911[[i]]$mn <- amcmc_death[[i]]$mn <- 
    amcmc_mu[[i]]$mn <- matrix(0, ncol = 1, nrow = cpl)
  amcmc_911[[i]]$var <- amcmc_death[[i]]$var <-
    amcmc_mu[[i]]$var <- matrix(0, ncol = cpl, nrow = cpl)
}
amcmc_it <- 100

# Function used to calculate DIC
D <- function(lstar_911, lstar_death) {
  -2 * (LogLike(lstar_911, Nk_911, E_911) + 
          LogLike(lstar_death, Nk_death, E_death))
}

DIC <- function(d_vals, lstar) {
  d_bar <- mean(d_vals)
  l911_bar <- apply(lstar$`911`, 2, mean)
  ldeath_bar <- apply(lstar$death, 2, mean)
  d_theta_bar <- D(l911_bar, ldeath_bar)
  2*d_bar - d_theta_bar
}

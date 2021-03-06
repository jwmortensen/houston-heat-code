library(Rcpp)
library(RcppArmadillo)
library(FNN)
library(fields)
library(gtools)

load("../data/Spatial911PtPtrn.RData") # calls, hrdlas.grid, pp.grid, houston, kp.gp
load("../data/AllTempData.RData") # temp.data.nomiss
load("../data/RevisedInterceptData.RData") # intercept.df
sourceCpp("./estimateParametersLowerMu.cpp")
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

############ Gibbs Sampler ######################

# Initialize vars
n_draws <- 100
thin_factor <- 1
burn <- 100
X <- model.matrix(~1 + HI_MIN + over65PCT + under5PCT + NOACPCT + BlackPCT +
                    HispanicPCT + disabledPCT + alonePCT + povertyPCT,
                  data = intercept_df)

# delta is conjugate so we draw it up front
# think of the total number of 911 calls / deaths as one draw from a poisson distribution 
delta_911 <- rgamma(n_draws, shape = nrow(calls) + 0.001, rate = 1 + 0.001)
delta_death <- rgamma(n_draws, shape = nrow(death) + 0.001, rate = 1 + 0.001)

# Conjugate draws for age probabilities
alpha_hp <- 2
count_ages_911 <- data.frame(table(round(calls$Age)))
names(count_ages_911) <- c("age", "count")
count_ages_911$age <- fac_to_num(count_ages_911$age)
count_ages_911$count <- fac_to_num(count_ages_911$count)
n_ages_911 <- merge(data.frame(age = 0:113), count_ages_911, all.x = T)
n_ages_911$count[is.na(n_ages_911$count)] <- 0
age_draws_911 <- rdirichlet(n_draws, n_ages_911$count + alpha_hp)

count_ages_death <- data.frame(table(round(death$Age_Presume)))
names(count_ages_death) <- c("age", "count")
count_ages_death$age <- fac_to_num(count_ages_death$age)
count_ages_death$count <- fac_to_num(count_ages_death$count)
n_ages_death <- merge(data.frame(age = 0:113), count_ages_death, all.x = T)
n_ages_death$count[is.na(n_ages_death$count)] <- 0
age_draws_death <- rdirichlet(n_draws, n_ages_death$count + alpha_hp)

# Conjugate draws for gender probabilities
n_gender_911 <- data.frame(table(calls$Gender))
names(n_gender_911) <- c("gender", "count")
gender_draws_911 <- rdirichlet(n_draws, n_gender_911$count + alpha_hp)

death$gender_num <- ifelse(death$Gender == "Female", 1, 0)
n_gender_death <- data.frame(table(death$gender_num))
names(n_gender_death) <- c("gender", "count")
gender_draws_death <- rdirichlet(n_draws, n_gender_death$count + alpha_hp)

# Conjugate draws for ethnicity
n_eth_911 <- data.frame(table(calls$Eth))
names(n_eth_911) <- c("race", "count")
eth_draws_911 <- rdirichlet(n_draws, n_eth_911$count + alpha_hp)

n_eth_death <- data.frame(table(death$RACETH))
names(n_eth_death) <- c("ethnicity", "count")
eth_draws_death <- rdirichlet(n_draws, n_eth_death$count + alpha_hp)

# Containers to hold draws 
lstar <- list(`911` = matrix(NA, nrow = n_draws, ncol = n_predlocs),
              death = matrix(NA, nrow = n_draws, ncol = n_predlocs),
              mu = matrix(NA, nrow = n_draws, ncol = n_predlocs))
lvar <- list(`911` = numeric(n_draws), 
             death = numeric(n_draws),
             mu = numeric(n_draws))
beta_draws <- matrix(NA, nrow = n_draws, ncol = ncol(X))
d_vals <- numeric(n_draws)

lstar_911 <- log((0.001 + Nk_911) / sum(Nk_911))
lstar_death <- log((0.001 + Nk_death) / sum(Nk_death))
lstar_mu <- log((0.001 + Nk_mu) / sum(Nk_mu))
beta <- rep(0, ncol(X))
lvar_911 <- lvar_death <- lvar_mu <- 50
lvar_a <- 2.001
lvar_b <- 1.001
lam_a <- lvar_a + n_predlocs / 2


pb <- txtProgressBar(min = 0, max = n_draws * thin_factor, style = 3)
n <- 0#ifelse(thin_factor == 1, 1, 0)

for (i in 1:(n_draws * thin_factor + burn)) {
  # Get draws for lvar_911 using the complete conditional
  lam_b <- lvar_b + 0.5 * t(lstar_911)%*%lambda_invmat%*%lstar_911
  lvar_911 <- 1 / rgamma(1, shape = lam_a, rate = lam_b)

  # Use MH to get draws for lstar_911
  for (j in 1:n_blocks) {
    # If enough iterations have been completed, begin using AMCMC
    scale_911 <- 1e-3
    if (i < amcmc_it) {
      prop_var <- diag(scale_911, n_per_block)
    } else {
      prop_var <- (2.4^2 / (n_per_block)) * (diag(scale_911, n_per_block) +
        amcmc_911[[j]]$var)
    }

    prop_lstar_911 <- lstar_911
    prop_lstar_911[cp_ind[[j]]] <- t(mvrnormC(1, lstar_911[cp_ind[[j]]], prop_var))

    log_MH <- LogLike(prop_lstar_911, Nk_911, E_911) -
      LogLike(lstar_911, Nk_911, E_911) +
      LogLambdaPrior(prop_lstar_911, lstar_mu, lvar_911, lambda_invmat) -
      LogLambdaPrior(lstar_911, lstar_mu, lvar_911, lambda_invmat)

    if (log(runif(1)) < log_MH) {
      lstar_911 <- prop_lstar_911
    }

    amcmc_911[[j]] <- amcmcUpdate(lstar_911[cp_ind[[j]]], amcmc_911[[j]]$mn, 
      amcmc_911[[j]]$var, i - 1)
  }

  # Get draws for lvar_death using the complete conditional
  lam_b <- lvar_b + 0.5 * t(lstar_death)%*%lambda_invmat%*%lstar_death
  lvar_death <- 1 / rgamma(1, shape = lam_a, rate = lam_b)

  # Use MH to get draws for lstar_death
  for (j in 1:n_blocks) {
    scale_death <- 2e-4
    if (i < amcmc_it) {
      prop_var <- diag(scale_death, n_per_block)
    } else {
      prop_var <- (2.4^2 / (n_per_block)) * (diag(scale_death, n_per_block) +
        amcmc_death[[j]]$var)
    }

    prop_lstar_death <- lstar_death
    prop_lstar_death[cp_ind[[j]]] <- t(mvrnormC(1, lstar_death[cp_ind[[j]]], 
      prop_var))

    log_MH <- LogLike(prop_lstar_death, Nk_death, E_death) -
      LogLike(lstar_death, Nk_death, E_death) +
      LogLambdaPrior(prop_lstar_death, lstar_mu, lvar_death, lambda_invmat) -
      LogLambdaPrior(lstar_death, lstar_mu, lvar_death, lambda_invmat)

    if (log(runif(1)) < log_MH) {
      lstar_death <- prop_lstar_death
    }

    amcmc_death[[j]] <- amcmcUpdate(lstar_death[cp_ind[[j]]], 
      amcmc_death[[j]]$mn, amcmc_death[[j]]$var, i - 1)
  }

  # Draw lvar_mu using the complete conditional
  lam_b <- lvar_b + 0.5 * t(lstar_mu - X%*%beta)%*%lambda_invmat%*%
    (lstar_mu - X%*%beta)
  lvar_mu <- 1 / rgamma(1, shape = lam_a, rate = lam_b)

  # Use MH to draw lstar_mu
  for (j in 1:n_blocks) {
    scale_mu <- 2e-4
    if (i < amcmc_it) {
      prop_var <- diag(scale_mu, n_per_block)
    } else {
      prop_var <- (2.4^2 / n_per_block) * (diag(scale_mu, n_per_block) + amcmc_mu[[j]]$var)
    }

    prop_lstar_mu <- lstar_mu
    prop_lstar_mu[cp_ind[[j]]] <- t(mvrnormC(1, lstar_mu[cp_ind[[j]]], 
      prop_var))

    log_MH <- LogLambdaPrior(lstar_911, prop_lstar_mu, lvar_911, lambda_invmat) - 
      LogLambdaPrior(lstar_911, lstar_mu, lvar_911, lambda_invmat) +
      LogLambdaPrior(lstar_death, prop_lstar_mu, lvar_death, lambda_invmat) -
      LogLambdaPrior(lstar_death, lstar_mu, lvar_death, lambda_invmat) +
      LogLambdaMuPrior(prop_lstar_mu, lvar_mu, lambda_invmat, X, beta) -
      LogLambdaMuPrior(lstar_mu, lvar_mu, lambda_invmat, X, beta)

    if (log(runif(1)) < log_MH) {
      lstar_mu <- prop_lstar_mu
    }

    amcmc_mu[[j]] <- amcmcUpdate(lstar_mu[cp_ind[[j]]], amcmc_mu[[j]]$mn,
      amcmc_mu[[j]]$var, i - 1)
  }

  # Complete conditional draws for beta
  s2 <- 50
  beta_var <- solve(diag(1 / s2, ncol(X)) + 
                      (1 / lvar_mu) * t(X)%*%lambda_invmat%*%X)
  beta_mean <- (1/lvar_mu)*beta_var%*%t(X)%*%lambda_invmat%*%lstar_mu
  beta <- t(mvrnormC(1, beta_mean, beta_var))

  if (i %% thin_factor == 0 & i > burn) {
    n <- n + 1
    lvar$`911`[n] <- lvar_911
    lstar$`911`[n, ] <- lstar_911
    lvar$death[n] <- lvar_death
    lstar$death[n, ] <- lstar_death
    lvar$mu[n] <- lvar_mu
    lstar$mu[n, ] <- lstar_mu
    beta_draws[n, ] <- beta
    d_vals[n] <- D(lstar_911, lstar_death)
    
  }
  setTxtProgressBar(pb, i)
}





################################################
#
# An example of an auto-logistic occupancy model
#
# Written by M. Fidino
#
################################################

# load packaages
library(runjags)
library(mcmcplots)

# Step 1. Simulate the data.

# Genearal bookkeeping
nsite <- 25
nyear <- 8
ncovar <- 3
nrepeats <- 5

# for covariates
X <- matrix(
  NA,
  ncol = ncovar,
  nrow = nsite
)

set.seed(333)
# Create covariates
X <- cbind(1,apply(X, 2, function(x) rnorm(nsite)))

# Occupancy coefficients, +1 for intercept
psi_betas <- rnorm(ncovar + 1)

# auto-logistic term
theta <- 0.75

# Detection coefficients, decreasing magnitude here
rho_betas <- rnorm(ncovar + 1, 0, 0.5)

# latent state, give same dimensions as X
z <- matrix(NA, ncol = nyear, nrow = nsite)

# Do first year
psi <- plogis(X %*% psi_betas)
z[,1] <- rbinom(nsite, 1, psi)

# And then the rest, which also uses the theta term
for(year in 2:nyear){
  psi <- plogis(X %*% psi_betas + theta * z[,year-1])
  z[,year] <- rbinom(nsite,1,psi)
}

# Add imperfect detection, make it a matrix with
#  the same dimensions as z. Then multiply by z.
rho <- matrix(
  plogis(X %*% rho_betas),
  ncol = nyear,
  nrow = nsite
) * z

# Create the observed data. Again, same dimensions.
y <- matrix(
  rbinom(
    length(rho),
    nrepeats,
    rho
  ),
  ncol = nyear,
  nrow = nsite
)

# Step 2. Fit the model
data_list <- list(
  J = nrepeats,
  nsite = nsite,
  nyear = nyear,
  X = X,
  ncovar = ncovar + 1, # for intercept
  y = y
)

my_inits <- function(chain){
  gen_list <- function(chain = chain){
    list(
      z = matrix(1, ncol = data_list$nyear, nrow = data_list$nsite),
      psi_beta = rnorm(data_list$ncovar),
      rho_beta = rnorm(data_list$ncovar),
      theta = rnorm(1),
      .RNG.name = switch(chain,
                         "1" = "base::Wichmann-Hill",
                         "2" = "base::Wichmann-Hill",
                         "3" = "base::Super-Duper",
                         "4" = "base::Mersenne-Twister",
                         "5" = "base::Wichmann-Hill",
                         "6" = "base::Marsaglia-Multicarry",
                         "7" = "base::Super-Duper",
                         "8" = "base::Mersenne-Twister"),
      .RNG.seed = sample(1:1e+06, 1)
    )
  }
  return(switch(chain,
                "1" = gen_list(chain),
                "2" = gen_list(chain),
                "3" = gen_list(chain),
                "4" = gen_list(chain),
                "5" = gen_list(chain),
                "6" = gen_list(chain),
                "7" = gen_list(chain),
                "8" = gen_list(chain)
  )
  )
}

my_mod <- runjags::run.jags(
  model = "./jags/autologistic_model.R",
  monitor = c("theta", "psi_beta", "rho_beta"),
  data = data_list,
  inits = my_inits,
  n.chains = 4,
  adapt = 1000,
  burnin = 10000,
  sample = 10000,
  method = "parallel"
)

my_sum <- summary(my_mod)

round(cbind(c(theta, psi_betas, rho_betas),my_sum),2 )

caterplot(
  my_mod,
  collapse = TRUE,
  reorder = FALSE
)

points(
  x =  c(theta, psi_betas, rho_betas),
  y = rev(1:9),
  pch = 21,
  bg = "yellow",
  cex = 1.2
)

# try with random effect oo.

my_mod2 <- runjags::run.jags(
  model = "./jags/ranef_model.R",
  monitor = c("psi_sd", "psi_beta", "rho_beta"),
  data = data_list,
  inits = my_inits,
  n.chains = 4,
  adapt = 1000,
  burnin = 10000,
  sample = 10000,
  method = "parallel"
)

caterplot(
  my_mod2,
  collapse = TRUE,
  reorder = FALSE
)

points(
  x =  c(NA, psi_betas, rho_betas),
  y = rev(1:9),
  pch = 21,
  bg = "yellow",
  cex = 1.2
)

m1_mcmc <- do.call("rbind", my_mod$mcmc)
m2_mcmc <- do.call("rbind", my_mod2$mcmc)

par(mar = c(3,7,1,1))
plot(1~1, xlim = c(-5,5), ylim = c(0,9), type = "n", bty = "l",
     xlab = "Coefficient estimate",
     ylab = "", yaxt = "n")
axis(2, rev(1:8), colnames(m1_mcmc)[-1], las = 2)
for(i in 1:8){
  # add 95% CI and median
  lines(
    x = c(
      min(m1_mcmc[,i+1]),
      max(m1_mcmc[,i+1])
    ),
    y = rep(rev(1:8)[i] + 0.15,2),
    col = "blue",
    lwd = 3
  )
  points(
    x = median(m1_mcmc[,i+1]),
    y = rev(1:8)[i] +0.15,
    pch = 21,
    bg = "blue",
    cex = 2
    )
  lines(
    x = c(
      min(m2_mcmc[,i+1]),
      max(m2_mcmc[,i+1])
    ),
    y = rep(rev(1:8)[i] - 0.15,2),
    col = "gray50",
    lwd = 3
  )
  points(
    x = median(m2_mcmc[,i+1]),
    y = rev(1:8)[i] -0.15,
    pch = 21,
    bg = "gray50",
    cex = 2
  )
 
}
legend(
  "topright",
  c("Auto-logistic model", "Random site-effect model"),
  pch = 21,
  pt.bg = c("blue", "gray50"),
  pt.cex = 2,
  bty = "n"
)

abline(v = 0, lty = 2)

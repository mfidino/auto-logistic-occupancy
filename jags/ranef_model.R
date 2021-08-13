model{
  for(site in 1:nsite){
    for(year in 1:nyear){
      #
      # latent state model, with random effect term now.
      #
      logit(psi[site,year]) <- inprod(psi_beta, X[site,]) + sre[site]
      z[site,year] ~ dbern(psi[site,year])
      #
      # data model
      #
      logit(rho[site,year]) <- inprod(rho_beta, X[site,])
      y[site,year] ~ dbin(rho[site,year] * z[site,year], J)
    }
  }
  #
  # Priors
  #
  # Intercept and slope terms
  for(covar in 1:ncovar){
    psi_beta[covar] ~ dlogis(0,1)
    rho_beta[covar] ~ dlogis(0,1)
  }
  # random effect term for site on psi
  for(site in 1:nsite){
    sre[site] ~ dnorm(0, tau)
  }
  # random effect precision
  tau ~ dgamma(1,1)
  psi_sd <- 1 / sqrt(tau)
}
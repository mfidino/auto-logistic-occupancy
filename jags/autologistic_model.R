model{
  for(site in 1:nsite){
    # This is for the first year
    #
    # latent state model
    #
    logit(psi[site,1]) <- inprod(psi_beta, X[site,])
    z[site,1] ~ dbern(psi[site,1])
    #
    # data model
    #
    logit(rho[site,1]) <- inprod(rho_beta, X[site,])
    y[site,1] ~ dbin(rho[site,1] * z[site,1], J)
    #
    # For remaining years of sampling
    #
    for(year in 2:nyear){
      #
      # latent state model, has theta term now.
      #
      logit(psi[site,year]) <- inprod(psi_beta, X[site,]) +
        theta * z[site,year-1]
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
  # First-order autoregressive term
  theta ~ dlogis(0,1)
}
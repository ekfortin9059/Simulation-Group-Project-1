#===============================================================================
#================= Gibbs Sampling for Driver Initial Locations =================
#===============================================================================
#' This R file contains code to generate samples from the posterior distribution
#' for the mean vector and covariance matrix parameters corresponding to the 
#' Bayesian model that assumes the driver initial locations are distributed from 
#' a multivariate normal distribution (where the mean vector and covariance 
#' matrix are to be estimated from the data). 
#===============================================================================
#===============================================================================


# Install Required Packages
require(rjags)
require(MASS)

# Read the data for rider drop-off locations
data <- read.csv("driver_oloc.csv")

# Extract coordinate data
coords <- cbind(data$X.Coordinate, data$Y.Coordinate)

# Number of data points
n <- nrow(coords)

# Define the JAGS model
model.string <- "
model{
  # Covariance matrix prior
  mu[1:2] ~ dmnorm(mu0[], Omega0[,])
  
  # Precision matrix prior
  Omega[1:2,1:2] ~ dwish(R[,], nu)
  
  # Convert precision matrix to covariance matrix
  Sigma[1:2,1:2] <- inverse(Omega[,])
   
  # Extract covariance matrix components to track
  Sigma11 <- Sigma[1,1]
  Sigma12 <- Sigma[1,2]
  Sigma22 <- Sigma[2,2]
  
  # Calculate Covariance
  rho <- Sigma12 / sqrt(Sigma11 * Sigma22)
  
  # Likelihood
  for (i in 1:n){
    coords[i, 1:2] ~ dmnorm(mu[], Omega[,])
  }
}
"

# Priors for mean vector and precision matrix
mu0 <- c(10, 11)
Omega0 <- matrix(c(0.1,0,0,0.1), nrow=2, ncol=2)

# Number of chains
n.chains <- 3

# Define the data for the JAGS model
jags.data <- list(coords=coords, mu0=mu0, Omega0=Omega0, n=n, nu=3, R=diag(2))

# Define initial values for the 3 chains
inits <- list(list(mu=c(5, 5), Omega=matrix(c(0.01,0,0,0.01), nrow=2, ncol=2)), 
              list(mu=c(10, 10), Omega=matrix(c(0.05,0,0,0.05), nrow=2, ncol=2)), 
              list(mu=c(15, 15), Omega=matrix(c(0.1,0,0,0.1), nrow=2, ncol=2)))

# Define parameters to track
params <- c("mu", "Sigma11", "Sigma12", "Sigma22", "rho")

# Compile the JAGS model
jags.mod <- jags.model(textConnection(model.string), data=jags.data, 
                                    n.chains=n.chains, inits=inits)

# Burn-in
update(jags.mod, 1000, progress.bar="none")

# Run 10000 MCMC samples for each chain
samp <- coda.samples(model=jags.mod, variable.names=params, 
                     n.iter=10000)

# Traceplots and density estimates for mean and covariance matrix components
plot(samp)

# Summary statistics
summary(samp)


# Posterior mean for the multivariate normal distribution parameters
mu_pred <- c(summary(samp)$statistics["mu[1]", "Mean"], summary(samp)$statistics["mu[2]", "Mean"])
Sigma_pred <- matrix(c(summary(samp)$statistics["Sigma11", "Mean"], 
                  summary(samp)$statistics["Sigma12", "Mean"],
                  summary(samp)$statistics["Sigma12", "Mean"],
                  summary(samp)$statistics["Sigma22", "Mean"]), nrow=2, ncol=2)

# Create replicate dataset and plot
yrep <- mvrnorm(n=n, mu=mu_pred, Sigma=Sigma_pred)
plot(yrep, cex=0.2, col="blue", xlim=c(0,20), ylim=c(0,20))
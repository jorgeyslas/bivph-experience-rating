##############################################
############### Prior regimes ################
######### Furrer, Yslas & Soerensen ##########
################### 2025 #####################
##############################################

source("generate_oe.R")

library(copula)

set.seed(555)

Ng <- read.csv("data/group_sizes.csv")[, 1]

regimes <- NULL

# Fixed variance of group effects
var_theta1 <- 0.3
var_theta2 <- 0.3

# Target values for Kendall's tau
tau_neg <- -0.5
tau_ind <- 0.0
tau_pos <- 0.5

# Generate group effects under a regime (tau, prior distribution)
r_regime <- function(id, tau, param1, param2, param3, param4, q_func, data) {
  # Copula parameter theta
  # Clayton: tau = theta / (theta + 2)
  theta <- 2 * tau / (1 - tau)
  cc <- claytonCopula(theta, dim = 2)

  # Draw uniform rvs
  u <- rCopula(G, copula = cc)

  # Group effect using quantile function
  theta1 <- q_func(u[, 1], param1, param2) # it has to be divided by the mean
  theta2 <- q_func(u[, 2], param3, param4)

  tau_ <- cor(theta1, theta2, method = "kendall")

  # Print Kendall's tau
  print(data.frame(id = id, tau = tau_, var_theta1 = var(theta1), var_theta2 = var(theta2)))

  # Left join group effects to data
  rbind(data, cbind(rep(id, G), 1:G, theta1, theta2))
}

###########
## Gamma ##
###########

# Parameters of the Gamma-distribution
psi1 <- 1 / var_theta1
psi2 <- 1 / var_theta2

regimes <- r_regime("gamma_neg", tau_neg, psi1, psi1, psi2, psi2, qgamma, regimes)
regimes <- r_regime("gamma_ind", tau_ind, psi1, psi1, psi2, psi2, qgamma, regimes)
regimes <- r_regime("gamma_pos", tau_pos, psi1, psi1, psi2, psi2, qgamma, regimes)


#######################
## Mixture of gammas ##
#######################

dgammamix <- function(x, shape1, rate1, shape2, rate2, prob) {
  prob * dgamma(x, shape1, rate1) + (1 - prob) * dgamma(x, shape2, rate2)
}

pgammamix <- function(x, shape1, rate1, shape2, rate2, prob) {
  prob * pgamma(x, shape1, rate1) + (1 - prob) * pgamma(x, shape2, rate2)
}

qgammamix <- function(p, shape1, rate1, shape2, rate2, prob) {
  L2 <- function(q, p) {
    (p - pgammamix(q, shape1, rate1, shape2, rate2, prob))^2
  }
  sapply(p, function(p) optimize(L2, c(0, 50), p = p)$minimum)
}

rgammamix <- function(n, shape1, rate1, shape2, rate2, prob) {
  mix <- rbinom(n, 1, prob)
  mix * rgamma(n, shape1, rate1) + (1 - mix) * rgamma(n, shape2, rate2)
}

mgammamix <- function(k, shape1, rate1, shape2, rate2, prob) {
  integrate(function(x) x^k * dgammamix(x, shape1, rate1, shape2, rate2, prob), 0, Inf)$value
}

x <- seq(0, 5, by = 0.01)
plot(x, dgammamix(x, 5, 2, 2, 6, 0.85), col = "red", type = "l")
mean1 <- mgammamix(1, 5, 2, 2, 6, 0.85)


# Generate group effects under a regime (tau, prior distribution)
r_regime_mg <- function(id, tau, param1, param2, param3, param4, param5, param6, param7, q_func1, q_func2, mean1, data) {
  # Copula parameter theta
  # Clayton: tau = theta / (theta + 2)
  theta <- 2 * tau / (1 - tau)
  cc <- claytonCopula(theta, dim = 2)

  # Draw uniform rvs
  u <- rCopula(G, copula = cc)

  # Group effect using quantile function
  theta1 <- q_func1(u[, 1], param1, param2, param3, param4, param5) / mean1 # it has to be divided by the mean
  theta2 <- q_func2(u[, 2], param6, param7)

  tau_ <- cor(theta1, theta2, method = "kendall")

  # Print Kendall's tau
  print(data.frame(id = id, tau = tau_, var_theta1 = var(theta1), var_theta2 = var(theta2)))

  # Left join group effects to data
  rbind(data, cbind(rep(id, G), 1:G, theta1, theta2))
}

regimes <- r_regime_mg("mgamma_ind", tau_ind, 5, 2, 2, 6, 0.85, psi2, psi2, qgammamix, qgamma, mean1, regimes)
regimes <- r_regime_mg("mgamma_pos", tau_pos, 5, 2, 2, 6, 0.85, psi2, psi2, qgammamix, qgamma, mean1, regimes)
regimes <- r_regime_mg("mgamma_neg", tau_neg, 5, 2, 2, 6, 0.85, psi2, psi2, qgammamix, qgamma, mean1, regimes)

regimes <- data.frame(
  id = regimes[, 1],
  group = as.integer(regimes[, 2]),
  dis_theta = as.numeric(regimes[, 3]),
  rec_theta = as.numeric(regimes[, 4])
)


# Save as .csv
write.csv(regimes, "data/regimes.csv", row.names = FALSE)

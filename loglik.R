##############################################
########### Loglikelihood functions ##########
######### Furrer, Yslas & Soerensen ##########
################### 2025 #####################
##############################################

# Simple model
simple_ll <- function(oe_data, oe_data_fit) {
  psi1 <- oe_data_fit[1, 10]
  psi2 <- oe_data_fit[1, 11]
  oe_data$o1_hat <- exp(oe_data_fit[1, 12] + oe_data_fit[1, 14] * oe_data$age + oe_data_fit[1, 16] * oe_data$age^2) * (oe_data$e1 + 1e-15)
  oe_data$o2_hat <- exp(oe_data_fit[1, 13] + oe_data_fit[1, 15] * oe_data$age) * (oe_data$e2 + 1e-15)

  simple_ag <- aggregate(cbind(o1, o1_hat, o2, o2_hat) ~ group + no + id, data = oe_data, FUN = sum)

  ll1 <- sum(lgamma(simple_ag$o1 + psi1) + psi1 * log(psi1) - lgamma(psi1) - (simple_ag$o1 + psi1) * log(psi1 + simple_ag$o1_hat)) + sum(oe_data[, 6] * log(oe_data[, 9])) - sum(lgamma(oe_data[, 6] + 1))

  ll2 <- sum(lgamma(simple_ag$o2 + psi2) + psi2 * log(psi2) - lgamma(psi2) - (simple_ag$o2 + psi2) * log(psi2 + simple_ag$o2_hat)) + sum(oe_data[, 8] * log(oe_data[, 10])) - sum(lgamma(oe_data[, 8] + 1))

  return(ll1 + ll2)
}

# Hierarchical model
hier_ll <- function(oe_data, oe_data_fit) {
  nu <- oe_data_fit[1, 10]
  eta <- oe_data_fit[1, 11]
  cons <- eta / nu
  oe_data$o1_hat <- exp(oe_data_fit[1, 12] + oe_data_fit[1, 14] * oe_data$age + oe_data_fit[1, 16] * oe_data$age^2) * (oe_data$e1 + 1e-15)
  oe_data$o2_hat <- exp(oe_data_fit[1, 13] + oe_data_fit[1, 15] * oe_data$age) * (oe_data$e2 + 1e-15)

  hier_ag <- aggregate(cbind(o1, o1_hat, o2, o2_hat) ~ group + no + id, data = oe_data, FUN = sum)

  comp1 <- sum(oe_data[, 6] * log(oe_data[, 9])) - sum(lgamma(oe_data[, 6] + 1)) + sum(oe_data[, 8] * log(oe_data[, 10])) - sum(lgamma(oe_data[, 8] + 1))

  aux_int <- function(theta, y1, y2, f1, f2) {
    exp(lgamma(theta + y1) - 2 * lgamma(theta) + 2 * theta * log(cons) - (theta + y1) * log(cons + f1) + lgamma(theta + y2) - (theta + y2) * log(cons + f2))
  }

  comp2 <- 0
  for (i in 1:(length(hier_ag[, 1]))) {
    comp2 <- comp2 + log(integrate(function(theta) aux_int(theta, hier_ag$o1[i], hier_ag$o2[i], hier_ag$o1_hat[i], hier_ag$o2_hat[i]) * dgamma(theta, eta, nu), 0, Inf)$value)
  }

  return(comp1 + comp2)
}

# Phase-type model
ph_ll <- function(oe_data, oe_data_fit, alpha, S11, S12, S22) {
  oe_data$o1_hat <- exp(oe_data_fit[1, 12] + oe_data_fit[1, 14] * oe_data$age + oe_data_fit[1, 16] * oe_data$age^2) * (oe_data$e1 + 1e-15)
  oe_data$o2_hat <- exp(oe_data_fit[1, 13] + oe_data_fit[1, 15] * oe_data$age) * (oe_data$e2 + 1e-15)

  aggreated_data <- aggregate(cbind(o1, o1_hat, o2, o2_hat) ~ group + no + id, data = oe_data, FUN = sum)

  s1 <- sum(lgamma(aggreated_data[, 4] + 1)) + sum(lgamma(aggreated_data[, 6] + 1)) - sum(lgamma(oe_data[, 6] + 1)) - sum(lgamma(oe_data[, 8] + 1))
  s2 <- sum(oe_data[, 6] * log(oe_data[, 9])) + sum(oe_data[, 8] * log(oe_data[, 10]))
  s3 <- sum((aggreated_data[, 4]) * log(aggreated_data[, 5])) + sum((aggreated_data[, 6]) * log(aggreated_data[, 7]))
  s4 <- sum(log(mp_cor_dens_cov(as.matrix(aggreated_data[, c(4, 6)]), as.matrix(aggreated_data[, c(5, 7)]), alpha, S11, S12, S22) / (aggreated_data[, 5] * aggreated_data[, 7])))

  return(s1 + s2 - s3 + s4)
}

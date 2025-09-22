##############################################
############## Estimation under ##############
######### Hierarchical Gamma Priors ##########
######### Furrer, Yslas & Soerensen ##########
################### 2025 #####################
##############################################

source("get_mixing_coef.R")

RIGHT <- function(x, n) {
  substr(x, nchar(x) - n + 1, nchar(x))
}

sens <- 1e-4 # sensivity hyperhyperparameter

run_glm <- 10 # skip Poisson regression for every "run_glm"'th step

# ECM-algorithm
estimate_hierarchical <- function(oe_data) {
  G <- max(oe_data$group)

  fn_min <- function(nu, eta, tmean0, smean0, smean1, smean2) {
    -eta * log(nu) + lgamma(eta) - eta * tmean0 + nu * smean0 - 2 * log(eta / nu) * smean0 + (eta / nu) * (smean1 + smean2)
  }

  hierarchical <- NULL

  # Univariate Poisson for initial values of regression parameters
  m1 <- glm(o1 ~ age + I(age^2) + offset(log(e1 + 1e-10)), family = poisson(link = "log"), data = oe_data)
  m2 <- glm(o2 ~ age + offset(log(e2 + 1e-10)), family = poisson(link = "log"), data = oe_data)

  oe_data$o1.hat <- fitted(m1)
  oe_data$o2.hat <- fitted(m2)

  # Initial values
  nu0 <- 1
  eta0 <- 1
  eta.0 <- nu.0 <- 0
  k <- 1 * run_glm
  stop <- F

  hierarchical <- aggregate(cbind(o1, o1.hat, o2, o2.hat) ~ group + no + id, data = oe_data, FUN = sum)

  hierarchical$N <- hierarchical$o1 + hierarchical$o2

  lim <- 80 # use proxy if no. occurrences > lim -- to avoid numerical overflow

  hierarchical$N_norm <- hierarchical$N / lim

  hierarchical$o1_norm <- round(hierarchical$o1 / hierarchical$N_norm)
  hierarchical$o1_norm <- ifelse(hierarchical$N > lim, hierarchical$o1_norm, hierarchical$o1)

  hierarchical$o2_norm <- round(hierarchical$o2 / hierarchical$N_norm)
  hierarchical$o2_norm <- ifelse(hierarchical$N > lim, hierarchical$o2_norm, hierarchical$o2)

  hierarchical$N_alt <- hierarchical$o1_norm + hierarchical$o2_norm

  while (abs(eta.0 - eta0) > sens | abs(nu.0 - nu0) > sens) {
    tmp <- aggregate(cbind(o1.hat, o2.hat) ~ group + no + id, data = oe_data, FUN = sum)

    hierarchical$o1.hat <- tmp$o1.hat
    hierarchical$o2.hat <- tmp$o2.hat

    hierarchical$o1.hat_norm <- round(hierarchical$o1.hat / hierarchical$N_norm)
    hierarchical$o1.hat_norm <- ifelse(hierarchical$N > lim, hierarchical$o1.hat_norm, hierarchical$o1.hat)

    hierarchical$o2.hat_norm <- round(hierarchical$o2.hat / hierarchical$N_norm)
    hierarchical$o2.hat_norm <- ifelse(hierarchical$N > lim, hierarchical$o2.hat_norm, hierarchical$o2.hat)

    hierarchical$p <- log(1 + (nu0 / eta0) * hierarchical$o1.hat) + log(1 + (nu0 / eta0) * hierarchical$o2.hat)
    hierarchical$p_norm <- log(1 + (nu0 / eta0) * hierarchical$o1.hat_norm) + log(1 + (nu0 / eta0) * hierarchical$o2.hat_norm)
    hierarchical$p_norm <- ifelse(hierarchical$N > lim, hierarchical$p_norm, hierarchical$p)

    for (i in 1:G) {
      w <- NULL

      arg1 <- min(hierarchical$o1_norm[hierarchical$group == i], hierarchical$o2_norm[hierarchical$group == i])
      arg2 <- max(hierarchical$o1_norm[hierarchical$group == i], hierarchical$o2_norm[hierarchical$group == i])
      arg3 <- hierarchical$p_norm[hierarchical$group == i]
      arg4 <- hierarchical$N_alt[hierarchical$group == i]

      w <- get_mixing_coef(arg1, arg2, arg3, arg4, eta0, nu0)
      w. <- sum(w)

      t0 <- s0 <- 0
      for (j in 1:length(w)) {
        t0 <- t0 + (w[j] / w.) * (digamma(j - 1 + eta0) - log(arg3 + nu0))
        s0 <- s0 + (w[j] / w.) * ((j - 1 + eta0) / (arg3 + nu0))
      }

      if (is.na(s0) | is.na(t0)) {
        stop <- T
        break
      }

      oe_data$t0[oe_data$group == i] <- hierarchical$t0[hierarchical$group == i] <- t0
      oe_data$s0[oe_data$group == i] <- hierarchical$s0[hierarchical$group == i] <- s0
      oe_data$s1[oe_data$group == i] <- hierarchical$s1[hierarchical$group == i] <- (hierarchical$o1[hierarchical$group == i] + s0) / (hierarchical$o1.hat[hierarchical$group == i] + eta0 / nu0)
      oe_data$s2[oe_data$group == i] <- hierarchical$s2[hierarchical$group == i] <- (hierarchical$o2[hierarchical$group == i] + s0) / (hierarchical$o2.hat[hierarchical$group == i] + eta0 / nu0)
    }

    if (stop) {
      break
    }

    if (k %% run_glm == 0 | k < 2 * run_glm) {
      m1 <- glm(o1 ~ age + I(age^2) + offset(log(e1 * s1 + 1e-10)), family = poisson(link = "log"), data = oe_data)
      m2 <- glm(o2 ~ age + offset(log(e2 * s2 + 1e-10)), family = poisson(link = "log"), data = oe_data)

      beta_fit1 <- unname(coef(m1))
      oe_data$o1.hat <- exp(beta_fit1[1] + beta_fit1[2] * oe_data$age + beta_fit1[3] * oe_data$age^2) * oe_data$e1

      beta_fit2 <- unname(coef(m2))
      oe_data$o2.hat <- exp(beta_fit2[1] + beta_fit2[2] * oe_data$age) * oe_data$e2
    }

    # nu0 <- optim(nu.0, fn_min, eta = eta0, tmean0 = mean(hierarchical$t0), smean0 = mean(hierarchical$s0), smean1 = mean(hierarchical$s1), smean2 = mean(hierarchical$s2), method = "Brent", lower = 0, upper = 200)$par
    # eta0 <- optim(eta.0, fn_min, nu = nu0, tmean0 = mean(hierarchical$t0), smean0 = mean(hierarchical$s0), smean1 = mean(hierarchical$s1), smean2 = mean(hierarchical$s2), method = "Brent", lower = 0, upper = 200)$par

    fn_tmp <- function(p) {
      nu_ <- p[1]
      eta_ <- p[2]
      fn_min(nu_, eta_, mean(hierarchical$t0), mean(hierarchical$s0), mean(hierarchical$s1), mean(hierarchical$s2))
    }

    nu.0 <- nu0
    eta.0 <- eta0

    op <- NULL
    if (inherits(try(optim(c(nu.0, eta.0), fn_tmp, method = "L-BFGS-B", lower = 0, upper = 200), silent = TRUE), "try-error")) {
      break
    } else {
      op <- optim(c(nu.0, eta.0), fn_tmp, method = "L-BFGS-B", lower = 0, upper = 200)
    }

    nu0 <- op$par[1]
    eta0 <- op$par[2]

    # Alternative: Newton-Rhapson
    # nu0 <- (-mean(hierarchical$t0) + eta0 / 2 + sqrt((eta0 / 2 - mean(hierarchical$s0))^2 + mean(hierarchical$s0) * eta0 * ((mean(hierarchical$s1) + mean(hierarchical$s2))))) / (mean(hierarchical$s0))
    # eta0 <- eta0 - (1 / 2 * (log(nu0) - digamma(eta0) + mean(hierarchical$t0)) + mean(hierarchical$s0) / eta0 - ((mean(hierarchical$s1) + mean(hierarchical$s2)) / 2) / nu0) / (-mean(hierarchical$s0) / (eta0^2) - trigamma(eta0) / 2)

    # Log-likelihood
    old <- fn_min(nu.0, eta.0, mean(hierarchical$t0), mean(hierarchical$s0), mean(hierarchical$s1), mean(hierarchical$s2))
    new <- fn_min(nu0, eta0, mean(hierarchical$t0), mean(hierarchical$s0), mean(hierarchical$s1), mean(hierarchical$s2))

    k <- k + 1
  }

  oe_out <- hierarchical[, 1:7]

  oe_out$theta1 <- (hierarchical$o1 + hierarchical$s0) / (hierarchical$o1.hat + eta0 / nu0)
  oe_out$theta2 <- (hierarchical$o2 + hierarchical$s0) / (hierarchical$o2.hat + eta0 / nu0)

  oe_out$nu <- nu0
  oe_out$eta <- eta0

  m1 <- glm(o1 ~ age + I(age^2) + offset(log(e1 * s1 + 1e-10)), family = poisson(link = "log"), data = oe_data)
  m2 <- glm(o2 ~ age + offset(log(e2 * s2 + 1e-10)), family = poisson(link = "log"), data = oe_data)

  c1 <- coef(m1)
  c2 <- coef(m2)

  # Intercept
  oe_out$int1 <- c1[1]
  oe_out$int2 <- c2[1]

  # (log) Linear term
  oe_out$lin1 <- c1[2]
  oe_out$lin2 <- c2[2]

  # (log) Second order term
  oe_out$sec1 <- c1[3]

  print(c(oe_out$id[1], oe_out$no[1]))

  oe_out
}

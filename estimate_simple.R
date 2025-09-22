##############################################
############## Estimation under ##############
############ Simple gamma priors #############
######### Furrer, Yslas & Soerensen ##########
################### 2025 #####################
##############################################

sens1 <- 1e-5 # sensivity psi1
sens2 <- 1e-5 # sensivity psi2
# h <- 0.05 # Newton-Rhapson step size for psi1 and psi2

run_glm <- 4 # skip Poisson regression for every "run_glm"'th step

estimate_simple <- function(oe_data) {
  G <- max(oe_data$group)

  fn_min <- function(psi, tmean, smean) {
    lgamma(psi) - psi * log(psi) - psi * (tmean - smean)
  }

  simple_gamma <- NULL

  # Initial values
  m1 <- glm(o1 ~ age + I(age^2) + offset(log(e1 + 1e-10)), family = poisson(link = "log"), data = oe_data)
  m2 <- glm(o2 ~ age + offset(log(e2 + 1e-10)), family = poisson(link = "log"), data = oe_data)

  oe_data$o1_hat <- fitted(m1)
  oe_data$o2_hat <- fitted(m2)

  psi1 <- 5
  psi2 <- 5
  psi.1 <- psi.2 <- 0

  i <- 1 * run_glm

  while (abs(psi.1 - psi1) > sens1 | abs(psi.2 - psi2) > sens2) {
    # simple_gamma <- aggregate(cbind(o1, o1_hat, e2, o2, o2_hat, lgamma(o1 + 1), lgamma(o2 + 1), o1 * log(o1_hat), o2 * log(o2_hat)) ~ group + no + id, data = oe_data, FUN = sum)
    simple_gamma <- aggregate(cbind(o1, o1_hat, o2, o2_hat) ~ group + no + id, data = oe_data, FUN = sum)

    # log-likelihood
    # print(sum(lgamma(simple_gamma$o1 + psi1) + psi1 * log(psi1) + simple_gamma$V8 - lgamma(psi1) - simple_gamma$V6 - (simple_gamma$o1 + psi1) * log(psi1 + simple_gamma$o1_hat)))

    if (abs(psi.1 - psi1) > sens1) {
      # Posterior means : E[Theta] & E[log(Theta)]
      simple_gamma$s1 <- (simple_gamma$o1 + psi1) / (simple_gamma$o1_hat + psi1)
      simple_gamma$t1 <- digamma(simple_gamma$o1 + psi1) - log(simple_gamma$o1_hat + psi1)

      oe_data$s1 <- rep(simple_gamma$s1[1:G], sapply(1:G, function(g) {
        length(oe_data$age[oe_data$group == g])
      }))
      oe_data$t1 <- rep(simple_gamma$t1[1:G], sapply(1:G, function(g) {
        length(oe_data$age[oe_data$group == g])
      }))

      if (i %% run_glm == 0) {
        m1 <- glm(o1 ~ age + I(age^2) + offset(log(e1 * s1 + 1e-10)), family = poisson(link = "log"), data = oe_data)
        beta_fit1 <- unname(coef(m1))
        oe_data$o1_hat <- exp(beta_fit1[1] + beta_fit1[2] * oe_data$age + beta_fit1[3] * oe_data$age^2) * oe_data$e1
      }

      psi.1 <- psi1
      psi1 <- optim(psi.1, fn_min, tmean = mean(simple_gamma$t1), smean = mean(simple_gamma$s1), method = "Brent", lower = 0, upper = 20)$par

      # Alternative - Newton-Rhapson
      # psi1 <- psi1 - h * (digamma(psi1) - log(psi1) - mean(simple_gamma$t1) + mean(simple_gamma$s1) - 1) / (trigamma(psi1) - 1 / psi1)

      # print(psi1)
    }

    if (abs(psi.2 - psi2) > sens2) {
      # log-likelihood
      # print(sum(lgamma(simple_gamma$o2 + psi2) + psi2 * log(psi2) + simple_gamma$V9 - lgamma(psi2) - simple_gamma$V7 - (simple_gamma$o2 + psi2) * log(psi2 + simple_gamma$o2_hat)))

      # Posterior means : E[Theta] & E[log(Theta)]
      simple_gamma$s2 <- (simple_gamma$o2 + psi2) / (simple_gamma$o2_hat + psi2)
      simple_gamma$t2 <- digamma(simple_gamma$o2 + psi2) - log(simple_gamma$o2_hat + psi2)

      oe_data$s2 <- rep(simple_gamma$s2[1:G], sapply(1:G, function(g) {
        length(oe_data$age[oe_data$group == g])
      }))
      oe_data$t2 <- rep(simple_gamma$t2[1:G], sapply(1:G, function(g) {
        length(oe_data$age[oe_data$group == g])
      }))

      if (i %% run_glm == 0) {
        m2 <- glm(o2 ~ age + offset(log(e2 * s2 + 1e-10)), family = poisson(link = "log"), data = oe_data)
        beta_fit2 <- unname(coef(m2))
        oe_data$o2_hat <- exp(beta_fit2[1] + beta_fit2[2] * oe_data$age) * oe_data$e2
      }

      psi.2 <- psi2
      psi2 <- optim(psi.2, fn_min, tmean = mean(simple_gamma$t2), smean = mean(simple_gamma$s2), method = "Brent", lower = 0, upper = 20)$par

      # Alternative - Newton-Rhapson
      # psi2 <- psi2 - h * (digamma(psi2) - log(psi2) - mean(oe_data$t2) + mean(oe_data$s2) - 1) / (trigamma(psi2) - 1 / psi2)
      # print(psi2)
    }

    i <- i + 1

    # print(c(psi1,psi2))
  }

  oe_out <- simple_gamma[, 1:7]

  oe_out$theta1 <- (simple_gamma$o1 + psi1) / (simple_gamma$o1_hat + psi1)
  oe_out$theta2 <- (simple_gamma$o2 + psi2) / (simple_gamma$o2_hat + psi2)

  oe_out$psi1 <- psi1
  oe_out$psi2 <- psi2

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

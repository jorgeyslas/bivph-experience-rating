##############################################
############## Estimation under ##############
################# no mixing ##################
######### Furrer, Yslas & Soerensen ##########
################### 2025 #####################
##############################################

estimate_std <- function(oe_data) {
  m1 <- glm(o1 ~ age + I(age^2) + offset(log(e1 + 1e-6)), family = poisson(link = "log"), data = oe_data)
  m2 <- glm(o2 ~ age + offset(log(e2 + 1e-6)), family = poisson(link = "log"), data = oe_data)

  oe_data$o1_hat <- fitted(m1)
  oe_data$o2_hat <- fitted(m2)

  oe_out <- aggregate(cbind(o1, o1_hat, o2, o2_hat) ~ group + no + id, data = oe_data, FUN = sum)

  oe_out$theta1 <- 1
  oe_out$theta2 <- 1

  oe_out$psi1 <- 0
  oe_out$psi2 <- 0

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

  oe_out
}

estimate_fac <- function(oe_data) {
  m1 <- glm(o1 ~ as.factor(group) + age + I(age^2) + offset(log(e1 + 1e-6)), family = poisson(link = "log"), data = oe_data)
  m2 <- glm(o2 ~ as.factor(group) + age + offset(log(e2 + 1e-6)), family = poisson(link = "log"), data = oe_data)

  oe_data$o1_hat <- fitted(m1)
  oe_data$o2_hat <- fitted(m2)

  oe_out <- aggregate(cbind(o1, o1_hat, o2, o2_hat) ~ group + no + id, data = oe_data, FUN = sum)

  oe_out$theta1 <- exp(c(coef(m1)[1], coef(m1)[1] + coef(m1)[2:100]))
  oe_out$theta2 <- exp(c(coef(m2)[1], coef(m2)[1] + coef(m2)[2:100]))

  oe_out$psi1 <- Inf
  oe_out$psi2 <- Inf

  # Intercept
  oe_out$int1 <- 0
  oe_out$int2 <- 0

  # (log) Linear term
  oe_out$lin1 <- coef(m1)[101]
  oe_out$lin2 <- coef(m2)[101]

  # (log) Second order term
  oe_out$sec1 <- coef(m1)[102]

  oe_out
}

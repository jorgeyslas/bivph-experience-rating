##############################################
################# RK4 solver #################
######### Furrer, Yslas & Soerensen ##########
################### 2025 #####################
##############################################

source("thinning.R")

tau <- 3

# Thiele

R_ <- function(a, s, V, u01, u10, theta1, theta2) {
  return(c(
    {
      (r + mu02(s)) * V[1] - theta1 * u01(s) * (V[2] - V[1]) * ifelse(s - a <= tau, 1, 0)
    },
    {
      (r + mu12(s)) * V[2] - 1 - theta2 * u10(s) * (V[1] - V[2])
    }
  ))
}

RK4 <- function(h, a, b, uu01, uu10, ttheta1, ttheta2) {
  # Initialize
  G1 <- G2 <- G3 <- G4 <- NULL

  a_ <- round(a, 0)

  if (a_ == b) {
    return(0)
  }

  n <- round((a_ - b) / h, 0)
  n_pre <- round(ifelse(b - a_ <= tau, n, (a_ + tau - b) / h), 0) # steps to reach discontinuity at tau

  f <- matrix(NA, n + 1, 2)
  f[1, ] <- c(0, 0) # start

  for (i in 1:n_pre)
  {
    s <- b + h * (i - 1)
    # 4'th order Runge Kutta scheme
    G1 <- h * R_(a_, s, f[i, ], uu01, uu10, ttheta1, ttheta2)
    G2 <- h * R_(a_, s + h / 2, f[i, ] + G1 / 2, uu01, uu10, ttheta1, ttheta2)
    G3 <- h * R_(a_, s + h / 2, f[i, ] + G2 / 2, uu01, uu10, ttheta1, ttheta2)
    G4 <- h * R_(a_, s + h, f[i, ] + G3, uu01, uu10, ttheta1, ttheta2)

    f[i + 1, ] <- f[i, ] + (G1 + 2 * G2 + 2 * G3 + G4) / 6
  }

  if (n_pre != n) {
    f[n_pre + 2, ] <- f[n_pre + 1, ]

    for (i in (n_pre + 2):n)
    {
      s <- b + h * (i - 1)
      # 4'th order Runge Kutta scheme
      G1 <- h * R_(a_, s, f[i, ], uu01, uu10, ttheta1, ttheta2)
      G2 <- h * R_(a_, s + h / 2, f[i, ] + G1 / 2, uu01, uu10, ttheta1, ttheta2)
      G3 <- h * R_(a_, s + h / 2, f[i, ] + G2 / 2, uu01, uu10, ttheta1, ttheta2)
      G4 <- h * R_(a_, s + h, f[i, ] + G3, uu01, uu10, ttheta1, ttheta2)

      f[i + 1, ] <- f[i, ] + (G1 + 2 * G2 + 2 * G3 + G4) / 6
    }
  }

  f[ifelse(n_pre != n, n + 1, n_pre + 1), 1]
}

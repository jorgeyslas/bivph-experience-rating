##############################################
################# Validation #################
######### Furrer, Yslas & Soerensen ##########
################### 2025 #####################
##############################################

source("RK4_solver.R")
source("generate_oe.R")

portfoliodata <- read.csv("data/portfoliodata.csv")
regimes <- read.csv("data/regimes.csv")

true <- regimes
true$theta1 <- true$dis_theta
true$theta2 <- true$rec_theta
true$int1 <- -4.5
true$int2 <- 0.3
true$lin1 <- -0.018
true$lin2 <- -0.049
true$sec1 <- 0.00064
true <- subset(true, id == "gamma_ind" | id == "mgamma_neg" | id == "mgamma_pos")
true <- split(true, paste(true$id, true$no))

est_std <- read.csv("data/est_std.csv")
est_std <- split(est_std, paste(est_std$id, est_std$no))

est_fixed <- read.csv("data/est_fixed.csv")
est_fixed$theta2 <- ifelse(est_fixed$theta2 > 5, 5, est_fixed$theta2)
est_fixed <- split(est_fixed, paste(est_fixed$id, est_fixed$no))

est_simple <- read.csv("data/est_simple.csv")
est_simple <- split(est_simple, paste(est_simple$id, est_simple$no))

est_hierarchical <- read.csv("data/est_hierarchical.csv")
est_hierarchical <- split(est_hierarchical, paste(est_hierarchical$id, est_hierarchical$no))

est_ph_gamma_ind <- read.csv("data/est_ph_gamma_ind.csv")
est_ph_mgamma_neg <- read.csv("data/est_ph_mgamma_neg.csv")
est_ph_mgamma_pos <- read.csv("data/est_ph_mgamma_pos.csv")
est_ph <- rbind(est_ph_gamma_ind, est_ph_mgamma_pos, est_ph_mgamma_neg)
est_ph <- split(est_ph, paste(est_ph$id, est_ph$no))

r <- 0.01 # interest rate
h <- -1 / 12 # step length in years

G <- max(portfoliodata$group)
H <- length(portfoliodata$age_at_entry)

portfoliodata$age <- as.integer(portfoliodata$age_at_entry)
portfoliodata$Id <- portfoliodata$group * 1000 + portfoliodata$age
portfoliodata$interpol_w <- portfoliodata$age_at_entry - portfoliodata$age

# Reduce no of data points
pf1 <- pf2 <- portfoliodata[, -2]
pf2$age <- pf1$age + 1
pf <- rbind(pf1, pf2)
pf$Id <- pf$group * 1000 + pf$age
pf <- pf[!duplicated(pf$Id), ]
rownames(pf) <- NULL

# Calculate premiums per insured
calc_premium <- function(est_data) {
  # est_data <- true$`gamma_ind ` # test
  V0 <- NULL

  mu01_ <- function(x) {
    exp(est_data$int1[1] + est_data$lin1[1] * x + est_data$sec1[1] * x^2)
  }

  mu10_ <- function(x) {
    exp(est_data$int2[1] + est_data$lin2[1] * x)
  }

  RK4_ <- function(xh, a, b, ttheta1, ttheta2) {
    RK4(xh, a, b, mu01_, mu10_, ttheta1, ttheta2)
  }

  pf$theta1 <- est_data$theta1[match(pf$group, est_data$group)]
  pf$theta2 <- est_data$theta2[match(pf$group, est_data$group)]

  pf$V0 <- mapply(RK4_, rep(h, length(pf$group)), pf$age, rep(x_max, length(pf$group)), pf$theta1, pf$theta2)

  # Linear interpolation
  portfoliodata$V0_left <- pf$V0[match(portfoliodata$Id, pf$Id)]
  portfoliodata$V0_right <- pf$V0[match(portfoliodata$Id + 1, pf$Id)]
  portfoliodata$V0 <- portfoliodata$V0_left * (1 - portfoliodata$interpol_w) + portfoliodata$V0_right * portfoliodata$interpol_w

  print(est_data$no[1])

  portfoliodata$V0
}

prem_true <- lapply(true, calc_premium)
prem_true <- do.call(rbind, prem_true)
prem_true <- c(t(prem_true))
write.csv(prem_true, "data/premiums_true_incl_linear_interpolation.csv", row.names = FALSE)

prem_std <- lapply(est_std, calc_premium)
prem_std <- do.call(rbind, prem_std)
prem_std <- c(t(prem_std))
write.csv(prem_std, "data/prem_std.csv", row.names = FALSE)

prem_fixed <- lapply(est_fixed, calc_premium)
prem_fixed <- do.call(rbind, prem_fixed)
prem_fixed <- c(t(prem_fixed))
write.csv(prem_fixed, "data/prem_fixed.csv", row.names = FALSE)

prem_simple <- lapply(est_simple, calc_premium)
prem_simple <- do.call(rbind, prem_simple)
prem_simple <- c(t(prem_simple))
write.csv(prem_simple, "data/prem_simple.csv", row.names = FALSE)

prem_hierarchical <- lapply(est_hierarchical, calc_premium)
prem_hierarchical <- do.call(rbind, prem_hierarchical)
prem_hierarchical <- c(t(prem_hierarchical))
write.csv(prem_hierarchical, "data/prem_hierarchical.csv", row.names = FALSE)

prem_ph <- lapply(est_ph, calc_premium)
prem_ph <- do.call(rbind, prem_ph)
prem_ph <- c(t(prem_ph))
write.csv(prem_ph, "data/prem_ph.csv", row.names = FALSE)

prem_true <- scan("data/premiums_true_incl_linear_interpolation.csv", skip = 1, dec = ".", sep = ",")

prem_true_gamma_ind <- prem_true[1:50000]
prem_true_mgamma_neg <- prem_true[50001:100000]
prem_true_mgamma_pos <- prem_true[100001:150000]
write.csv(prem_true_gamma_ind, "data/prem_true_gamma_ind.csv", row.names = FALSE)
write.csv(prem_true_mgamma_neg, "data/prem_true_mgamma_neg.csv", row.names = FALSE)
write.csv(prem_true_mgamma_pos, "data/prem_true_mgamma_pos.csv", row.names = FALSE)

prem_std <- scan("data/prem_std_20240527.csv", skip = 1, dec = ".", sep = ",")

prem_std_gamma_ind <- prem_std[1:5000000]
prem_std_mgamma_neg <- prem_std[5000001:10000000]
prem_std_mgamma_pos <- prem_std[10000001:15000000]
write.csv(prem_std_gamma_ind, "data/prem_std_gamma_ind.csv", row.names = FALSE)
write.csv(prem_std_mgamma_neg, "data/prem_std_mgamma_neg.csv", row.names = FALSE)
write.csv(prem_std_mgamma_pos, "data/prem_std_mgamma_pos.csv", row.names = FALSE)

prem_fixed <- scan("data/prem_fixed_20240527.csv", skip = 1, dec = ".", sep = ",")

prem_fixed_gamma_ind <- prem_fixed[1:5000000]
prem_fixed_mgamma_neg <- prem_fixed[5000001:10000000]
prem_fixed_mgamma_pos <- prem_fixed[10000001:15000000]
write.csv(prem_fixed_gamma_ind, "data/prem_fixed_gamma_ind.csv", row.names = FALSE)
write.csv(prem_fixed_mgamma_neg, "data/prem_fixed_mgamma_neg.csv", row.names = FALSE)
write.csv(prem_fixed_mgamma_pos, "data/prem_fixed_mgamma_pos.csv", row.names = FALSE)

prem_simple <- scan("data/prem_simple_20240527.csv", skip = 1, dec = ".", sep = ",")

prem_simple_gamma_ind <- prem_simple[1:5000000]
prem_simple_mgamma_neg <- prem_simple[5000001:10000000]
prem_simple_mgamma_pos <- prem_simple[10000001:15000000]
write.csv(prem_simple_gamma_ind, "data/prem_simple_gamma_ind.csv", row.names = FALSE)
write.csv(prem_simple_mgamma_neg, "data/prem_simple_mgamma_neg.csv", row.names = FALSE)
write.csv(prem_simple_mgamma_pos, "data/prem_simple_mgamma_pos.csv", row.names = FALSE)

prem_hierarchical <- scan("data/prem_hierarchical_20240527.csv", skip = 1, dec = ".", sep = ",")

prem_hierarchical_gamma_ind <- prem_hierarchical[1:5000000]
prem_hierarchical_mgamma_neg <- prem_hierarchical[5000001:10000000]
prem_hierarchical_mgamma_pos <- prem_hierarchical[10000001:15000000]
write.csv(prem_hierarchical_gamma_ind, "data/prem_hierarchical_gamma_ind.csv", row.names = FALSE)
write.csv(prem_hierarchical_mgamma_neg, "data/prem_hierarchical_mgamma_neg.csv", row.names = FALSE)
write.csv(prem_hierarchical_mgamma_pos, "data/prem_hierarchical_mgamma_pos.csv", row.names = FALSE)

prem_ph <- scan("data/prem_ph.csv", skip = 1, dec = ".", sep = ",")

prem_ph_gamma_ind <- prem_ph[1:5000000]
prem_ph_mgamma_neg <- prem_ph[5000001:10000000]
prem_ph_mgamma_pos <- prem_ph[10000001:15000000]
write.csv(prem_ph_gamma_ind, "data/prem_ph_gamma_ind.csv", row.names = FALSE)
write.csv(prem_ph_mgamma_neg, "data/prem_ph_mgamma_neg.csv", row.names = FALSE)
write.csv(prem_ph_mgamma_pos, "data/prem_ph_mgamma_pos.csv", row.names = FALSE)

## MAE AND RMSE
options(digits = 16)

## True
prem_true_gamma_ind <- scan("data/prem_true_gamma_ind.csv", skip = 1, dec = ".", sep = ",")
prem_true_mgamma_neg <- scan("data/prem_true_mgamma_neg.csv", skip = 1, dec = ".", sep = ",")
prem_true_mgamma_pos <- scan("data/prem_true_mgamma_pos.csv", skip = 1, dec = ".", sep = ",")

## Standard
prem_std_gamma_ind <- scan("data/prem_std_gamma_ind.csv", skip = 1, dec = ".", sep = ",")
prem_std_mgamma_neg <- scan("data/prem_std_mgamma_neg.csv", skip = 1, dec = ".", sep = ",")
prem_std_mgamma_pos <- scan("data/prem_std_mgamma_pos.csv", skip = 1, dec = ".", sep = ",")

# Gamma ind
sqrt(mean((prem_std_gamma_ind - rep(prem_true_gamma_ind, 100))^2)) # RMSE
mean(abs(prem_std_gamma_ind - rep(prem_true_gamma_ind, 100))) # MAE

# Mixed gamma neg
sqrt(mean((prem_std_mgamma_neg - rep(prem_true_mgamma_neg, 100))^2)) # RMSE
mean(abs(prem_std_mgamma_neg - rep(prem_true_mgamma_neg, 100))) # MAE

# Mixed gamma pos
sqrt(mean((prem_std_mgamma_pos - rep(prem_true_mgamma_pos, 100))^2)) # RMSE
mean(abs(prem_std_mgamma_pos - rep(prem_true_mgamma_pos, 100))) # MAE

## Fixed
prem_fixed_gamma_ind <- scan("data/prem_fixed_gamma_ind.csv", skip = 1, dec = ".", sep = ",")
prem_fixed_mgamma_neg <- scan("data/prem_fixed_mgamma_neg.csv", skip = 1, dec = ".", sep = ",")
prem_fixed_mgamma_pos <- scan("data/prem_fixed_mgamma_pos.csv", skip = 1, dec = ".", sep = ",")

# Gamma ind
sqrt(mean((prem_fixed_gamma_ind - rep(prem_true_gamma_ind, 100))^2)) # RMSE
mean(abs(prem_fixed_gamma_ind - rep(prem_true_gamma_ind, 100))) # MAE

# Mixed gamma neg
sqrt(mean((prem_fixed_mgamma_neg - rep(prem_true_mgamma_neg, 100))^2)) # RMSE
mean(abs(prem_fixed_mgamma_neg - rep(prem_true_mgamma_neg, 100))) # MAE

# Mixed gamma pos
sqrt(mean((prem_fixed_mgamma_pos - rep(prem_true_mgamma_pos, 100))^2)) # RMSE
mean(abs(prem_fixed_mgamma_pos - rep(prem_true_mgamma_pos, 100))) # MAE

## Simple
prem_simple_gamma_ind <- scan("data/prem_simple_gamma_ind.csv", skip = 1, dec = ".", sep = ",")
prem_simple_mgamma_neg <- scan("data/prem_simple_mgamma_neg.csv", skip = 1, dec = ".", sep = ",")
prem_simple_mgamma_pos <- scan("data/prem_simple_mgamma_pos.csv", skip = 1, dec = ".", sep = ",")

# Gamma ind
sqrt(mean((prem_simple_gamma_ind - rep(prem_true_gamma_ind, 100))^2)) # RMSE
mean(abs(prem_simple_gamma_ind - rep(prem_true_gamma_ind, 100))) # MAE

# Mixed gamma neg
sqrt(mean((prem_simple_mgamma_neg - rep(prem_true_mgamma_neg, 100))^2)) # RMSE
mean(abs(prem_simple_mgamma_neg - rep(prem_true_mgamma_neg, 100))) # MAE

# Mixed gamma pos
sqrt(mean((prem_simple_mgamma_pos - rep(prem_true_mgamma_pos, 100))^2)) # RMSE
mean(abs(prem_simple_mgamma_pos - rep(prem_true_mgamma_pos, 100))) # MAE

## Hierarchical
prem_hierarchical_gamma_ind <- scan("data/prem_hierarchical_gamma_ind.csv", skip = 1, dec = ".", sep = ",")
prem_hierarchical_mgamma_neg <- scan("data/prem_hierarchical_mgamma_neg.csv", skip = 1, dec = ".", sep = ",")
prem_hierarchical_mgamma_pos <- scan("data/prem_hierarchical_mgamma_pos.csv", skip = 1, dec = ".", sep = ",")

# Gamma ind
sqrt(mean((prem_hierarchical_gamma_ind - rep(prem_true_gamma_ind, 100))^2)) # RMSE
mean(abs(prem_hierarchical_gamma_ind - rep(prem_true_gamma_ind, 100))) # MAE

# Mixed gamma neg
sqrt(mean((prem_hierarchical_mgamma_neg - rep(prem_true_mgamma_neg, 100))^2)) # RMSE
mean(abs(prem_hierarchical_mgamma_neg - rep(prem_true_mgamma_neg, 100))) # MAE

# Mixed gamma pos
sqrt(mean((prem_hierarchical_mgamma_pos - rep(prem_true_mgamma_pos, 100))^2)) # RMSE
mean(abs(prem_hierarchical_mgamma_pos - rep(prem_true_mgamma_pos, 100))) # MAE

## Phase-type
prem_ph_gamma_ind <- scan("data/prem_ph_gamma_ind.csv", skip = 1, dec = ".", sep = ",")
prem_ph_mgamma_neg <- scan("data/prem_ph_mgamma_neg.csv", skip = 1, dec = ".", sep = ",")
prem_ph_mgamma_pos <- scan("data/prem_ph_mgamma_pos.csv", skip = 1, dec = ".", sep = ",")

# Gamma ind
sqrt(mean((prem_ph_gamma_ind - rep(prem_true_gamma_ind, 100))^2)) # RMSE
mean(abs(prem_ph_gamma_ind - rep(prem_true_gamma_ind, 100))) # MAE

# Mixed gamma neg
sqrt(mean((prem_ph_mgamma_neg - rep(prem_true_mgamma_neg, 100))^2)) # RMSE
mean(abs(prem_ph_mgamma_neg - rep(prem_true_mgamma_neg, 100))) # MAE

# Mixed gamma pos
sqrt(mean((prem_ph_mgamma_pos - rep(prem_true_mgamma_pos, 100))^2)) # RMSE
mean(abs(prem_ph_mgamma_pos - rep(prem_true_mgamma_pos, 100))) # MAE

#############################################
######### Estimation of phase-type ##########
######### Furrer, Yslas & Soerensen #########
################## 2025 #####################
#############################################

# library(devtools)
# install_github("https://github.com/jorgeyslas/phfrailty.git")

library(phfrailty)
source("ph_fitting.R")

### Fitting for scenario 22

# Case A: Independent Gamma
oe_data <- read.csv("data/oe_gamma.csv")

oe <- oe_data[oe_data$no == 22 & oe_data$id == "gamma_ind", ]

# Initial PH
set.seed(1)
bivph_ini <- bivphasetype(dimensions = c(3, 3))

mpph_fitA <- experience_rating_bph(bivph_ini, oe, stepsEM = 20, stepsPH = 20)


# Case B: Positively correlated scaled mixture
oe_data <- read.csv("data/oe_mgamma.csv")

oe <- oe_data[oe_data$no == 22 & oe_data$id == "mgamma_pos", ]

# Initial PH
set.seed(1)
bivph_ini <- bivphasetype(dimensions = c(6, 6))

mpph_fitB <- experience_rating_bph(bivph_ini, oe, stepsEM = 30, stepsPH = 20)


# Case C: Negatively correlated scaled mixture
oe <- oe_data[oe_data$no == 22 & oe_data$id == "mgamma_neg", ]

mpph_fitC <- experience_rating_bph(bivph_ini, oe, stepsEM = 30, stepsPH = 20, delta1 = 0.025, delta2 = 0.025)


### Computation for different scenarios
oe_data <- read.csv("data/oe_gamma.csv")

scenarios <- 1:100
case <- "gamma_ind"

# Initial PH
set.seed(1)
bivph_ini <- bivphasetype(dimensions = c(3, 3))

ph_par <- list()
est_ph <- data.frame()

for (i in scenarios) {
  oe <- oe_data[oe_data$no == i & oe_data$id == case, ]
  ph_fit_aux <- experience_rating_bph(bivph_ini, oe, stepsEM = 20, stepsPH = 20)
  est_ph <- rbind(est_ph, ph_fit_aux$oe_out)
  ph_par[[paste(case, i)]] <- ph_fit_aux$bph_fit
}

# write.csv(est_ph, "data/est_ph_gamma_ind.csv", row.names = FALSE)
# saveRDS(ph_par, file = "data/phfit_gamma_ind.rds")

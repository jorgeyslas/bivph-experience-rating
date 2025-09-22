##############################################
######### Loglikelihood calculation ##########
######### Furrer, Yslas & Soerensen ##########
################### 2025 #####################
##############################################

# library(devtools)
# install_github("https://github.com/jorgeyslas/phfrailty.git")

library(phfrailty)
source("loglik.R")

oe <- read.csv("data/oe_gamma.csv")
simple_fit <- read.csv("data/est_simple.csv")
hier_fit <- read.csv("data/est_hierarchical.csv")

scenarios <- 1:100

# Case A: Independent gamma
ph_fit <- read.csv("data/est_ph_gamma_ind.csv")
ph_fit_par <- readRDS("data/phfit_gamma_ind.rds")

case <- "gamma_ind"
ll_simple <- rep(0, 100)
ll_hier <- rep(0, 100)
ll_ph <- rep(0, 100)
for (i in scenarios) {
  oe_data <- oe[oe$id == case & oe$no == i, ]
  simple_fit_fil <- simple_fit[simple_fit$no == i & simple_fit$id %in% c(case), ]
  hier_fit_fil <- hier_fit[hier_fit$no == i & hier_fit$id %in% c(case), ]
  ph_fit_fil <- ph_fit[ph_fit$no == i & ph_fit$id %in% c(case), ]
  ph_fit_par_fil <- ph_fit_par[[paste(case, i)]]

  ll_simple[i] <- simple_ll(oe_data, simple_fit_fil)
  ll_hier[i] <- hier_ll(oe_data, hier_fit_fil)
  ll_ph[i] <- ph_ll(oe_data, ph_fit_fil, ph_fit_par_fil@pars$alpha, ph_fit_par_fil@pars$S11, ph_fit_par_fil@pars$S12, ph_fit_par_fil@pars$S22)
}
pdf("plots/ll_gamma_ind_all.pdf", width = 12, height = 6)
par(mar = c(5, 4.1, 4, 2) + 0.1)
plot(scenarios, ll_simple, xlab = "Scenario", ylab = "Loglikelihood", main = paste0("Loglikelihood for scenario (Case A)"), ylim = c(-4200, -3800), cex.main = 1.8, cex.lab = 1.5, cex.axis = 1.5)
points(scenarios, ll_hier, col = "red", pch = 2)
points(scenarios, ll_ph, col = "blue", pch = 3)
legend("bottomleft", inset = 0.05, c("Simple", "Hierarchical", "Phase-type"), col = c(1, "red", "blue"), pch = c(1, 2, 3), cex = 1.2)
dev.off()


oe <- read.csv("data/oe_mgamma.csv")

# Case B: Positive correlated scaled mixture
ph_fit <- read.csv("data/est_ph_mgamma_pos.csv")
ph_fit_par <- readRDS("data/phfit_mgamma_pos.rds")

case <- "mgamma_pos"
ll_simple <- rep(0, 100)
ll_hier <- rep(0, 100)
ll_ph <- rep(0, 100)
for (i in scenarios) {
  oe_data <- oe[oe$id == case & oe$no == i, ]
  simple_fit_fil <- simple_fit[simple_fit$no == i & simple_fit$id %in% c(case), ]
  hier_fit_fil <- hier_fit[hier_fit$no == i & hier_fit$id %in% c(case), ]
  ph_fit_fil <- ph_fit[ph_fit$no == i & ph_fit$id %in% c(case), ]
  ph_fit_par_fil <- ph_fit_par[[paste(case, i)]]

  ll_simple[i] <- simple_ll(oe_data, simple_fit_fil)
  ll_hier[i] <- hier_ll(oe_data, hier_fit_fil)
  ll_ph[i] <- ph_ll(oe_data, ph_fit_fil, ph_fit_par_fil@pars$alpha, ph_fit_par_fil@pars$S11, ph_fit_par_fil@pars$S12, ph_fit_par_fil@pars$S22)
}

pdf("plots/ll_mgamma_pos.pdf", width = 12, height = 6)
par(mar = c(5, 4.1, 4, 2) + 0.1)
plot(scenarios, ll_simple, xlab = "Scenario", ylab = "Loglikelihood", main = "Loglikelihood for scenario (Case B)", cex.main = 1.8, cex.lab = 1.5, cex.axis = 1.5)
points(scenarios, ll_hier, col = "red", pch = 2)
points(scenarios, ll_ph, col = "blue", pch = 3)
legend("bottomleft", inset = 0.05, c("Simple", "Hierarchical", "Phase-type"), col = c(1, "red", "blue"), pch = c(1, 2, 3), cex = 1.2)
dev.off()


# Case C: Negative correlated scaled mixture
ph_fit <- read.csv("data/est_ph_mgamma_neg.csv")
ph_fit_par <- readRDS("data/phfit_mgamma_neg.rds")

case <- "mgamma_neg"
ll_simple <- rep(0, 100)
ll_hier <- rep(0, 100)
ll_ph <- rep(0, 100)
for (i in scenarios) {
  oe_data <- oe[oe$id == case & oe$no == i, ]
  simple_fit_fil <- simple_fit[simple_fit$no == i & simple_fit$id %in% c(case), ]
  hier_fit_fil <- hier_fit[hier_fit$no == i & hier_fit$id %in% c(case), ]
  ph_fit_fil <- ph_fit[ph_fit$no == i & ph_fit$id %in% c(case), ]
  ph_fit_par_fil <- ph_fit_par[[paste(case, i)]]

  ll_simple[i] <- simple_ll(oe_data, simple_fit_fil)
  ll_hier[i] <- hier_ll(oe_data, hier_fit_fil)
  ll_ph[i] <- ph_ll(oe_data, ph_fit_fil, ph_fit_par_fil@pars$alpha, ph_fit_par_fil@pars$S11, ph_fit_par_fil@pars$S12, ph_fit_par_fil@pars$S22)
}
pdf("plots/ll_mgamma_neg.pdf", width = 12, height = 6)
par(mar = c(5, 4.1, 4, 2) + 0.1)
plot(scenarios, ll_simple, xlab = "Scenario", ylab = "Loglikelihood", main = paste0("Loglikelihood for scenario (Case C)"), ylim = c(-4050, -3600), cex.main = 1.8, cex.lab = 1.5, cex.axis = 1.5)
points(scenarios, ll_ph, col = "blue", pch = 3)
legend("bottomleft", inset = 0.05, c("Simple", "Phase-type"), col = c(1, "blue"), pch = c(1, 3), cex = 1.2)
dev.off()

pdf("plots/ll_mgamma_neg_all.pdf", width = 12, height = 6)
par(mar = c(5, 4.1, 4, 2) + 0.1)
plot(scenarios, ll_simple, xlab = "Scenario", ylab = "Loglikelihood", main = paste0("Loglikelihood for scenario (Case C)"), ylim = c(-4050, -3600), cex.main = 1.8, cex.lab = 1.5, cex.axis = 1.5)
points(scenarios, ll_hier, col = "red", pch = 2)
points(scenarios, ll_ph, col = "blue", pch = 3)
legend("bottomleft", inset = 0.05, c("Simple", "Hierarchical", "Phase-type"), col = c(1, "red", "blue"), pch = c(1, 2, 3), cex = 1.2)
dev.off()

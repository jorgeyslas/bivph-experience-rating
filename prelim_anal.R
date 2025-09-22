##############################################
########### Preliminary analysis #############
######### Furrer, Yslas & Soerensen ##########
################### 2025 #####################
##############################################

# True parameters
regimes <- read.csv("data/regimes.csv")

# Output for cases gamma_ind/mgamma_pos and scenario 22
oe_gamma <- read.csv("data/oe_gamma.csv")
oe_mgamma <- read.csv("data/oe_mgamma.csv")

oe <- rbind(oe_gamma, oe_mgamma)
oe <- split(oe, paste(oe$id, oe$no))

cases <- list(oe$`gamma_ind 22`, oe$`mgamma_pos 22`, oe$`mgamma_neg 22`)

source("estimate_nomix.R")

std_fit <- do.call(rbind, lapply(cases, estimate_std))
fac_fit <- do.call(rbind, lapply(cases, estimate_fac))

simple_fit <- read.csv("data/est_simple.csv")
simple_fit <- simple_fit[simple_fit$no == "22" & simple_fit$id %in% c("gamma_ind", "mgamma_pos", "mgamma_neg"), ]

hier_fit <- read.csv("data/est_hierarchical.csv")
hier_fit <- hier_fit[hier_fit$no == "22" & hier_fit$id %in% c("gamma_ind", "mgamma_pos", "mgamma_neg"), ]

ph_fit <- rbind(rbind(read.csv("data/est_ph_gamma_ind.csv"), read.csv("data/est_ph_mgamma_pos.csv")), read.csv("data/est_ph_mgamma_neg.csv"))
ph_fit <- ph_fit[ph_fit$no == "22" & ph_fit$id %in% c("gamma_ind", "mgamma_pos", "mgamma_neg"), ]

# Rates
true_dis <- function(x, g, case) {
  regimes[regimes$id == case & regimes$group == g, 3] * exp(-4.5 - 0.018 * x + 0.00064 * x^2)
}

true_dis_help <- function(x) {
  exp(-4.5 - 0.018 * x + 0.00064 * x^2)
}

std_dis <- function(x, case) {
  dat <- std_fit[std_fit$id == case & std_fit$group == 1, ]
  as.numeric(dat[8]) * exp(cbind(rep(1, length(x)), x, x^2) %*% t(dat[c(12, 14, 16)]))
}

fac_dis <- function(x, g, case) {
  dat <- fac_fit[fac_fit$id == case & fac_fit$group == g, ]
  as.numeric(dat[8]) * exp(cbind(rep(1, length(x)), x, x^2) %*% t(dat[c(12, 14, 16)]))
}

simple_dis <- function(x, g, case) {
  dat <- simple_fit[simple_fit$id == case & simple_fit$group == g, ]
  as.numeric(dat[8]) * exp(cbind(rep(1, length(x)), x, x^2) %*% t(dat[c(12, 14, 16)]))
}

hier_dis <- function(x, g, case) {
  dat <- hier_fit[hier_fit$id == case & hier_fit$group == g, ]
  as.numeric(dat[8]) * exp(cbind(rep(1, length(x)), x, x^2) %*% t(dat[c(12, 14, 16)]))
}

ph_dis <- function(x, g, case) {
  dat <- ph_fit[ph_fit$id == case & ph_fit$group == g, ]
  as.numeric(dat[8]) * exp(cbind(rep(1, length(x)), x, x^2) %*% t(dat[c(12, 14, 16)]))
}

true_rec <- function(x, g, case) {
  regimes[regimes$id == case & regimes$group == g, 4] * exp(0.3 - 0.049 * x)
}

true_rec_help <- function(x) {
  exp(0.3 - 0.049 * x)
}

std_rec <- function(x, case) {
  dat <- std_fit[std_fit$id == case & std_fit$group == 1, ]
  as.numeric(dat[9]) * exp(cbind(rep(1, length(x)), x) %*% t(dat[c(13, 15)]))
}

fac_rec <- function(x, g, case) {
  dat <- fac_fit[fac_fit$id == case & fac_fit$group == g, ]
  as.numeric(dat[9]) * exp(cbind(rep(1, length(x)), x) %*% t(dat[c(13, 15)]))
}

simple_rec <- function(x, g, case) {
  dat <- simple_fit[simple_fit$id == case & simple_fit$group == g, ]
  as.numeric(dat[9]) * exp(cbind(rep(1, length(x)), x) %*% t(dat[c(13, 15)]))
}

hier_rec <- function(x, g, case) {
  dat <- hier_fit[hier_fit$id == case & hier_fit$group == g, ]
  as.numeric(dat[9]) * exp(cbind(rep(1, length(x)), x) %*% t(dat[c(13, 15)]))
}

ph_rec <- function(x, g, case) {
  dat <- ph_fit[ph_fit$id == case & ph_fit$group == g, ]
  as.numeric(dat[9]) * exp(cbind(rep(1, length(x)), x) %*% t(dat[c(13, 15)]))
}

# Help
dis_weights_gamma_ind_22_fct <- function(g) {
  dat <- oe$`gamma_ind 22`
  dat <- dat[dat$group == g, 5]
  dat
}

rec_weights_gamma_ind_22_fct <- function(g) {
  dat <- oe$`gamma_ind 22`
  dat <- dat[dat$group == g, 7]
  dat
}

dis_weights_mgamma_pos_22_fct <- function(g) {
  dat <- oe$`mgamma_pos 22`
  dat <- dat[dat$group == g, 5]
  dat
}

rec_weights_mgamma_pos_22_fct <- function(g) {
  dat <- oe$`mgamma_pos 22`
  dat <- dat[dat$group == g, 7]
  dat
}

dis_weights_mgamma_neg_22_fct <- function(g) {
  dat <- oe$`mgamma_neg 22`
  dat <- dat[dat$group == g, 5]
  dat
}

rec_weights_mgamma_neg_22_fct <- function(g) {
  dat <- oe$`mgamma_neg 22`
  dat <- dat[dat$group == g, 7]
  dat
}

ages <- seq(20.5, 66.5, 1)

dis_norm_gamma_ind_22 <- NULL
for (g in 1:100) {
  z <- t(true_dis_help(ages)) %*% dis_weights_gamma_ind_22_fct(g)
  dis_norm_gamma_ind_22 <- c(dis_norm_gamma_ind_22, z)
}

rec_norm_gamma_ind_22 <- NULL
for (g in 1:100) {
  z <- t(true_rec_help(ages)) %*% rec_weights_gamma_ind_22_fct(g)
  rec_norm_gamma_ind_22 <- c(rec_norm_gamma_ind_22, z)
}

dis_norm_mgamma_pos_22 <- NULL
for (g in 1:100) {
  z <- t(true_dis_help(ages)) %*% dis_weights_mgamma_pos_22_fct(g)
  dis_norm_mgamma_pos_22 <- c(dis_norm_mgamma_pos_22, z)
}

rec_norm_mgamma_pos_22 <- NULL
for (g in 1:100) {
  z <- t(true_rec_help(ages)) %*% rec_weights_mgamma_pos_22_fct(g)
  rec_norm_mgamma_pos_22 <- c(rec_norm_mgamma_pos_22, z)
}

dis_norm_mgamma_neg_22 <- NULL
for (g in 1:100) {
  z <- t(true_dis_help(ages)) %*% dis_weights_mgamma_neg_22_fct(g)
  dis_norm_mgamma_neg_22 <- c(dis_norm_mgamma_neg_22, z)
}

rec_norm_mgamma_neg_22 <- NULL
for (g in 1:100) {
  z <- t(true_rec_help(ages)) %*% rec_weights_mgamma_neg_22_fct(g)
  rec_norm_mgamma_neg_22 <- c(rec_norm_mgamma_neg_22, z)
}

# Plots
ages <- seq(20, 67, 0.1)
groups <- c(23, 16, 68)

# Independent gammas (case A)
cas <- "gamma_ind"

pdf("plots/dis_thetas_scen22_gamma_ind.pdf", width = 6, height = 6)
par(mar = c(5, 4.1, 4, 2) + 0.1)
plot(regimes[regimes$id == cas, 3], simple_fit[simple_fit$id == cas, 8] * simple_fit[simple_fit$id == cas, 5] / dis_norm_gamma_ind_22, col = "blue", ylim = c(0, 4), xlim = c(0, 4), xlab = "True", ylab = "Estimate", main = "Comparison of disability Bayes \n estimates (Case A)", cex.main = 1.8, cex.lab = 1.5, cex.axis = 1.5)
points(regimes[regimes$id == cas, 3], fac_fit[fac_fit$id == cas, 5] / dis_norm_gamma_ind_22)
points(regimes[regimes$id == cas, 3], ph_fit[ph_fit$id == cas, 8] * ph_fit[ph_fit$id == cas, 5] / dis_norm_gamma_ind_22, col = "blue", pch = 6)
legend("bottomright", inset = 0.05, c("Simple", "Fixed", "Phase-type"), col = c("blue", 1, "blue"), pch = c(1, 1, 6), cex = 1.2)
# Add hierarchical
abline(0, 1)
dev.off()

# OBS: Lots of NaN's (for all groups with zero occurrences -- so the full power of shrinkage is not showcased)
pdf("plots/rec_thetas_scen22_gamma_ind.pdf", width = 6, height = 6)
par(mar = c(5, 4.1, 4, 2) + 0.1)
plot(regimes[regimes$id == cas, 4], simple_fit[simple_fit$id == cas, 9] * simple_fit[simple_fit$id == cas, 7] / rec_norm_gamma_ind_22, col = "blue", ylim = c(0, 5), xlim = c(0, 5), xlab = "True", ylab = "Estimate", main = "Comparison of recovery Bayes \n estimates (Case A)", cex.main = 1.8, cex.lab = 1.5, cex.axis = 1.5)
points(regimes[regimes$id == cas, 4], fac_fit[fac_fit$id == cas, 6] / rec_norm_gamma_ind_22)
points(regimes[regimes$id == cas, 4], ph_fit[ph_fit$id == cas, 9] * ph_fit[ph_fit$id == cas, 7] / rec_norm_gamma_ind_22, col = "blue", pch = 6)
legend("bottomright", inset = 0.05, c("Simple", "Fixed", "Phase-type"), col = c("blue", 1, "blue"), pch = c(1, 1, 6), cex = 1.2)
# Add hierarchical
abline(0, 1)
dev.off()

for (gr in groups) {
  pdf(paste0("plots/dis_rates_scen22_gamma_ind_group", gr, ".pdf"), width = 6, height = 6)
  plot(ages, true_dis(ages, gr, cas), type = "l", ylim = c(0, 0.12), xlab = "Age", ylab = "", main = paste0("Disability rates for group ", gr, "\n(Case A)"), xlim = c(20, 67), xaxt = "n", cex.main = 1.8, cex.lab = 1.5, cex.axis = 1.5)
  axis(side = 1, at = c(seq(20, 60, 10), 67), cex.axis = 1.5)
  lines(ages, std_dis(ages, cas), lty = 3)
  lines(ages, fac_dis(ages, gr, cas), lty = 3)
  lines(ages, simple_dis(ages, gr, cas), col = "blue")
  # lines(ages, hier_dis(ages, gr, cas), col = "blue", lty = 3)
  lines(ages, ph_dis(ages, gr, cas), col = "blue", lty = 5)
  legend("topleft", inset = 0.05, c("True", "Standard", "Fixed", "Simple", "Phase-type"), col = c(1, 1, 1, "blue", "blue"), lty = c(1, 3, 5, 1, 5), cex = 1.2)
  dev.off()

  pdf(paste0("plots/rec_rates_scen22_gamma_ind_group", gr, ".pdf"), width = 6, height = 6)
  plot(ages, true_rec(ages, gr, cas), type = "l", ylim = c(0, 0.8), xlab = "Age", ylab = "", main = paste0("Recovery rates for group ", gr, "\n(Case A)"), xlim = c(20, 67), xaxt = "n", cex.main = 1.8, cex.lab = 1.5, cex.axis = 1.5)
  axis(side = 1, at = c(seq(20, 60, 10), 67), cex.axis = 1.5)
  lines(ages, std_rec(ages, cas), lty = 3)
  lines(ages, fac_rec(ages, gr, cas), lty = 5)
  lines(ages, simple_rec(ages, gr, cas), col = "blue")
  # lines(ages, hier_rec(ages, gr, cas), col = "blue", lty = 3)
  lines(ages, ph_rec(ages, gr, cas), col = "blue", lty = 5)
  legend("topright", inset = 0.05, c("True", "Standard", "Fixed", "Simple", "Phase-type"), col = c(1, 1, 1, "blue", "blue"), lty = c(1, 3, 5, 1, 5), cex = 1.2)
  dev.off()
}


# Independent gammas (case A) - Including hierarchical
cas <- "gamma_ind"

pdf("plots/dis_thetas_scen22_gamma_ind_all.pdf", width = 6, height = 6)
par(mar = c(5, 4.1, 4, 2) + 0.1)
plot(regimes[regimes$id == cas, 3], simple_fit[simple_fit$id == cas, 8] * simple_fit[simple_fit$id == cas, 5] / dis_norm_gamma_ind_22, col = "blue", ylim = c(0, 4), xlim = c(0, 4), xlab = "True", ylab = "Estimate", main = "Comparison of disability Bayes \n estimates (Case A)", cex.main = 1.8, cex.lab = 1.5, cex.axis = 1.5)
points(regimes[regimes$id == cas, 3], fac_fit[fac_fit$id == cas, 5] / dis_norm_gamma_ind_22)
points(regimes[regimes$id == cas, 3], hier_fit[hier_fit$id == cas, 8] * hier_fit[hier_fit$id == cas, 5] / dis_norm_gamma_ind_22, col = "blue", pch = 2)
points(regimes[regimes$id == cas, 3], ph_fit[ph_fit$id == cas, 8] * ph_fit[ph_fit$id == cas, 5] / dis_norm_gamma_ind_22, col = "blue", pch = 6)
legend("bottomright", inset = 0.05, c("Simple", "Fixed", "Hierarchical", "Phase-type"), col = c("blue", 1, "blue", "blue"), pch = c(1, 1, 2, 6), cex = 1.2)
abline(0, 1)
dev.off()

pdf("plots/rec_thetas_scen22_gamma_ind_all.pdf", width = 6, height = 6)
par(mar = c(5, 4.1, 4, 2) + 0.1)
plot(regimes[regimes$id == cas, 4], simple_fit[simple_fit$id == cas, 9] * simple_fit[simple_fit$id == cas, 7] / rec_norm_gamma_ind_22, col = "blue", ylim = c(0, 5), xlim = c(0, 5), xlab = "True", ylab = "Estimate", main = "Comparison of recovery Bayes \n estimates (Case A)", cex.main = 1.8, cex.lab = 1.5, cex.axis = 1.5)
points(regimes[regimes$id == cas, 4], fac_fit[fac_fit$id == cas, 6] / rec_norm_gamma_ind_22)
points(regimes[regimes$id == cas, 4], hier_fit[hier_fit$id == cas, 9] * hier_fit[hier_fit$id == cas, 7] / rec_norm_gamma_ind_22, col = "blue", pch = 2)
points(regimes[regimes$id == cas, 4], ph_fit[ph_fit$id == cas, 9] * ph_fit[ph_fit$id == cas, 7] / rec_norm_gamma_ind_22, col = "blue", pch = 6)
legend("bottomright", inset = 0.05, c("Simple", "Fixed", "Hierarchical", "Phase-type"), col = c("blue", 1, "blue", "blue"), pch = c(1, 1, 2, 6), cex = 1.2)
abline(0, 1)
dev.off()

norm_sim_fit <- cbind(simple_fit[simple_fit$id == cas, 8], simple_fit[simple_fit$id == cas, 9])
cairo_pdf("plots/est_thetas_scen22_gamma_ind_sim.pdf", width = 6, height = 6)
par(mar = c(5, 4.1, 4, 2) + 0.1)
plot(norm_sim_fit, pch = 4,  main = "Estimated group effects (Case A)\n Simple model", xlab = "Disability", ylab = "Recovery", xlim = c(0, 4), ylim = c(0, 4), cex.main = 1.8, cex.lab = 1.5, cex.axis = 1.5)
points(norm_sim_fit[c(23, 16, 68), 1], norm_sim_fit[c(23, 16, 68), 2], col = "red", pch = 4)
text(norm_sim_fit[c(23, 16, 68), 1], norm_sim_fit[c(23, 16, 68), 2], labels = c("23", "16", "68"), pos = 3, offset = 0.35, cex = 1.2, col = adjustcolor("red", alpha.f = 0.7))
dev.off()

norm_hier_fit <- cbind(hier_fit[hier_fit$id == cas, 8], hier_fit[hier_fit$id == cas, 9])
cairo_pdf("plots/est_thetas_scen22_gamma_ind_hier.pdf", width = 6, height = 6)
par(mar = c(5, 4.1, 4, 2) + 0.1)
plot(norm_hier_fit, pch = 4,  main = "Estimated group effects (Case A)\n Hierarchical model", xlab = "Disability", ylab = "Recovery", xlim = c(0, 4), ylim = c(0, 4), cex.main = 1.8, cex.lab = 1.5, cex.axis = 1.5)
points(norm_hier_fit[c(23, 16, 68), 1], norm_hier_fit[c(23, 16, 68), 2], col = "red", pch = 4)
text(norm_hier_fit[c(23, 16, 68), 1], norm_hier_fit[c(23, 16, 68), 2], labels = c("23", "16", "68"), pos = 3, offset = 0.35, cex = 1.2, col = adjustcolor("red", alpha.f = 0.7))
dev.off()

norm_ph_fit <- cbind(ph_fit[ph_fit$id == cas, 8], ph_fit[ph_fit$id == cas, 9])
cairo_pdf("plots/est_thetas_scen22_gamma_ind_ph.pdf", width = 6, height = 6)
par(mar = c(5, 4.1, 4, 2) + 0.1)
plot(norm_ph_fit, pch = 4,  main = "Estimated group effects (Case A)\n Phase-type model", xlab = "Disability", ylab = "Recovery", xlim = c(0, 4), ylim = c(0, 4), cex.main = 1.8, cex.lab = 1.5, cex.axis = 1.5)
points(norm_ph_fit[c(23, 16, 68), 1], norm_ph_fit[c(23, 16, 68), 2], col = "red", pch = 4)
text(norm_ph_fit[c(23, 16, 68), 1], norm_ph_fit[c(23, 16, 68), 2], labels = c("23", "16", "68"), pos = 3, offset = 0.35, cex = 1.2, col = adjustcolor("red", alpha.f = 0.7))
dev.off()

for (gr in groups) {
  pdf(paste0("plots/dis_rates_scen22_gamma_ind_all_group", gr, ".pdf"), width = 6, height = 6)
  plot(ages, true_dis(ages, gr, cas), type = "l", ylim = c(0, 0.12), xlab = "Age", ylab = "", main = paste0("Disability rates for group ", gr, "\n(Case A)"), xlim = c(20, 67), xaxt = "n", cex.main = 1.8, cex.lab = 1.5, cex.axis = 1.5)
  axis(side = 1, at = c(seq(20, 60, 10), 67), cex.axis = 1.5)
  lines(ages, std_dis(ages, cas), lty = 3)
  lines(ages, fac_dis(ages, gr, cas), lty = 5)
  lines(ages, simple_dis(ages, gr, cas), col = "blue")
  lines(ages, hier_dis(ages, gr, cas), col = "blue", lty = 3)
  lines(ages, ph_dis(ages, gr, cas), col = "blue", lty = 5)
  legend("topleft", inset = 0.05, c("True", "Standard", "Fixed", "Simple", "Hierarchical", "Phase-type"), col = c(1, 1, 1, "blue", "blue", "blue"), lty = c(1, 3, 5, 1, 3, 5), cex = 1.2)
  dev.off()

  pdf(paste0("plots/rec_rates_scen22_gamma_ind_all_group", gr, ".pdf"), width = 6, height = 6)
  plot(ages, true_rec(ages, gr, cas), type = "l", ylim = c(0, 0.8), xlab = "Age", ylab = "", main = paste0("Recovery rates for group ", gr, "\n(Case A)"), xlim = c(20, 67), xaxt = "n", cex.main = 1.8, cex.lab = 1.5, cex.axis = 1.5)
  axis(side = 1, at = c(seq(20, 60, 10), 67), cex.axis = 1.5)
  lines(ages, std_rec(ages, cas), lty = 3)
  lines(ages, fac_rec(ages, gr, cas), lty = 5)
  lines(ages, simple_rec(ages, gr, cas), col = "blue")
  lines(ages, hier_rec(ages, gr, cas), col = "blue", lty = 3)
  lines(ages, ph_rec(ages, gr, cas), col = "blue", lty = 5)
  legend("topright", inset = 0.05, c("True", "Standard", "Fixed", "Simple", "Hierarchical", "Phase-type"), col = c(1, 1, 1, "blue", "blue", "blue"), lty = c(1, 3, 5, 1, 3, 5), cex = 1.2)
  dev.off()
}



# Positively correlated scaled mixture (case B)
cas <- "mgamma_pos"

pdf("plots/dis_thetas_scen22_mgamma_pos.pdf", width = 6, height = 6)
par(mar = c(5, 4.1, 4, 2) + 0.1)
plot(regimes[regimes$id == cas, 3], simple_fit[simple_fit$id == cas, 8] * simple_fit[simple_fit$id == cas, 5] / dis_norm_mgamma_pos_22, col = "blue", ylim = c(0, 4), xlim = c(0, 4), xlab = "True", ylab = "Estimate", main = "Comparison of disability Bayes \n estimates (Case B)", cex.main = 1.8, cex.lab = 1.5, cex.axis = 1.5)
points(regimes[regimes$id == cas, 3], fac_fit[fac_fit$id == cas, 5] / dis_norm_mgamma_pos_22)
points(regimes[regimes$id == cas, 3], hier_fit[hier_fit$id == cas, 8] * hier_fit[hier_fit$id == cas, 5] / dis_norm_mgamma_pos_22, col = "blue", pch = 2)
points(regimes[regimes$id == cas, 3], ph_fit[ph_fit$id == cas, 8] * ph_fit[ph_fit$id == cas, 5] / dis_norm_mgamma_pos_22, col = "blue", pch = 6)
legend("bottomright", inset = 0.05, c("Simple", "Fixed", "Hierarchical", "Phase-type"), col = c("blue", 1, "blue", "blue"), pch = c(1, 1, 2, 6), cex = 1.2)
abline(0, 1)
dev.off()

pdf("plots/rec_thetas_scen22_mgamma_pos.pdf", width = 6, height = 6)
par(mar = c(5, 4.1, 4, 2) + 0.1)
plot(regimes[regimes$id == cas, 4], simple_fit[simple_fit$id == cas, 9] * simple_fit[simple_fit$id == cas, 7] / rec_norm_mgamma_pos_22, col = "blue", ylim = c(0, 5), xlim = c(0, 5), xlab = "True", ylab = "Estimate", main = "Comparison of recovery Bayes \n estimates (Case B)", cex.main = 1.8, cex.lab = 1.5, cex.axis = 1.5)
points(regimes[regimes$id == cas, 4], fac_fit[fac_fit$id == cas, 6] / rec_norm_mgamma_pos_22)
points(regimes[regimes$id == cas, 4], hier_fit[hier_fit$id == cas, 9] * hier_fit[hier_fit$id == cas, 7] / rec_norm_mgamma_pos_22, col = "blue", pch = 2)
points(regimes[regimes$id == cas, 4], ph_fit[ph_fit$id == cas, 9] * ph_fit[ph_fit$id == cas, 7] / rec_norm_mgamma_pos_22, col = "blue", pch = 6)
legend("bottomright", inset = 0.05, c("Simple", "Fixed", "Hierarchical", "Phase-type"), col = c("blue", 1, "blue", "blue"), pch = c(1, 1, 2, 6), cex = 1.2)
abline(0, 1)
dev.off()

norm_sim_fit <- cbind(simple_fit[simple_fit$id == cas, 8], simple_fit[simple_fit$id == cas, 9])
cairo_pdf("plots/est_thetas_scen22_gamma_pos_sim.pdf", width = 6, height = 6)
par(mar = c(5, 4.1, 4, 2) + 0.1)
plot(norm_sim_fit, pch = 4,  main = "Estimated group effects (Case B)\n Simple model", xlab = "Disability", ylab = "Recovery", xlim = c(0, 3), ylim = c(0, 3), cex.main = 1.8, cex.lab = 1.5, cex.axis = 1.5)
points(norm_sim_fit[c(23, 16, 68), 1], norm_sim_fit[c(23, 16, 68), 2], col = "red", pch = 4)
text(norm_sim_fit[c(23, 16, 68), 1], norm_sim_fit[c(23, 16, 68), 2], labels = c("23", "16", "68"), pos = 3, offset = 0.35, cex = 1.2, col = adjustcolor("red", alpha.f = 0.7))
dev.off()

norm_hier_fit <- cbind(hier_fit[hier_fit$id == cas, 8], hier_fit[hier_fit$id == cas, 9])
cairo_pdf("plots/est_thetas_scen22_gamma_pos_hier.pdf", width = 6, height = 6)
par(mar = c(5, 4.1, 4, 2) + 0.1)
plot(norm_hier_fit, pch = 4,  main = "Estimated group effects (Case B)\n Hierarchical model", xlab = "Disability", ylab = "Recovery", xlim = c(0, 3), ylim = c(0, 3), cex.main = 1.8, cex.lab = 1.5, cex.axis = 1.5)
points(norm_hier_fit[c(23, 16, 68), 1], norm_hier_fit[c(23, 16, 68), 2], col = "red", pch = 4)
text(norm_hier_fit[c(23, 16, 68), 1], norm_hier_fit[c(23, 16, 68), 2], labels = c("23", "16", "68"), pos = 3, offset = 0.35, cex = 1.2, col = adjustcolor("red", alpha.f = 0.7))
dev.off()

norm_ph_fit <- cbind(ph_fit[ph_fit$id == cas, 8], ph_fit[ph_fit$id == cas, 9])
cairo_pdf("plots/est_thetas_scen22_gamma_pos_ph.pdf", width = 6, height = 6)
par(mar = c(5, 4.1, 4, 2) + 0.1)
plot(norm_ph_fit, pch = 4,  main = "Estimated group effects (Case B)\n Phase-type model", xlab = "Disability", ylab = "Recovery", xlim = c(0, 3), ylim = c(0, 3), cex.main = 1.8, cex.lab = 1.5, cex.axis = 1.5)
points(norm_ph_fit[c(23, 16, 68), 1], norm_ph_fit[c(23, 16, 68), 2], col = "red", pch = 4)
text(norm_ph_fit[c(23, 16, 68), 1], norm_ph_fit[c(23, 16, 68), 2], labels = c("23", "16", "68"), pos = 3, offset = 0.35, cex = 1.2, col = adjustcolor("red", alpha.f = 0.7))
dev.off()

for (gr in groups) {
  pdf(paste0("plots/dis_rates_scen22_mgamma_pos_group", gr, ".pdf"), width = 6, height = 6)
  plot(ages, true_dis(ages, gr, cas), type = "l", ylim = c(0, 0.12), xlab = "Age", ylab = "", main = paste0("Disability rates for group ", gr, "\n(Case B)"), xlim = c(20, 67), xaxt = "n", cex.main = 1.8, cex.lab = 1.5, cex.axis = 1.5)
  axis(side = 1, at = c(seq(20, 60, 10), 67), cex.axis = 1.5)
  lines(ages, std_dis(ages, cas), lty = 3)
  lines(ages, fac_dis(ages, gr, cas), lty = 5)
  lines(ages, simple_dis(ages, gr, cas), col = "blue")
  lines(ages, hier_dis(ages, gr, cas), col = "blue", lty = 3)
  lines(ages, ph_dis(ages, gr, cas), col = "blue", lty = 5)
  legend("topleft", inset = 0.05, c("True", "Standard", "Fixed", "Simple", "Hierarchical", "Phase-type"), col = c(1, 1, 1, "blue", "blue", "blue"), lty = c(1, 3, 5, 1, 3, 5), cex = 1.2)
  dev.off()

  pdf(paste0("plots/rec_rates_scen22_mgamma_pos_group", gr, ".pdf"), width = 6, height = 6)
  plot(ages, true_rec(ages, gr, cas), type = "l", ylim = c(0, 0.8), xlab = "Age", ylab = "", main = paste0("Recovery rates for group ", gr, "\n(Case B)"), xlim = c(20, 67), xaxt = "n", cex.main = 1.8, cex.lab = 1.5, cex.axis = 1.5)
  axis(side = 1, at = c(seq(20, 60, 10), 67), cex.axis = 1.5)
  lines(ages, std_rec(ages, cas), lty = 3)
  lines(ages, fac_rec(ages, gr, cas), lty = 5)
  lines(ages, simple_rec(ages, gr, cas), col = "blue")
  lines(ages, hier_rec(ages, gr, cas), col = "blue", lty = 3)
  lines(ages, ph_rec(ages, gr, cas), col = "blue", lty = 5)
  legend("topright", inset = 0.05, c("True", "Standard", "Fixed", "Simple", "Hierarchical", "Phase-type"), col = c(1, 1, 1, "blue", "blue", "blue"), lty = c(1, 3, 5, 1, 3, 5), cex = 1.2)
  dev.off()
}

# Negatively correlated scaled mixture (case C)
cas <- "mgamma_neg"

pdf("plots/dis_thetas_scen22_mgamma_neg.pdf", width = 6, height = 6)
par(mar = c(5, 4.1, 4, 2) + 0.1)
plot(regimes[regimes$id == cas, 3], simple_fit[simple_fit$id == cas, 8] * simple_fit[simple_fit$id == cas, 5] / dis_norm_mgamma_neg_22, col = "blue", ylim = c(0, 4), xlim = c(0, 4), xlab = "True", ylab = "Estimate", main = "Comparison of disability Bayes \n estimates (Case C)", cex.main = 1.8, cex.lab = 1.5, cex.axis = 1.5)
points(regimes[regimes$id == cas, 3], fac_fit[fac_fit$id == cas, 5] / dis_norm_mgamma_neg_22)
points(regimes[regimes$id == cas, 3], ph_fit[ph_fit$id == cas, 8] * ph_fit[ph_fit$id == cas, 5] / dis_norm_mgamma_neg_22, col = "blue", pch = 6)
legend("bottomright", inset = 0.05, c("Simple", "Fixed", "Phase-type"), col = c("blue", 1, "blue", "blue"), pch = c(1, 1, 6), cex = 1.2)
abline(0, 1)
dev.off()

pdf("plots/rec_thetas_scen22_mgamma_neg.pdf", width = 6, height = 6)
par(mar = c(5, 4.1, 4, 2) + 0.1)
plot(regimes[regimes$id == cas, 4], simple_fit[simple_fit$id == cas, 9] * simple_fit[simple_fit$id == cas, 7] / rec_norm_mgamma_neg_22, col = "blue", ylim = c(0, 5), xlim = c(0, 5), xlab = "True", ylab = "Estimate", main = "Comparison of recovery Bayes \n estimates (Case C)", cex.main = 1.8, cex.lab = 1.5, cex.axis = 1.5)
points(regimes[regimes$id == cas, 4], fac_fit[fac_fit$id == cas, 6] / rec_norm_mgamma_neg_22)
points(regimes[regimes$id == cas, 4], ph_fit[ph_fit$id == cas, 9] * ph_fit[ph_fit$id == cas, 7] / rec_norm_mgamma_neg_22, col = "blue", pch = 6)
legend("bottomright", inset = 0.05, c("Simple", "Fixed", "Phase-type"), col = c("blue", 1, "blue", "blue"), pch = c(1, 1, 6), cex = 1.2)
abline(0, 1)
dev.off()


for (gr in groups) {
  pdf(paste0("plots/dis_rates_scen22_mgamma_neg_group", gr, ".pdf"), width = 6, height = 6)
  plot(ages, true_dis(ages, gr, cas), type = "l", ylim = c(0, 0.12), xlab = "Age", ylab = "", main = paste0("Disability rates for group ", gr, "\n(Case C)"), xlim = c(20, 67), xaxt = "n", cex.main = 1.8, cex.lab = 1.5, cex.axis = 1.5)
  axis(side = 1, at = c(seq(20, 60, 10), 67), cex.axis = 1.5)
  lines(ages, std_dis(ages, cas), lty = 3)
  lines(ages, fac_dis(ages, gr, cas), lty = 5)
  lines(ages, simple_dis(ages, gr, cas), col = "blue")
  #lines(ages, hier_dis(ages, gr, cas), col = "blue", lty = 3)
  lines(ages, ph_dis(ages, gr, cas), col = "blue", lty = 5)
  legend("topleft", inset = 0.05, c("True", "Standard", "Fixed", "Simple", "Phase-type"), col = c(1, 1, 1, "blue", "blue", "blue"), lty = c(1, 3, 5, 1, 5), cex = 1.2)
  dev.off()
  
  pdf(paste0("plots/rec_rates_scen22_mgamma_neg_group", gr, ".pdf"), width = 6, height = 6)
  plot(ages, true_rec(ages, gr, cas), type = "l", ylim = c(0, 0.8), xlab = "Age", ylab = "", main = paste0("Recovery rates for group ", gr, "\n(Case C)"), xlim = c(20, 67), xaxt = "n", cex.main = 1.8, cex.lab = 1.5, cex.axis = 1.5)
  axis(side = 1, at = c(seq(20, 60, 10), 67), cex.axis = 1.5)
  lines(ages, std_rec(ages, cas), lty = 3)
  lines(ages, fac_rec(ages, gr, cas), lty = 5)
  lines(ages, simple_rec(ages, gr, cas), col = "blue")
  #lines(ages, hier_rec(ages, gr, cas), col = "blue", lty = 3)
  lines(ages, ph_rec(ages, gr, cas), col = "blue", lty = 5)
  legend("topright", inset = 0.05, c("True", "Standard", "Fixed", "Simple", "Phase-type"), col = c(1, 1, 1, "blue", "blue", "blue"), lty = c(1, 3, 5, 1, 5), cex = 1.2)
  dev.off()
}


# Negatively correlated scaled mixture (case C) - including Hierarchical
cas <- "mgamma_neg"

pdf("plots/dis_thetas_scen22_mgamma_neg_all.pdf", width = 6, height = 6)
par(mar = c(5, 4.1, 4, 2) + 0.1)
plot(regimes[regimes$id == cas, 3], simple_fit[simple_fit$id == cas, 8] * simple_fit[simple_fit$id == cas, 5] / dis_norm_mgamma_neg_22, col = "blue", ylim = c(0, 4), xlim = c(0, 4), xlab = "True", ylab = "Estimate", main = "Comparison of disability Bayes \n estimates (Case C)", cex.main = 1.8, cex.lab = 1.5, cex.axis = 1.5)
points(regimes[regimes$id == cas, 3], fac_fit[fac_fit$id == cas, 5] / dis_norm_mgamma_neg_22)
points(regimes[regimes$id == cas, 3], hier_fit[hier_fit$id == cas, 8] * hier_fit[hier_fit$id == cas, 5] / dis_norm_mgamma_neg_22, col = "blue", pch = 2)
points(regimes[regimes$id == cas, 3], ph_fit[ph_fit$id == cas, 8] * ph_fit[ph_fit$id == cas, 5] / dis_norm_mgamma_neg_22, col = "blue", pch = 6)
legend("bottomright", inset = 0.05, c("Simple", "Fixed", "Hierarchical", "Phase-type"), col = c("blue", 1, "blue", "blue"), pch = c(1, 1, 2, 6), cex = 1.2)
abline(0, 1)
dev.off()

pdf("plots/rec_thetas_scen22_mgamma_neg_all.pdf", width = 6, height = 6)
par(mar = c(5, 4.1, 4, 2) + 0.1)
plot(regimes[regimes$id == cas, 4], simple_fit[simple_fit$id == cas, 9] * simple_fit[simple_fit$id == cas, 7] / rec_norm_mgamma_neg_22, col = "blue", ylim = c(0, 5), xlim = c(0, 5), xlab = "True", ylab = "Estimate", main = "Comparison of recovery Bayes \n estimates (Case C)", cex.main = 1.8, cex.lab = 1.5, cex.axis = 1.5)
points(regimes[regimes$id == cas, 4], fac_fit[fac_fit$id == cas, 6] / rec_norm_mgamma_neg_22)
points(regimes[regimes$id == cas, 4], hier_fit[hier_fit$id == cas, 9] * hier_fit[hier_fit$id == cas, 7] / rec_norm_mgamma_neg_22, col = "blue", pch = 2)
points(regimes[regimes$id == cas, 4], ph_fit[ph_fit$id == cas, 9] * ph_fit[ph_fit$id == cas, 7] / rec_norm_mgamma_neg_22, col = "blue", pch = 6)
legend("bottomright", inset = 0.05, c("Simple", "Fixed", "Hierarchical", "Phase-type"), col = c("blue", 1, "blue", "blue"), pch = c(1, 1, 2, 6), cex = 1.2)
abline(0, 1)
dev.off()

norm_sim_fit <- cbind(simple_fit[simple_fit$id == cas, 8], simple_fit[simple_fit$id == cas, 9])
cairo_pdf("plots/est_thetas_scen22_gamma_neg_sim.pdf", width = 6, height = 6)
par(mar = c(5, 4.1, 4, 2) + 0.1)
plot(norm_sim_fit, pch = 4,  main = "Estimated group effects (Case C)\n Simple model", xlab = "Disability", ylab = "Recovery", xlim = c(0, 2.5), ylim = c(0, 3.5), cex.main = 1.8, cex.lab = 1.5, cex.axis = 1.5)
points(norm_sim_fit[c(23, 16, 68), 1], norm_sim_fit[c(23, 16, 68), 2], col = "red", pch = 4)
text(norm_sim_fit[c(23, 16, 68), 1], norm_sim_fit[c(23, 16, 68), 2], labels = c("23", "16", "68"), pos = 3, offset = 0.35, cex = 1.2, col = adjustcolor("red", alpha.f = 0.7))
dev.off()

norm_hier_fit <- cbind(hier_fit[hier_fit$id == cas, 8], hier_fit[hier_fit$id == cas, 9])
cairo_pdf("plots/est_thetas_scen22_gamma_neg_hier.pdf", width = 6, height = 6)
par(mar = c(5, 4.1, 4, 2) + 0.1)
plot(norm_hier_fit, pch = 4,  main = "Estimated group effects (Case C)\n Hierarchical model", xlab = "Disability", ylab = "Recovery", xlim = c(0, 2.5), ylim = c(0, 3.5), cex.main = 1.8, cex.lab = 1.5, cex.axis = 1.5)
points(norm_hier_fit[c(23, 16, 68), 1], norm_hier_fit[c(23, 16, 68), 2], col = "red", pch = 4)
text(norm_hier_fit[c(23, 16, 68), 1], norm_hier_fit[c(23, 16, 68), 2], labels = c("23", "16", "68"), pos = 3, offset = 0.35, cex = 1.2, col = adjustcolor("red", alpha.f = 0.7))
dev.off()

norm_ph_fit <- cbind(ph_fit[ph_fit$id == cas, 8], ph_fit[ph_fit$id == cas, 9])
cairo_pdf("plots/est_thetas_scen22_gamma_neg_ph.pdf", width = 6, height = 6)
par(mar = c(5, 4.1, 4, 2) + 0.1)
plot(norm_ph_fit, pch = 4,  main = "Estimated group effects (Case C)\n Phase-type model", xlab = "Disability", ylab = "Recovery", xlim = c(0, 2.5), ylim = c(0, 3.5), cex.main = 1.8, cex.lab = 1.5, cex.axis = 1.5)
points(norm_ph_fit[c(23, 16, 68), 1], norm_ph_fit[c(23, 16, 68), 2], col = "red", pch = 4)
text(norm_ph_fit[c(23, 16, 68), 1], norm_ph_fit[c(23, 16, 68), 2], labels = c("23", "16", "68"), pos = 3, offset = 0.35, cex = 1.2, col = adjustcolor("red", alpha.f = 0.7))
dev.off()


for (gr in groups) {
  pdf(paste0("plots/dis_rates_scen22_mgamma_neg_group_all", gr, ".pdf"), width = 6, height = 6)
  plot(ages, true_dis(ages, gr, cas), type = "l", ylim = c(0, 0.12), xlab = "Age", ylab = "", main = paste0("Disability rates for group ", gr, "\n(Case C)"), xlim = c(20, 67), xaxt = "n", cex.main = 1.8, cex.lab = 1.5, cex.axis = 1.5)
  axis(side = 1, at = c(seq(20, 60, 10), 67), cex.axis = 1.5)
  lines(ages, std_dis(ages, cas), lty = 3)
  lines(ages, fac_dis(ages, gr, cas), lty = 5)
  lines(ages, simple_dis(ages, gr, cas), col = "blue")
  lines(ages, hier_dis(ages, gr, cas), col = "blue", lty = 3)
  lines(ages, ph_dis(ages, gr, cas), col = "blue", lty = 5)
  legend("topleft", inset = 0.05, c("True", "Standard", "Fixed", "Simple", "Hierarchical", "Phase-type"), col = c(1, 1, 1, "blue", "blue", "blue"), lty = c(1, 3, 5, 1, 3, 5), cex = 1.2)
  dev.off()
  
  pdf(paste0("plots/rec_rates_scen22_mgamma_neg_group_all", gr, ".pdf"), width = 6, height = 6)
  plot(ages, true_rec(ages, gr, cas), type = "l", ylim = c(0, 0.8), xlab = "Age", ylab = "", main = paste0("Recovery rates for group ", gr, "\n(Case C)"), xlim = c(20, 67), xaxt = "n", cex.main = 1.8, cex.lab = 1.5, cex.axis = 1.5)
  axis(side = 1, at = c(seq(20, 60, 10), 67), cex.axis = 1.5)
  lines(ages, std_rec(ages, cas), lty = 3)
  lines(ages, fac_rec(ages, gr, cas), lty = 5)
  lines(ages, simple_rec(ages, gr, cas), col = "blue")
  lines(ages, hier_rec(ages, gr, cas), col = "blue", lty = 3)
  lines(ages, ph_rec(ages, gr, cas), col = "blue", lty = 5)
  legend("topright", inset = 0.05, c("True", "Standard", "Fixed", "Simple", "Hierarchical", "Phase-type"), col = c(1, 1, 1, "blue", "blue", "blue"), lty = c(1, 3, 5, 1, 3, 5), cex = 1.2)
  dev.off()
}

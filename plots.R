##############################################
################# Plotting ###################
######### Furrer, Yslas & Soerensen ##########
################### 2025 #####################
##############################################

# Load
portfoliodata <- read.csv("data/portfoliodata.csv")
regimes <- read.csv("data/regimes.csv")
Ng <- read.csv("data/group_sizes.csv")$Ng

# Source
source("generate_oe.R")
source("portfolio.R")

# Plots
# Log-scaled group sizes
pdf("plots/group_sizes.pdf", width = 12, height = 6)
par(mar = c(5, 4.1, 4, 2) + 0.1)
plot(Ng, log = "y", pch = 4, main = "Group sizes", xlab = "Group", ylab = "Size", cex.main = 1.8, cex.lab = 1.5, cex.axis = 1.5)
text(c(23, 16, 68), Ng[c(23, 16, 68)], labels = c("23", "16", "68"), pos = 3, offset = 0.35, cex = 1.2)
dev.off()

# Age distributions
pdf("plots/age_dist_true.pdf", width = 12, height = 6)
par(mar = c(5, 4.1, 4, 2) + 0.1)
plot(seq(20, 67, 0.1), f(seq(20, 67, 0.1)), type = "l", ylim = c(0, 0.08), xlab = "Age", ylab = "Density", main = "Distributions of initial ages", xlim = c(20, 67), xaxt = "n", cex.main = 1.8, cex.lab = 1.5, cex.axis = 1.5)
axis(side = 1, at = c(seq(20, 60, 10), 67), cex.axis = 1.5)
lines(seq(20, 67, 0.1), f1(seq(20, 67, 0.1)), type = "l", lty = 2)
lines(seq(20, 67, 0.1), f2(seq(20, 67, 0.1)), type = "l", lty = 3)
lines(seq(20, 67, 0.1), f3(seq(20, 67, 0.1)), type = "l", lty = 5)
legend("topright", inset = 0.05, c("Weibull", "Normal", "Gompertz", "Finite mixture"), lty = c(2, 3, 5, 1), col = c(1, 1, 1, 1), cex = 1.2)
dev.off()

# Age distribution of portfolio
pdf("plots/age_dist_portfolio.pdf", width = 6, height = 6)
par(mar = c(5, 4.1, 4, 2) + 0.1)
hist(portfoliodata$age_at_entry, prob = TRUE, breaks = 67 - 20 - 1, main = "Initial ages in the portfolio", xlim = c(20, 67), xlab = "Age", xaxt = "n", cex.main = 1.8, cex.lab = 1.5, cex.axis = 1.5)
axis(side = 1, at = c(seq(20, 60, 10), 67), cex.axis = 1.5)
lines(seq(20, 67, 0.1), f(seq(20, 67, 0.1)), lty = 2)
legend("topright", inset = 0.05, "Finite mixture (true)", lty = 2, col = 1, cex = 1.2)
dev.off()

# Variability in ages across groups
pdf("plots/age_dist_groups.pdf", width = 6, height = 6)
par(mar = c(5, 4.1, 4, 2) + 0.1)
mean_age <- aggregate(portfoliodata$age_at_entry, list(portfoliodata$group), FUN = mean)$x
plot(Ng, mean_age, log = "x", pch = 4, main = "Initial ages across groups", xlab = "Group size", ylab = "Average age", cex.main = 1.8, cex.lab = 1.5, cex.axis = 1.5)
text(Ng[c(23, 16, 68)], mean_age[c(23, 16, 68)], labels = c("23", "16", "68"), pos = 3, offset = 0.35, cex = 1.2)
dev.off()

# True rates
pdf("plots/activerates_true.pdf", width = 6, height = 6)
plot(seq(20, 67, 0.1), mu01(seq(20, 67, 0.1)), type = "l", xlab = "Age", main = "Transition rates from state active", ylab = "", ylim = c(0, 0.06), xlim = c(20, 67), xaxt = "n", cex.main = 1.8, cex.lab = 1.5, cex.axis = 1.5)
axis(side = 1, at = c(seq(20, 60, 10), 67), cex.axis = 1.5)
lines(seq(20, 67, 0.1), mu02(seq(20, 67, 0.1)), lty = 2)
legend("topleft", inset = 0.05, c("Disability rate", "Mortality rate"), lty = c(1, 2), col = c(1, 1), cex = 1.2)
dev.off()

pdf("plots/disabledrates_true.pdf", width = 6, height = 6)
plot(seq(20, 67, 0.1), mu10(seq(20, 67, 0.1)), type = "l", xlab = "Age", main = "Transition rates from state disabled", ylim = c(0, 0.5), ylab = "", xlim = c(20, 67), xaxt = "n", cex.main = 1.8, cex.lab = 1.5, cex.axis = 1.5)
axis(side = 1, at = c(seq(20, 60, 10), 67), cex.axis = 1.5)
lines(seq(20, 67, 0.1), mu12(seq(20, 67, 0.1)), lty = 2)
legend("topright", inset = 0.05, c("Recovery rate", "Mortality rate"), lty = c(1, 2), col = c(1, 1), cex = 1.2)
dev.off()

# Prior distributions
pdf("plots/priors.pdf", width = 12, height = 6)
plot(seq(0, 4, 0.01), dgamma(seq(0, 4, 0.01), 10 / 3, 10 / 3), type = "l", xlab = "", main = "Mixing densities", ylab = "", cex.main = 1.8, cex.lab = 1.5, cex.axis = 1.5)
lines(seq(0, 4, 0.01), 0.85 * 2.175 * dgamma(seq(0, 4, 0.01) * 2.175, 5, 2) + 2.175 * (1 - 0.85) * dgamma(seq(0, 4, 0.01) * 2.175, 2, 6), lty = 2)
legend("topright", inset = 0.05, c("Gamma", "Scaled mixture"), lty = c(1, 2), col = c(1, 1), cex = 1.2)
dev.off()

# Simulated random effects
sim_effect <- subset(regimes, id == "gamma_ind")
pdf("plots/hist_dis_ind.pdf", width = 6, height = 6)
par(mar = c(5, 4.1, 4, 2) + 0.1)
hist(sim_effect[, "dis_theta"], prob = TRUE, breaks = 20, main = "Case A: Disability group effect", xlim = c(0, 4), xlab = "", cex.main = 1.8, cex.lab = 1.5, cex.axis = 1.5)
dev.off()

pdf("plots/hist_rec_ind.pdf", width = 6, height = 6)
par(mar = c(5, 4.1, 4, 2) + 0.1)
hist(sim_effect[, "rec_theta"], prob = TRUE, breaks = 20, main = "Case A: Recovery group effect", xlim = c(0, 4.8), xlab = "", cex.main = 1.8, cex.lab = 1.5, cex.axis = 1.5)
dev.off()

cairo_pdf("plots/scatter_group_ind.pdf", width = 6, height = 6)
par(mar = c(5, 4.1, 4, 2) + 0.1)
plot(sim_effect[, 3], sim_effect[, 4], pch = 4, main = "Simulated group effects (Case A)", xlab = "Disability", ylab = "Recovery", xlim = c(0, 4), ylim = c(0, 4), cex.main = 1.8, cex.lab = 1.5, cex.axis = 1.5)
points(sim_effect[c(23, 16, 68), 3], sim_effect[c(23, 16, 68), 4], col = "red", pch = 4)
text(sim_effect[c(23, 16, 68), 3], sim_effect[c(23, 16, 68), 4], labels = c("23", "16", "68"), pos = 3, offset = 0.35, cex = 1.2, col = adjustcolor("red", alpha.f = 0.7))
dev.off()

round(mean(sim_effect[, "dis_theta"]), 6)
round(mean(sim_effect[, "rec_theta"]), 6)
round(sd(sim_effect[, "dis_theta"]), 6)
round(sd(sim_effect[, "rec_theta"]), 6)
round(cor(sim_effect[, "rec_theta"], sim_effect[, "dis_theta"], method = "kendall"), 6)

sim_effect <- subset(regimes, id == "mgamma_pos")
pdf("plots/hist_dis_pos.pdf", width = 6, height = 6)
par(mar = c(5, 4.1, 4, 2) + 0.1)
hist(sim_effect[, "dis_theta"], prob = TRUE, breaks = 20, main = "Case B: Disability group effect", xlim = c(0, 3), xlab = "", cex.main = 1.8, cex.lab = 1.5, cex.axis = 1.5)
dev.off()

pdf("plots/hist_rec_pos.pdf", width = 6, height = 6)
par(mar = c(5, 4.1, 4, 2) + 0.1)
hist(sim_effect[, "rec_theta"], prob = TRUE, breaks = 20, main = "Case B: Recovery group effect", xlim = c(0, 3), xlab = "", cex.main = 1.8, cex.lab = 1.5, cex.axis = 1.5)
dev.off()

cairo_pdf("plots/scatter_group_pos.pdf", width = 6, height = 6)
par(mar = c(5, 4.1, 4, 2) + 0.1)
plot(sim_effect[, 3], sim_effect[, 4], pch = 4, main = "Simulated group effects (Case B)", xlab = "Disability", ylab = "Recovery", xlim = c(0, 3), ylim = c(0, 3), cex.main = 1.8, cex.lab = 1.5, cex.axis = 1.5)
points(sim_effect[c(23, 16, 68), 3], sim_effect[c(23, 16, 68), 4], col = "red", pch = 4)
text(sim_effect[c(23, 16, 68), 3], sim_effect[c(23, 16, 68), 4], labels = c("23", "16", "68"), pos = 3, offset = 0.35, cex = 1.2, col = adjustcolor("red", alpha.f = 0.7))
dev.off()

round(mean(sim_effect[, "dis_theta"]), 6)
round(mean(sim_effect[, "rec_theta"]), 6)
round(sd(sim_effect[, "dis_theta"]), 6)
round(sd(sim_effect[, "rec_theta"]), 6)
round(cor(sim_effect[, "rec_theta"], sim_effect[, "dis_theta"], method = "kendall"), 6)

sim_effect <- subset(regimes, id == "mgamma_neg")
pdf("plots/hist_dis_neg.pdf", width = 6, height = 6)
par(mar = c(5, 4.1, 4, 2) + 0.1)
hist(sim_effect[, "dis_theta"], prob = TRUE, breaks = 20, main = "Case C: Disability group effect", xlim = c(0, 3), xlab = "", cex.main = 1.8, cex.lab = 1.5, cex.axis = 1.5)
dev.off()

pdf("plots/hist_rec_neg.pdf", width = 6, height = 6)
par(mar = c(5, 4.1, 4, 2) + 0.1)
hist(sim_effect[, "rec_theta"], prob = TRUE, breaks = 20, main = "Case C: Recovery group effect", xlim = c(0, 3), xlab = "", cex.main = 1.8, cex.lab = 1.5, cex.axis = 1.5)
dev.off()

cairo_pdf("plots/scatter_group_neg.pdf", width = 6, height = 6)
par(mar = c(5, 4.1, 4, 2) + 0.1)
plot(sim_effect[, 3], sim_effect[, 4], pch = 4, main = "Simulated group effects (Case C)", xlim = c(0, 2.5), ylim = c(0, 3.5), xlab = "Disability", ylab = "Recovery", cex.main = 1.8, cex.lab = 1.5, cex.axis = 1.5)
points(sim_effect[c(23, 16, 68), 3], sim_effect[c(23, 16, 68), 4], col = "red", pch = 4)
text(sim_effect[c(23, 16, 68), 3], sim_effect[c(23, 16, 68), 4], labels = c("23", "16", "68"), pos = 3, offset = 0.35, cex = 1.2, col = adjustcolor("red", alpha.f = 0.7))
dev.off()

round(mean(sim_effect[, "dis_theta"]), 6)
round(mean(sim_effect[, "rec_theta"]), 6)
round(sd(sim_effect[, "dis_theta"]), 6)
round(sd(sim_effect[, "rec_theta"]), 6)
round(cor(sim_effect[, "rec_theta"], sim_effect[, "dis_theta"], method = "kendall"), 6)

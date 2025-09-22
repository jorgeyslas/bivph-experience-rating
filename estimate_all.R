###################################################
######### Estimation of standard, fixed, ##########
######### simple and hierarchical models ##########
############ Furrer, Yslas & Soerensen ############
###################### 2025 #######################
###################################################

# Model comparison - scenarios

source("estimate_nomix.R")
source("estimate_simple.R")
source("estimate_hierarchical.R")

oe_gamma <- read.csv("data/oe_gamma.csv")
oe_mgamma <- read.csv("data/oe_mgamma.csv")

oe <- rbind(rbind(oe_gamma, oe_lnorm), oe_mgamma)
oe <- subset(oe, id == "gamma_ind" | id == "mgamma_neg" | id == "mgamma_pos")
oe <- split(oe, paste(oe$id, oe$no))

simple <- lapply(oe, estimate_simple)
simple <- do.call(rbind, simple)

hierarchical <- lapply(oe, estimate_hierarchical)
hierarchical <- do.call(rbind, hierarchical)

std <- lapply(oe, estimate_std)
std <- do.call(rbind, std)

fixed <- lapply(oe, estimate_fac)
fixed <- do.call(rbind, fixed)

# write.csv(std, "data/est_std.csv", row.names = FALSE)
# write.csv(fixed, "data/est_fixed.csv", row.names = FALSE)
# write.csv(simple, "data/est_simple.csv", row.names = FALSE)
# write.csv(hierarchical, "data/est_hierarchical.csv", row.names = FALSE)

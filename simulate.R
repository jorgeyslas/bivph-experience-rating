##############################################
################ Simulation ##################
######### Furrer, Yslas & Soerensen ##########
################### 2025 #####################
##############################################

# Load group effects/sizes
regimes <- read.csv("data/regimes.csv")
Ng <- read.csv("data/group_sizes.csv")[, 1]
G <- length(Ng)

# Source O/E-generator
source("generate_oe.R")

set.seed(9999)

# Number of simulations
N <- 100

group_effects <- split(regimes, regimes$id)
group_effects <- rep(group_effects, N)

r_oedata <- function(group_effects) {
  generate_oe(group_effects$id, group_effects$dis_theta, group_effects$rec_theta)
}

oe_data <- lapply(group_effects, r_oedata)
oe_data <- do.call(rbind, oe_data)

oe_data$no <- rep(1:N, each = length(oe_data$id) / N)

oe_data <- oe_data[, c("id", "no", "group", "age", "e1", "o1", "e2", "o2")]

agg <- aggregate(cbind(e1, o1, e2, o2) ~ age, data = oe_data, FUN = sum)

plot(agg$age, agg$o1 / agg$e1)
lines(agg$age, mu01(agg$age))
plot(agg$age, agg$o2 / agg$e2)
lines(agg$age, mu10(agg$age))

# Save as .csv
oe_gamma <- subset(oe_data, id == "gamma_ind" | id == "gamma_pos" | id == "gamma_neg")
oe_mgamma <- subset(oe_data, id == "mgamma_ind" | id == "mgamma_pos" | id == "mgamma_neg")

# write.csv(oe_gamma, "data/oe_gamma.csv", row.names = FALSE)
# write.csv(oe_mgamma, "data/oe_mgamma.csv", row.names = FALSE)

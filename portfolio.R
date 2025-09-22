#############################################
############ Portfolio generator ############
######### Furrer, Yslas & Soerensen #########
################## 2025 #####################
#############################################

# Preamble
rm(list = ls())
set.seed(6) # fixed seed
library(survival)
library(flexsurv) # needed for gompertz

#####################
## Age composition ##
#####################

xmin <- 20 # youngest age at entry
xmax <- 67 # oldest age at entry

# Young group - Weibull
wshape <- 3
wscale <- 15
wnorm <- pweibull(xmax - xmin, wshape, wscale) - pweibull(xmin - xmin, wshape, wscale)

# Middle-aged - Normal
nmean <- 42
nvar <- 7
nnorm <- pnorm(xmax, nmean, nvar) - pnorm(xmin, nmean, nvar)

# Old group - Gompertz
gshape <- 0.1
grate <- 1 / 500
gnorm <- pgompertz(xmax - xmin, gshape, grate) - pgompertz(xmin - xmin, gshape, grate)

# Weights
ww <- 0.25
wn <- 0.35
wg <- 0.4

# Normalization & new PDF
f1 <- function(x) {
  ifelse(xmin < x & x <= xmax, dweibull(x - xmin, wshape, wscale) / wnorm, 0)
}
f2 <- function(x) {
  ifelse(xmin < x & x <= xmax, dnorm(x, nmean, nvar) / nnorm, 0)
}
f3 <- function(x) {
  ifelse(xmin < x & x <= xmax, dgompertz(x - xmin, gshape, grate) / gnorm, 0)
}
f <- function(x) {
  ww * f1(x) + wn * f2(x) + wg * f3(x)
}

integrate(f, xmin, xmax)

#################
## Group specs ##
#################

G <- 100 # number of groups

# Group sizes, Ng ~ Gamma(0.05,0.0001)
# shape <- 0.05
# rate <- 0.0001
# Ng <- sample(5:5000, replace = TRUE, size = G, prob = dgamma(5:5000, shape, rate))
# write.csv(Ng, "data/group_sizes.csv", row.names = FALSE)

Ng <- read.csv("data/group_sizes.csv")[, 1]
H <- sum(Ng) # portfolio size

# Simulation of ages
ages <- NULL
for (i in 1:G) {
  age_profile <- sample(x = c("young", "middle", "old"), size = 1, prob = c(ww, wn, wg))
  h <- 0
  temp <- NULL
  if (age_profile == "young") {
    while (h < Ng[i]) {
      sim <- rweibull(1, wshape, wscale) + xmin
      if (sim > xmin & sim <= xmax) {
        h <- h + 1
        ages <- c(ages, sim)
      } else {
        h <- h
      }
    }
  }
  if (age_profile == "middle") {
    while (h < Ng[i]) {
      sim <- rnorm(1, nmean, nvar)
      if (sim > xmin & sim <= xmax) {
        h <- h + 1
        ages <- c(ages, sim)
      } else {
        h <- h
      }
    }
  }
  if (age_profile == "old") {
    while (h < Ng[i]) {
      sim <- rgompertz(1, gshape, grate) + xmin
      if (sim > xmin & sim <= xmax) {
        h <- h + 1
        ages <- c(ages, sim)
      } else {
        h <- h
      }
    }
  }
  ages
}

# Collecting portfolio data
portfoliodata <- data.frame(group = rep(1:G, Ng[1:G]), age_at_entry = ages)

# Save as csv
# write.csv(portfoliodata, "data/portfoliodata.csv", row.names = FALSE)

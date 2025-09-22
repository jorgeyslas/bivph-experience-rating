##############################################
############### O/E Generator ################
######### Furrer, Yslas & Soerensen ##########
################### 2025 #####################
##############################################

# Load portfolio data
portfoliodata <- read.csv("data/portfoliodata.csv")

x_min <- 20
x_max <- 67
cov_expiry <- 3

G <- max(portfoliodata$group)
H <- length(portfoliodata$age_at_entry)
Ng <- read.csv("data/group_sizes.csv")[, 1]

ages <- portfoliodata$age_at_entry
ages_max <- rep(x_max, H)

# Source thinning algorithm
source("thinning.R")

# Auxiliary function - transforms jumps & marks to init- age and E per state
get_exposure_per_state <- function(jumps_and_marks) {
  # group ids
  g_ <- unname(unlist(jumps_and_marks[1]))

  # ages at entry
  a_ <- unname(unlist(jumps_and_marks[2]))

  # jump times
  t_ <- unname(unlist(jumps_and_marks[3]))

  endpt_a <- ifelse(x_max < a_ + cov_expiry, x_max, a_ + cov_expiry)

  # if no jump times return active exposure accordingly
  if (is.null(t_)) {
    return(c(g_, a_, endpt_a - a_, 0, 0, 0))
  } else {
    # jump marks
    m_ <- unname(unlist(jumps_and_marks[4]))

    # transform jump times to (age at entry, age at jumps)
    t_ <- c(a_, a_ + t_)

    # initial mark active
    m_ <- c("Active", m_)

    # number of jump times + 1
    index <- length(t_)

    endpt_i <- ifelse(m_[index] == "Disabled", x_max, endpt_a)

    # age at next jump or censoring time
    t_endpts <- c(t_[2:index], ifelse(m_[index] == "ActiveDead" | m_[index] == "DisabledDead", t_[index], endpt_i))

    # exposures per state
    e_ <- t_endpts - t_

    # mark at right-censoring
    m_ <- c(m_, ifelse(m_[index] == "ActiveDead" | m_[index] == "DisabledDead", "Dead", "Expired"))

    # count disablements/recoveries and active/disabled deaths, respectively
    o_dis <- ifelse(m_[2:(index + 1)] == "Disabled" | m_[2:(index + 1)] == "Recovered", 1, 0)
    o_dea <- ifelse(m_[2:(index + 1)] == "ActiveDead" | m_[2:(index + 1)] == "DisabledDead", 1, 0)

    # state: active = 0, disabled = 1
    i_ <- ifelse(m_[1:index] == "Disabled" | m_[1:index] == "DisabledDead", 1, 0)

    # collect data on format: group | age | E | O_dis | O_dea | State i
    tmp_matrix <- cbind(g_, t_, e_, o_dis, o_dea, i_)

    # truncate if observed death
    if (e_[index] == 0) {
      tmp_matrix <- tmp_matrix[-index, ]
    }

    tmp_matrix
  }
}

# Auxiliary function 2 - transforms init. age and E per state to O/E-form
get_oe_inner <- function(init_age, exposure, o_dis, state, x) {
  # Indicator for initial being leq to x
  lt_age <- ifelse(init_age <= x, 1, 0)
  ct_age <- ifelse(init_age > x & init_age <= x + 1, x + 1 - init_age, lt_age)

  # Diff between age at end of observations and x if positive
  in_age <- ifelse(init_age + exposure - x < 0, 0, init_age + exposure - x)

  # Truncated age at end of observations
  rc_age <- as.integer(init_age + exposure)

  # Raw exposure
  e <- ifelse(init_age + exposure - x > 1, ct_age, in_age)

  # State-differentiated E & O
  e1 <- sum(ifelse(state == 0, e, 0))
  o1 <- sum(ifelse(rc_age == x & state == 0, o_dis, 0))
  e2 <- sum(ifelse(state == 1, e, 0))
  o2 <- sum(ifelse(rc_age == x & state == 1, o_dis, 0))

  # Return on format: x.5 | E1 | O1 | E2 | O2
  c(x + 0.5, e1, o1, e2, o2)
}

generate_oe <- function(id, theta1vec, theta2vec) {
  # Rep latent variables according to group sizes
  t1 <- rep(theta1vec, Ng[1:G])
  t2 <- rep(theta2vec, Ng[1:G])

  # Temporarily cache data s.t. thinning
  tmp <- t(mapply(thinning, ages, t1, t2, ages_max - ages, cov_expiry))

  # Right join group ids
  tmp <- cbind(rep(1:G, Ng[1:G]), tmp)

  # Reformat to initial age and E per state
  tmp <- apply(tmp, 1, get_exposure_per_state)
  tmp <- do.call(rbind, tmp)

  # Auxiliary function to generate O/E-data per group ...
  generate_oe_inner <- function(g) {
    OE_tmp <- tmp[tmp[, 1] == g, ]

    init_x <- OE_tmp[, 2]
    e <- OE_tmp[, 3]
    o <- OE_tmp[, 4]
    i <- OE_tmp[, 6]

    get_oe <- function(x) {
      get_oe_inner(init_x, e, o, i, x)
    }

    # ... and for all possible truncated ages
    OE_tmp <- lapply(x_min:(x_max - 1), get_oe)
    OE_tmp <- do.call(rbind, OE_tmp)

    OE_tmp
  }

  # Reformat to O/E
  tmp <- lapply(1:G, generate_oe_inner)
  tmp <- do.call(rbind, tmp)

  # Collect data
  oe_data <- data.frame(
    id = id,
    group = rep(1:G, each = x_max - x_min),
    age = tmp[, 1],
    e1 = tmp[, 2],
    o1 = tmp[, 3],
    e2 = tmp[, 4],
    o2 = tmp[, 5]
  )

  oe_data
}

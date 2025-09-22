##############################################
############# Thinning algorithm #############
######### Furrer, Yslas & Soerensen ##########
################### 2025 #####################
##############################################

#########################
## Markov jump process ##
#########################

# State space
# {0} : Active
# {1} : Disabled
# {2} : Dead

# Transition intensities
# Retrieved from: https://arxiv.org/pdf/2007.04051.pdf (Furrer 2021)
mu01 <- function(x) {
  exp(-4.5 - 0.018 * x + 0.00064 * x^2)
}

mu02 <- function(x, x0) {
  0.0005 + 10^(5.88 + 0.038 * x - 10)
}

mu10 <- function(x) {
  exp(0.3 - 0.049 * x)
}

mu12 <- function(x) {
  exp(-7.25 + 0.07 * x)
}

# Thinning algorithm
thinning <- function(age_at_entry, theta1, theta2, time, coverage_expiry) {
  # Initialization
  x <- age_at_entry # running age
  y <- "Active" # running mark
  x_max <- x + time # age at right-censoring (e.g. expiry)

  times <- NULL # jump times
  marks <- NULL # jump marks

  # mu_i : sum of intensities emanating from state i
  mu0 <- function(x) {
    theta1 * mu01(x) + mu02(x, age_at_entry)
  }
  mu1 <- function(x) {
    theta2 * mu10(x) + mu12(x)
  }

  # Upper bound for mu0 (monotonously increasing)
  mu0_star <- mu0(x_max)

  # Upper bound for mu1 (m10/mu12 monotonously decreasing/increasing)
  mu1_star <- theta2 * mu10(age_at_entry) + mu12(x_max)

  # Sample increment and uniform rv
  s <- rexp(1, mu0_star)
  u <- runif(1)
  x <- x + s

  while (x <= x_max & (y == "Active" | y == "Recovered") & x - age_at_entry < coverage_expiry) {
    # Acceptance-rejection
    if (u < mu0(x) / mu0_star) {
      # Sample mark
      y <- sample(c("Disabled", "ActiveDead"),
        size = 1,
        prob = c(theta1 * mu01(x) / mu0(x), mu02(x, age_at_entry) / mu0(x))
      )

      # Update times and marks
      times <- c(times, x - age_at_entry)
      marks <- c(marks, y)

      # Disablement
      if (y == "Disabled") {
        # Sample increment and uniform rv
        s <- rexp(1, mu1_star)
        u <- runif(1)
        x <- x + s

        while (x <= x_max & y == "Disabled") {
          # Acceptance-rejection
          if (u < mu1(x) / mu1_star) {
            # Sample mark
            y <- sample(c("Recovered", "DisabledDead"),
              size = 1,
              prob = c(theta2 * mu10(x) / mu1(x), mu12(x) / mu1(x))
            )

            # Update times and marks
            times <- c(times, x - age_at_entry)
            marks <- c(marks, y)
          }

          # Sample increment and uniform rv
          s <- rexp(1, mu1_star)
          u <- runif(1)
          x <- x + s
        }
      }
    }

    # Sample increment and uniform rv
    s <- rexp(1, mu0_star)
    u <- runif(1)
    x <- x + s
  }

  list(age_at_entry = age_at_entry, times = times, marks = marks)
}

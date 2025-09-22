##############################################
############# Mixing coefficients ############
######### Hierarchical gamma priors ##########
######### Furrer, Yslas & Soerensen ##########
################### 2025 #####################
##############################################

# Determine mixing coefficients w_m by recursion

get_mixing_coef <- function(min, max, p, N, gam, nu) {
  p1 <- p
  # if min(O02,O01) > 1
  if (min > 1) {
    a <- c(rep(0, N + 1))
    a[2] <- 1 # <- this is a_2
    # Jgm = 2, for m = 1,...,min(O02,O01)-1
    for (i in 1:(min - 1)) {
      a[((i - 1) * 2 + 3):2] <- i * a[((i - 1) * 2 + 3):2] + a[((i - 1) * 2 + 2):1]
      a[((i - 1) * 2 + 4):2] <- i * a[((i - 1) * 2 + 4):2] + a[((i - 1) * 2 + 3):1]
    }
    if (min != max) {
      # Jgm = 1, for m = min(O02,O01),...,max(O02,O01)-1
      for (i in (min:(max - 1))) {
        a[(min + i + 1):2] <- i * a[(min + i + 1):2] + a[(min + i):1]
      }
    }
    b <- 1
    w <- c(rep(0, N + 1))
    w[1] <- 0 # <- w0
    for (i in 1:N) {
      w[i + 1] <- (nu^gam / gamma(gam) * gamma(i + gam) / (p + nu)^(i + gam)) * a[i]
    }
    return(w)
  }

  # if min(O02,O01) <= 1 and max(O02,O01) > 1
  if (min <= 1 & max > 1) {
    a <- c(rep(0, max + 1))
    a[2] <- 1 # <- a1 if min = 0 and a2 if min = 1
    # Jgm = 1, for m = 1,...,max(O02,O01)-1
    for (i in 1:(max - 1)) {
      a[(i + 2):2] <- i * a[(i + 2):2] + a[(i + 1):1]
    }
    if (min > 0) {
      b <- 1
      w <- c(rep(0, N + 1)) # min > 0: max+1 = N
      w[1] <- 0 # <- w0
      for (i in 1:N) {
        w[i + 1] <- (nu^gam / gamma(gam) * gamma(i + gam) / (p + nu)^(i + gam)) * a[i]
      }
      return(w)
    } else {
      b <- 1
      w <- c(rep(0, N + 1)) # min = 0: max = N
      w[1] <- 0 # <- w0
      for (i in 1:N) {
        w[i + 1] <- (nu^gam / gamma(gam) * gamma(i + gam) / (p + nu)^(i + gam)) * a[i + 1]
      }
      return(w)
    }
  }

  # if min(O02,O01) <= 1 and max(O02,O01) <= 1
  else {
    ifelse(max == 0 & min == 0, return(1), ifelse(max == 1 & min == 1, return(c(0, 0, (nu^gam / gamma(gam) * gamma(2 + gam) / (p + nu)^(2 + gam)))), return(c(0, (nu^gam / gamma(gam) * gamma(1 + gam) / (p + nu)^(1 + gam))))))
  }
}

## Load library
pkgs <- c("aftgee", "survival", "geepack", "SQUAREM", "daarem")
invisible(sapply(pkgs, require, character.only = TRUE))

#' Function to generate data with a cure fraction
#'
#' @param n sample size
#' @param rho correlation
#' @param c censoring parameter
#' @param datadep correlation structure
#' @param para1 regression coefficients in the incidence component
#' @param para2 regression coefficients in the latency component
#' 
dat_with_cure <- function(n = 200, rho = 0, c = 0.005,
                          datadep = "ind",
                          para1 = c(0.3, 0.4, -0.2, -0.3),
                          para2 = c(1, 0.6, 0.5, -0.1)) {
  dat <- NULL
  for (i in 1:n) {
    alpha <- rnorm(1, mean = 0, sd = sqrt(0.05))
    cluster_size <- rpois(1, exp(1 + 2 * alpha)) + 1
    x1 <- runif(cluster_size, min = 0, max = 2)
    x2 <- rbinom(cluster_size, 1, 0.5)
    x3 <- rnorm(cluster_size, mean = 1, sd = sqrt(0.3))
    x4 <- rnorm(cluster_size, mean = -1, sd = 1)
    if (datadep == "ind") {
      e <- rnorm(cluster_size, mean = 0, sd = sqrt(0.5))
    }
    if (datadep == "ex") {
      e <- mvrnorm(1,
                   mu = rep(0, cluster_size),
                   Sigma = 0.5 * (diag(1 - rho, cluster_size) +
                                  matrix(rho, cluster_size, cluster_size)))
    }
    if (datadep == "ar1") {
      e <- mvrnorm(1,
                   mu = rep(0, cluster_size),
                   Sigma = 0.5 * (outer(
                     1:cluster_size, 1:cluster_size,
                     function(x, y) rho^abs(x - y)
                   )))
    }
    tij <- exp(2 + para1[1] * x1 + para1[2] * x2 + para1[3] * x3 + para1[4] * x4 + alpha + e)
    cen <- rexp(n = cluster_size, rate = c)
    eta <- rbinom( ## 0 if cure
      n = cluster_size, size = 1,
      prob = 1 / (1 + 1 / exp(para2[1] * x1 + para2[2] * x2 + para2[3] * x3 + para2[4] * x4))
    )
    dat[[i]] <- data.frame(
      id = rep(i, cluster_size),
      yij = pmin(ifelse(eta == 0, Inf, tij), cen),
      delta = 1 * (ifelse(eta == 0, Inf, tij) <= cen),
      x1 = x1,
      x2 = x2,
      x3 = x3,
      x4 = x4,
      weight = rep(1 / cluster_size, cluster_size)
    )
  }
  dat <- do.call(rbind, dat)
  rownames(dat) <- NULL
  return(dat)
}


#' Function to compute survival function from an AFT fit
#'
#' @param y is the response variable created by Surv()
#' @param x is the covariate matrix
#' @param b is the regression coefficient, estimated by regression x to y in an AFT model.
#' @param weight1 is the sampling weight for informative cluster size
#' @param weight2 is psi, the EM weight
#' 
aftsurv <- function(y, x, b, weight1, weight2) {
  err <- c(log(y[, 1]) - x %*% b)
  km <- survfit(Surv(err, y[, 2]) ~ 1, weights = weight1)
  km2 <- survfit(Surv(err, y[, 2]) ~ 1, weights = weight1 * weight2)
  km$surv <- exp(-cumsum(km$n.event / km2$n.risk)) # Nelson Aalen
  # km$surv <- cumprod(1 - km$n.event / km2$n.risk) # Kaplan Meier
  km$surv[km$time > max(err[y[, 2] > 0])] <- 0
  km2 <- summary(km, err, extend = TRUE)
  return(km2$surv[match(err, km2$time)])
}

#' Function to carry out one iteration of the EM algorithm 
#'
#' @param para is a vector of all parameter of interest;
#' this includes the regression coefficient in both components and
#' the survival estiamtes for the non-cure proportion.
#' @param x is the covariate matrix for the latency component
#' @param z is the covariate matrix for the incident component
#' @param y is the observed survival time created by Surv()
#' @param id is the cluster id
#' @param weight is the sampling weight for informative cluster size
#' @param corstr is the correlation structure
#' @param formula specifies the formula for the AFT model in the latency component
#' @param data optional data frame 
one <- function(para, x, z, y, id, weight, corstr, formula, data) {
  a <- para[seq_len(NCOL(z))]
  b <- para[seq_len(NCOL(z)) + NCOL(x)]
  surv_prob <- para[-seq_len(NCOL(z) + NCOL(x))]
  pi <- 1 / (1 + c(exp(-z %*% a)))
  psi <- ifelse(y[, 2] == 1, 1, pi * surv_prob / (1 - pi + pi * surv_prob))
  fit_a <- geese.fit(
             y = psi, x = z, id = id, corstr = corstr,
             weights = weight,
             family = quasibinomial(link = "logit"),
             variance = "binomial"
           )
  data$weight_new <- psi * weight
  fit_b <- aftgee(formula,
                  id = id, B = 0, weights = weight_new, # nolint: object_usage_linter.
                  corstr = corstr, data = data, binit = b
                  )
  return(c(
    fit_a$beta, coef(fit_b),
    aftsurv(y = y, x = x, b = coef(fit_b), weight1 = weight, weight2 = psi)
  ))
}

#' Function to run the whole EM algorithm 
#'
#' @param xformula model formula for the latency component
#' @param zformula model formula for the incident component
#' @param data optional data frame 
#' @param id is the cluster id
#' @param weight is the sampling weight for informative cluster size
#' @param corstr is the correlation structure
#' @param em which EM accelerator to use? squarem or daarem?
cure_em <- function(xformula, zformula, data, id, weight,
                    corstr = c("ind", "ex", "ar1"),
                    em = c("squarem", "daarem")) {
  call <- match.call()
  corstr <- match.arg(corstr)
  em <- match.arg(em)
  dat <- data[complete.cases(data), ]
  mf <- model.frame(xformula, dat)
  y <- model.extract(mf, "response")
  x <- model.matrix(xformula, dat)
  z <- model.matrix(zformula, dat)
  id <- eval(call["id"][[1]], dat)
  if (is.null(id)) id <- seq_len(nrow(dat))
  weight <- eval(call["weight"][[1]], dat)
  if (is.null(weight)) weight <- rep(1, nrow(dat))
  fit_a <- geese.fit(
             y = y[, 2], x = z, id = id, corstr = corstr,
             family = binomial(link = "logit"), variance = "binomial"
           )
  fit_b <- aftgee(
             formula = xformula, id = id, weights = weight,
             B = 0, corstr = corstr, data = dat, binit = "lm"
           )
  if (em == "squarem") em_fun <- squarem
  if (em == "daarem") em_fun <- daarem
  emfit <- em_fun(
             par = c(
                     fit_a$beta, coef(fit_b),
                     aftsurv(y, x, coef(fit_b), weight, rep(1, length(y)))
                   ),
             x = x, z = z, y = y, id = id, weight = weight,
             corstr = corstr, formula = xformula, data = dat,
             fixptfn = one,
             control = list(maxiter = 20, tol = 1e-7)
           )
  return(list(
    coef_cure = emfit$par[seq_len(NCOL(z))],
    coef_aft = emfit$par[seq_len(NCOL(z)) + NCOL(x)],
    surv_prob = emfit$par[-seq_len(NCOL(z) + NCOL(x))],
    convergence = emfit$convergence,
    iter = emfit$iter,
    fpevals = emfit$fpevals,
    call = call
  ))
}

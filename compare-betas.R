require(R2jags)

summary <- function(x, name, n.sample, prob, digits) {
  if (bugs) {
    y <- x$sims.list[[name]]
  } else
    y <- x$BUGSoutput$sims.list[[name]]
  }
  cat(name, ": ", sep="")
  non.repens <- sample(y[,1], n.sample, replace=TRUE)
  repens <- sample(y[,2], n.sample, replace=TRUE)
  diff <- repens-non.repens
  cat(round(mean(diff), digits), " ", sep="")
  hpd <- HPDinterval(as.mcmc(diff), prob=prob)
  cat("(", round(hpd[,"lower"], digits), ",",
      round(hpd[,"upper"], digits), ")", sep="")
  if ((hpd[,"lower"] > 0) || (hpd[,"upper"] < 0)) {
    cat("*")
  }
  cat("\n")
}

compare.env.coefficients <- function(x, n.sample, prob=0.95, digits=3) {
  env <- c("elev", "long", "map")
  cov <- c("fecundity", "fl.length", "fl.per.head", "seed.mass")
  for (trait in cov) {
    for (e in env) {
      name <- paste("beta.", trait, ".", e, sep="")
      summary(x, name, n.sample, prob, digits)
    }
  }
}

compare.infest.coefficients <- function(x, n.sample, prob=0.95, digits=3) {
  cov <- c("elev", "fecundity", "fl.length", "fl.per.head", "long", "map",
           "seed.mass")
  for (trait in cov) {
    name <- paste("beta.infest.", trait, sep="")
    summary(x, name, n.sample, prob, digits)
  }
}

compare.pink.coefficients <- function(x, n.sample, prob=0.95, digits=3) {
  cov <- c("elev", "fecundity", "fl.length", "fl.per.head", "infest", "long",
           "map", "seed.mass")
  for (trait in cov) {
    name <- paste("beta.pink.", trait, sep="")
    if (bugs) {
      y <- x$sims.list[[name]]
    } else
      y <- x$BUGSoutput$sims.list[[name]]
    }
    for (i in 1:4) {
      cat(name, "[", i, "]: ", sep="")
      non.repens <- sample(y[,i,1], n.sample, replace=TRUE)
      repens <- sample(y[,i,2], n.sample, replace=TRUE)
      diff <- repens-non.repens
      cat(round(mean(diff), digits), " ", sep="")
      hpd <- HPDinterval(as.mcmc(diff), prob=prob)
      cat("(", round(hpd[,"lower"], digits), ",",
          round(hpd[,"upper"], digits), ")", sep="")
      if ((hpd[,"lower"] > 0) || (hpd[,"upper"] < 0)) {
          cat("*")
      }
      cat("\n")
    }
  }
}

compare.poly.coefficients <- function(x, n.sample, prob=0.95, digits=3) {
  cov <- c("elev", "fecundity", "fl.length", "fl.per.head", "infest", "long",
           "map", "seed.mass")
  for (trait in cov) {
    name <- paste("beta.poly.", trait, sep="")
    if (bugs) {
      y <- x$sims.list[[name]]
    } else {
      y <- x$BUGSoutput$sims.list[[name]]
    }
    for (i in 1:2) {
      cat(name, "[", i, "]: ", sep="")
      non.repens <- sample(y[,i,1], n.sample, replace=TRUE)
      repens <- sample(y[,i,2], n.sample, replace=TRUE)
      diff <- repens-non.repens
      cat(round(mean(diff), digits), " ", sep="")
      hpd <- HPDinterval(as.mcmc(diff), prob=prob)
      cat("(", round(hpd[,"lower"], digits), ",",
          round(hpd[,"upper"], digits), ")", sep="")
      if ((hpd[,"lower"] > 0) || (hpd[,"upper"] < 0)) {
          cat("*")
      }
      cat("\n")
    }
  }
}

compare <- function(x, prob=0.95) {
  cat("Posterior comparisons of regression coefficients...\n",
      "   repens - non.repens\n\n")
  flush.console()
  cat("Environmental coefficients...\n")
  flush.console()
  cat(compare.env.coefficients(x, 25000, prob), "\n")
  flush.console()
  cat("Infestation coefficients...\n")
  flush.console()
  cat(compare.infest.coefficients(x, 25000, prob), "\n")
  flush.console()
  cat("Pink coefficients...\n")
  flush.console()
  cat(compare.pink.coefficients(x, 25000, prob), "\n")
  cat("Poly coefficients...\n")
  flush.console()
  cat(compare.poly.coefficients(x, 25000, prob), "\n")
}



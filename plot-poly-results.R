require(MCMCglmm)
require(ggplot2)
require(reshape2)

rm(list=ls())

list.mean <- function(x) {
  y <- numeric(0)
  for (i in 1:length(x)) {
    y <- c(y, x[[i]])
  }
  mean(y)
}

inv.logit <- function(x) {
  p <- exp(x)/(1 + exp(x))
  p
}

get.probability <- function(scaled, beta, var, x, species) {
  ## calculate probabilities
  ##
  ## non-varying covariates are incorporated into intercept term b.0
  ## otherwise set slope b.1
  ##
  pi <- numeric(length(x))
  for (i in 1:length(x)) {
    b.0 <- beta$beta.0 + species[grep(scaled$species, names(species))]
    for (v in c("ele",
                "ele.range",
                "long",
                "map",
                "map.range",
                "area")) {
      beta.var <- paste("beta.", v, sep="")
      if (v == var) {
        b.1 <- beta[[beta.var]]
      } else {
        b.0 <- b.0 + beta[[beta.var]]*scaled[[v]]
      }
    }
    mu <- b.0 + b.1*x[i]
    pi[i] <- inv.logit(mu)
  }
  pi
}

plot.prediction <- function(beta, scaled, unscaled, var, species) {
  ## set minima and maxima for plots
  ##
  min <- min(scaled[[var]])
  max <- max(scaled[[var]])
  min.unscaled <- min(unscaled[[var]])
  max.unscaled <- max(unscaled[[var]])
  ## set up covariate to plot
  ##
  n.pts <- 1000
  x <- seq(from=min, to=max, by=(max-min)/n.pts)

  ## calculate probabilities for each species (ignoring species
  ## random effect)
  ##
  pi.k <- matrix(nrow=nrow(scaled), ncol=length(x))
  for (k in 1:nrow(scaled)) {
    pi.k[k,] <- get.probability(scaled[k,], beta, var, x, species)
  }
  ## average across species
  #
  pi <- numeric(0)
  for (j in 1:length(x)) {
    pi[j] <- mean(pi.k[,j])
  }
  
  x.unscaled <- seq(from=min.unscaled, to=max.unscaled,
                    by=(max.unscaled-min.unscaled)/n.pts)

  ret.val <- data.frame(x=x.unscaled,
                        y=pi)
  ret.val
}

## x - an mcmc.list returned by MCMCglmm
## scaled - scaled covariates
## unscaled - unscaled covariates
##
plot.poly.results <- function(x, scaled, unscaled, pdf=TRUE) {
  ## concatenate MCMC results into single vector
  ##
  beta.0 <- numeric(0)
  beta.ele <- numeric(0)
  beta.ele.range <- numeric(0)
  beta.long <- numeric(0)
  beta.map <- numeric(0)
  beta.map.range <- numeric(0)
  beta.area <- numeric(0)
  for (i in 1:length(x)) {
    beta.0 <- c(beta.0, x[[i]][,1])
    beta.ele <- c(beta.ele, x[[i]][,2])
    beta.ele.range <- c(beta.ele.range, x[[i]][,3])
    beta.long <- c(beta.long, x[[i]][,4])
    beta.map <- c(beta.map, x[[i]][,5])
    beta.map.range <- c(beta.map.range, x[[i]][,6])
    beta.area <- c(beta.area, x[[i]][,7])
  }
  ## get posterior mean for all coefficients in the model and
  ## put them in a data frame
  ##
  beta <- data.frame(beta.0=mean(beta.0),
                     beta.ele=mean(beta.ele),
                     beta.ele.range=mean(beta.ele.range),
                     beta.long=mean(beta.long),
                     beta.map=mean(beta.map),
                     beta.map.range=mean(beta.map.range),
                     beta.area=mean(beta.area))
  ## get posterior mean of species random effects
  ##
  idx <- grep("Pr", colnames(x[[1]]))
  species <- numeric(0)
  for (i in idx) {
    species[i] <- list.mean(x[,i])
  }
  species.names <- substring(colnames(results.mcmc[[1]])[idx], 8)
  names(species) <- c(rep("",29), species.names)
  ## loop through each covariate and produce a plot
  ##
  idx <- grep("animal", colnames(results.mcmc[[1]]), invert=TRUE)[-1]
  var <- colnames(x[[1]])[idx]
  for (v in var) {
    if (v == "ele") {
      label <- "Elevation"
    } else if (v == "ele.range") {
      label <- "Elevation range"
    } else if (v == "long") {
      label <- "Longitude"
    } else if (v == "map") {
      label <- "Mean annual precipitation"
    } else if (v == "map.range") {
      label <- "Mean annual precipitation range"
    } else if (v == "area") {
      label <- "Geographical range area"
    } else {
      stop("Should never get here")
    }
    comp <- plot.prediction(beta,
                            scaled,
                            unscaled,
                            v,
                            species)
    tmp.yend <- numeric(0)
    for (i in 1:nrow(unscaled)) {
      if (unscaled$poly[i] == 0) {
        tmp.yend[i] <- -0.025
      } else if (unscaled$poly[i] == 1) {
        tmp.yend[i] <- 1.025
      } else {
        stop("Should never get here")
      }
    }
    rug <- data.frame(x=unscaled[[v]],
                      xend=unscaled[[v]],
                      y=unscaled$poly,
                      yend=tmp.yend,
                      variable=v)
    dev.new()
    p <- ggplot(comp, aes(x=x, y=y)) +
         geom_line() +
         xlab(label) +
         ylab("Polymorphism probability") +
         geom_segment(mapping=aes(x=x, xend=xend, y=y, yend=yend), data=rug)
    print(p)
    if (pdf) {
      ggsave(file=paste("poly-MCMC-", v, ".pdf", sep=""))
    }
  }
}

load("poly-results.Rsave")

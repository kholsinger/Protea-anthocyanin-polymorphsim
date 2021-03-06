require(ggplot2)
require(reshape2)

debug <- FALSE

get.probabilities <- function(work, beta, var, x, verbose=FALSE) {
  ## set up the necessary matrices and vectors
  ##
  n.cat <- nrow(beta) + 1
  l <- matrix(nrow=n.cat-1, ncol=length(x))
  l.star <- matrix(nrow=n.cat-1, ncol=length(x))
  phi <- matrix(nrow=n.cat, ncol=length(x))
  pi <- matrix(nrow=n.cat, ncol=length(x))
  b.0 <- numeric(n.cat-1)
  b.1 <- numeric(n.cat-1)

  ## calculate probabilities
  ##
  ## non-varying covariates are incorporated into intercept term b.0
  ## otherwise set slope b.1
  ##
  for (i in 1:(n.cat-1)) {
    b.0[i] <- beta$beta.0[i]
    if (verbose) {
      cat("b.0[", i, "]: ", b.0[i], "\n", sep="")
    }
    for (v in c("fecundity", "fl.per.head", "infest", "map", "elev", "long",
                "seed.mass")) {
      beta.var <- paste("beta.", v, sep="")
      if (verbose) {
        cat("v: ", v, "\n", sep="")
      }
      if (v == var) {
        b.1[i] <- beta[[beta.var]][i]
        if (verbose) {
          cat("b.1[", i, "]: ", b.1[i], "\n", sep="")
        }
      } else {
        b.0[i] <- b.0[i] + beta[[beta.var]][i]*work[[v]]
        if (verbose) {
          cat("b.0[", i, "]: ", b.0[i], "\n", sep="")
        }
      }
    }
    for (j in 1:length(x)) {
      l[i,j] <- b.0[i] + b.1[i]*x[j]
    }
  }
  ## rows (i) are categories
  ## columns (j) are level of covariate being plotted
  ##
  for (i in 1:(n.cat-1)) {
    for (j in 1:length(x)) {
      l.star[i,j] <- sum(l[i:(n.cat-1),j])
      phi[i,j] <- exp(l.star[i,j])
    }
  }
  phi[n.cat,] <- rep(1.0, length(x))
  for (j in 1:length(x)) {
    sum.phi <- sum(phi[,j])
    for (i in 1:n.cat) {
      pi[i,j] <- phi[i,j]/sum.phi
    }
  }

  pi
}

plot.prediction <- function(beta, species, var) {
  ## set minima and maxima for plots
  ##
  min <- min(as.numeric(color[[var]]))
  max <- max(as.numeric(color[[var]]))
  min.unscaled <- min(unscaled[[var]])
  max.unscaled <- max(unscaled[[var]])
  ## set up covariate to plot
  ##
  n.cat <- nrow(beta) + 1
  n.pts <- 1000
  x <- seq(from=min, to=max, by=(max-min)/n.pts)

  ## get site-specific covariates for species
  ##
  work <- color[color$species==species,]

  ## calculate piulative probabilities for each site
  ##
  pi.k <- array(dim=c(nrow(work), n.cat, length(x)))
  for (k in 1:nrow(work)) {
    if (debug && (species == "repens") && (var == "infest")) {
      verbose <- TRUE
      cat("species: ", species, "\n", sep="")
      cat("work[k,]:\n")
      print(work[k,])
      cat("beta:\n")
      print(beta)
      cat("var: ", var, "\n", sep="")
    } else {
      verbose <- FALSE
    }
    pi.k[k,,] <- get.probabilities(work[k,], beta, var, x, verbose)
  }
  ## average across sites
  pi <- matrix(nrow=n.cat, ncol=length(x))
  for (i in 1:n.cat) {
    for(j in 1:length(x)) {
      pi[i,j] <- mean(pi.k[,i,j])
    }
  }

  if (n.cat==5) {
    rownames(pi) <- c("white", "skewed white", "even", "skewed pink",
                       "pink")
  } else {
    rownames(pi) <- c("monomorphic", "skewed", "polymorphic")
  }
  colnames(pi) <- x
  pi <- data.frame(t(pi))

  x.unscaled <- seq(from=min.unscaled, to=max.unscaled,
                    by=(max.unscaled-min.unscaled)/n.pts)
  pi$species <- species
  pi[[var]] <- x.unscaled
  one.species <- melt(pi, id=c("species", var))
  one.species
}

extract.pink.coefficients <- function(x) {
  beta.0 <- x$BUGSoutput$sims.list$beta.pink.0
  beta.fecundity <- x$BUGSoutput$sims.list$beta.pink.fecundity
  beta.fl.per.head <- x$BUGSoutput$sims.list$beta.pink.fl.per.head
  beta.infest <- x$BUGSoutput$sims.list$beta.pink.infest
  beta.map <- x$BUGSoutput$sims.list$beta.pink.map
  beta.elev <- x$BUGSoutput$sims.list$beta.pink.elev
  beta.long <- x$BUGSoutput$sims.list$beta.pink.long
  beta.seed.mass <- x$BUGSoutput$sims.list$beta.pink.seed.mass
  list(beta.0=beta.0, beta.fecundity=beta.fecundity,
       beta.fl.per.head=beta.fl.per.head,
       beta.infest=beta.infest, beta.map=beta.map, beta.elev=beta.elev,
       beta.long=beta.long, beta.seed.mass=beta.seed.mass)
}

extract.poly.coefficients <- function(x) {
  beta.0 <- x$BUGSoutput$sims.list$beta.poly.0
  beta.fecundity <- x$BUGSoutput$sims.list$beta.poly.fecundity
  beta.fl.per.head <- x$BUGSoutput$sims.list$beta.poly.fl.per.head
  beta.infest <- x$BUGSoutput$sims.list$beta.poly.infest
  beta.map <- x$BUGSoutput$sims.list$beta.poly.map
  beta.elev <- x$BUGSoutput$sims.list$beta.poly.elev
  beta.long <- x$BUGSoutput$sims.list$beta.poly.long
  beta.seed.mass <- x$BUGSoutput$sims.list$beta.poly.seed.mass
  list(beta.0=beta.0, beta.fecundity=beta.fecundity,
       beta.fl.per.head=beta.fl.per.head,
       beta.infest=beta.infest, beta.map=beta.map, beta.elev=beta.elev,
       beta.long=beta.long, beta.seed.mass=beta.seed.mass)
}

plot.var <- function(x, var, model) {
##  old <- par(mfrow=c(2,2))
  comp <- data.frame(species=NULL,
                     variable=NULL,
                     value=NULL)
  comp[[var]] <- NULL
  for (i in 1:4) {
    species.label <- as.character(with(species.table, species[number==i]))
    ## set repens/non-repens idx
    ##
    if (species.label == "repens") {
        idx <- 2
    } else {
        idx <- 1
    }
    beta <- data.frame(beta.0=apply(x$beta.0[,,i], 2, mean),
                       beta.fecundity=apply(x$beta.fecundity[,,idx], 2, mean),
                       beta.fl.per.head=apply(x$beta.fl.per.head[,,idx], 2, mean),
                       beta.infest=apply(x$beta.infest[,,idx], 2, mean),
                       beta.map=apply(x$beta.map[,,idx], 2, mean),
                       beta.elev=apply(x$beta.elev[,,idx], 2, mean),
                       beta.long=apply(x$beta.long[,,idx], 2, mean),
                       beta.seed.mass=apply(x$beta.seed.mass[,,idx], 2, mean))
    temp <- plot.prediction(beta,
                            species.label,
                            var)
    comp <- rbind(comp, temp)
  }
  if (model == "pink") {
    pal <- "Reds"
  } else if (model == "poly") {
    pal <- "Greys"
  } else {
    cat("Unrecognized model: ", model, "\n", sep="")
    cat("  Model must be \"pink\" or \"poly\" (case sensitive)\n")
    stop(call.=FALSE)
  }
  temp <- data.frame(x=unscaled[[var]],
                     xend=unscaled[[var]],
                     y=0.0,
                     yend=-0.025,
                     species=unscaled$species,
                     variable="white")
  if (var == "elev") {
    label <- "Elevation in meters"
    p <- ggplot(comp, aes(x=elev, y=value, fill=variable)) +
         geom_area(colour="black", size=0.2, alpha=0.4) +
         scale_fill_brewer(palette=pal,
                           breaks=rev(levels(comp$variable))) +
         xlab(label) +
         ylab("Cumulative probability") +
         labs(fill="Color category") +
         geom_segment(mapping=aes(x=x, xend=xend, y=y, yend=yend),
                      data=temp) +
         facet_wrap(~ species)
  } else if (var == "fecundity") {
    label <- "Fecundity"
    p <- ggplot(comp, aes(x=fecundity, y=value, fill=variable)) +
         geom_area(colour="black", size=0.2, alpha=0.4) +
         scale_fill_brewer(palette=pal,
                           breaks=rev(levels(comp$variable))) +
         xlab(label) +
         ylab("Cumulative probability") +
         labs(fill="Color category") +
         geom_segment(mapping=aes(x=x, xend=xend, y=y, yend=yend),
                      data=temp) +
         facet_wrap(~ species)
  } else if (var == "fl.per.head") {
    label <- "Flowers per head"
    p <- ggplot(comp, aes(x=fl.per.head, y=value, fill=variable)) +
         geom_area(colour="black", size=0.2, alpha=0.4) +
         scale_fill_brewer(palette=pal,
                           breaks=rev(levels(comp$variable))) +
         xlab(label) +
         ylab("Cumulative probability") +
         labs(fill="Color category") +
         geom_segment(mapping=aes(x=x, xend=xend, y=y, yend=yend),
                      data=temp) +
         facet_wrap(~ species)
  } else if (var == "infest") {
    label <- "Fraction of heads infested"
    p <- ggplot(comp, aes(x=infest, y=value, fill=variable)) +
         geom_area(colour="black", size=0.2, alpha=0.4) +
         scale_fill_brewer(palette=pal,
                           breaks=rev(levels(comp$variable))) +
         xlab(label) +
         ylab("Cumulative probability") +
         labs(fill="Color category") +
         geom_segment(mapping=aes(x=x, xend=xend, y=y, yend=yend),
                      data=temp) +
         facet_wrap(~ species)
  } else if (var == "long") {
    label <- "East longitude"
    p <- ggplot(comp, aes(x=long, y=value, fill=variable)) +
         geom_area(colour="black", size=0.2, alpha=0.4) +
         scale_fill_brewer(palette=pal,
                           breaks=rev(levels(comp$variable))) +
         xlab(label) +
         ylab("Cumulative probability") +
         labs(fill="Color category") +
         geom_segment(mapping=aes(x=x, xend=xend, y=y, yend=yend),
                      data=temp) +
         facet_wrap(~ species)
  } else if (var == "map") {
    label <- "Mean annual precipitation"
    p <- ggplot(comp, aes(x=map, y=value, fill=variable)) +
         geom_area(colour="black", size=0.2, alpha=0.4) +
         scale_fill_brewer(palette=pal,
                           breaks=rev(levels(comp$variable))) +
         xlab(label) +
         ylab("Cumulative probability") +
         labs(fill="Color category") +
         geom_segment(mapping=aes(x=x, xend=xend, y=y, yend=yend),
                      data=temp) +
         facet_wrap(~ species)
  } else if (var == "seed.mass") {
    label <- "Seed mass"
    p <- ggplot(comp, aes(x=seed.mass, y=value, fill=variable)) +
         geom_area(colour="black", size=0.2, alpha=0.4) +
         scale_fill_brewer(palette=pal,
                           breaks=rev(levels(comp$variable))) +
         xlab(label) +
         ylab("Cumulative probability") +
         labs(fill="Color category") +
         geom_segment(mapping=aes(x=x, xend=xend, y=y, yend=yend),
                      data=temp) +
         facet_wrap(~ species)
  }
  print(p)
}

## x - an object returned by R2jags
##
plot.pink.results <- function(x, pdf=TRUE) {
  var <- c("elev",
           "infest",
           "long",
           "fl.per.head",
           "seed.mass")
  for (v in var) {
    dev.new()
    plot.var(extract.pink.coefficients(x), v, "pink")
    if (pdf) {
      ggsave(file=paste("pink-", v, ".pdf", sep=""))
    }
  }
}

## x - an object returned by R2jags
##
plot.poly.results <- function(x, pdf=TRUE) {
  var <- c("infest",
           "fl.per.head",
           "long",
           "seed.mass")
  for (v in var) {
    dev.new()
    plot.var(extract.poly.coefficients(x), v, "poly")
    if (pdf) {
      ggsave(file=paste("poly-", v, ".pdf", sep=""))
    }
  }
}

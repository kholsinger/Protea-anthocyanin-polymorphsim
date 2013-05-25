require(plotrix)

get.probabilities <- function(work, beta, var, x) {
  ## set up the necessary matrices and vectors
  ##
  n.cat <- nrow(beta) + 1
  l <- matrix(nrow=n.cat-1, ncol=length(x))
  l.star <- matrix(nrow=n.cat-1, ncol=length(x))
  phi <- matrix(nrow=n.cat, ncol=length(x))
  pi <- matrix(nrow=n.cat, ncol=length(x))
  cum <- matrix(nrow=n.cat, ncol=length(x))
  b.0 <- numeric(n.cat-1)
  b.1 <- numeric(n.cat-1)

  ## calculate probabilities
  ##
  ## rows (i) are categories
  ## columns (j) are level of covariate being plotted
  ##
  for (i in 1:(n.cat-1)) {
    b.0 <- 0.0
    for (v in c("fecundity", "fl.per.head", "infest", "map", "elev", "long",
                "seed.mass")) {
      beta.var <- paste("beta.", v, sep="")
      if (v == var) {
        b.1[i] <- beta[[beta.var]][i]
      } else {
        b.0[i] <- beta[[beta.var]][i]*work[[var]]
      }
    }
    for (j in 1:length(x)) {
      l[i,j] <- b.0[i] + b.1[i]*x[j]
    }
  }
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

  ## calculate cumulative probabilities
  ##
  for (j in 1:length(x)) {
    cum[1,j] <- pi[1,j]
    for (i in 2:n.cat) {
      cum[i,j] <- pi[i,j] + cum[i-1,j]
    }
  }

  ## return cumulative probabilities
  ##
  cum
}

get.xaxlab <- function(min, max, x) {
  val <- min + (max - min)*x
  lab <- as.character(signif(val[1], digits=2))
  for (i in 2:length(x)) {
    lab <- c(lab, as.character(signif(val[i], digits=2)))
  }
  lab
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
  
  ## calculate cumulative probabilities for each site
  ##
  cum.k <- array(dim=c(nrow(work), n.cat, length(x)))
  for (k in 1:nrow(work)) {
    cum.k[k,,] <- get.probabilities(work[k,], beta, var, x)
  }
  ## average across sites
  cum <- matrix(nrow=n.cat, ncol=length(x))
  for (i in 1:n.cat) {
    for(j in 1:length(x)) {
      cum[i,j] <- mean(cum.k[,i,j])
    }
  }

  if (n.cat==5) {
    rownames(cum) <- c("white", "skewed white", "even", "skewed pink",
                       "pink")
  } else {
    rownames(cum) <- c("monomorphic", "skewed", "polymorphic")
  }
  colnames(cum) <- x

  cum <- data.frame(t(cum))
  xat <- seq(from=n.pts/4, to=3*n.pts/4, by=n.pts/4)
  xaxlab <- get.xaxlab(min.unscaled, max.unscaled, c(0.25, 0.50, 0.75))
  colors <- rgb(red=c(255, 255, 255, 255, 255),
                green=c(255, 207, 159, 111, 63),
                blue=c(255, 207, 159, 111, 63),
                maxColorValue=255)
  stackpoly(cum, col=colors, axis4=FALSE, main=species,
            xat=xat, xaxlab=xaxlab, xlab=var)
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
  list(beta.0=beta.0, beta.fecundity=beta.fecundity, beta.fl.per.head=beta.fl.per.head,
       beta.infest=beta.infest, beta.map=beta.map, beta.elev=beta.elev,
       beta.long=beta.long, beta.seed.mass=beta.seed.mass)
}

plot.var <- function(x, var) {
  old <- par(mfrow=c(2,2))
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
    plot.prediction(beta,
                    species.label,
                    var)
  }
  par(old)
}

## x - an object returned by R2jags
##
plot.results <- function(x) {
  var <- c("elev",
           "fecundity",
           "fl.per.head",
           "infest",
           "long",
           "seed.mass")
  for (v in var) {
    ## dev.new()
    ## plot.var(extract.pink.coefficients(x), v)
    pdf(file=paste(v, ".pdf", sep=""))
    plot.var(extract.pink.coefficients(x), v)
    dev.off()
  }
}

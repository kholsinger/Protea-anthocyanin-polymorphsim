require(plotrix)

plot.prediction <- function(b.0, b.1, min, max, label, varname) {
  ## set up matrices
  ##
  n.cat <- length(b.0) + 1
  n.pts <- 1000
  x <- seq(from=min, to=max, by=(max-min)/n.pts)
  l <- matrix(nrow=n.cat-1, ncol=length(x))
  l.star <- matrix(nrow=n.cat-1, ncol=length(x))
  phi <- matrix(nrow=n.cat, ncol=length(x))
  pi <- matrix(nrow=n.cat, ncol=length(x))
  cum <- matrix(nrow=n.cat, ncol=length(x))

  ## calculate probabilities
  ##
  ## rows are categories
  ## columns are level of covariate
  ##
  for (i in 1:(n.cat-1)) {
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
  if (n.cat==5) {
    rownames(cum) <- c("white", "skewed white", "even", "skewed pink",
                       "pink")
  } else {
    rownames(cum) <- c("monomorphic", "skewed", "polymorphic")
  }
  colnames(cum) <- x

  cum <- data.frame(t(cum))
  xat <- seq(from=n.pts/4, to=3*n.pts/4, by=n.pts/4)
  if (label == "lacticolor") {
    main.label <- paste(varname, "\n\n", label, sep="")
  } else {
    main.label <- label
  }
  if (n.cat==3) {
    val <- seq(from=5, to=255, by=250/(n.cat))/255
    stackpoly(cum, col=rgb(val, val, val), axis4=FALSE, main=main.label,
              xat=xat)
  } else {
    colors <- rgb(red=c(255, 255, 255, 255, 255),
                  green=c(255, 207, 159, 111, 63),
                  blue=c(255, 207, 159, 111, 63),
                  maxColorValue=255)
    stackpoly(cum, col=colors, axis4=FALSE, main=main.label,
              xat=xat)
  }
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
  par(mfrow=c(2,3))
  for (i in 1:5) {
    species.label <- as.character(with(species.table, species[number==i]))
    ## set repens/non-repens idx
    ##
    if (species.label == "repens") {
        idx <- 2
    } else {
        idx <- 1
    }
    flush.console()
    x.var <- paste("beta.", v, sep="")
    y <- x[[x.var]]
    plot.prediction(apply(x$beta.0[,,i], 2, mean),
                    apply(y[,,idx], 2, mean),
                    min(as.numeric(color[[var]])),
                    max(as.numeric(color[[var]])),
                    species.label,
                    var)
  }
}

plot.results <- function(x) {
  var <- c("elev",
           "fecundity",
           "fl.per.head",
           "infest",
           "long",
           "seed.mass")
  for (v in var) {
    dev.new()
    plot.var(x, v)
  }
}

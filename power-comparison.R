require(R2jags)

rm(list=ls())
options("width"=120)

ptm <- proc.time()

debug <- FALSE
betas.at.zero <- TRUE

## to seed RNG
##
runif(1)

n.cat <- 5
sample.size <- 20
n.simulations <- 1000
model.file <- "power-analysis.txt"
model.file.logistic <- "power-analysis-logistic.txt"

if (debug) {
  n.chains <- 1
  n.burnin <- 1000
  n.iter <- 6000
  n.thin <- 1
} else {
  n.chains <- 5
  n.burnin <- 5000
  n.iter <- 30000
  n.thin <- 5
}

drop.levels <- function(dat) {
  dat[] <- lapply(dat, function(x) x[,drop=TRUE])
  return(dat)
}

diversity.to.poly <- function(x) {
  y <- 2*x + 1
  y
}

get.pink.scale <- function(x, n.cat=7) {
  y <- numeric(n.cat)
  if (n.cat == 7) {
    y[x=="white"] <- 1
    y[x=="less_than_ten_pink"] <- 2
    y[x=="less_than_25_pink"] <- 3
    y[x=="equal"] <- 4
    y[x=="less_than_25_white"] <- 5
    y[x=="less_than_ten_white"] <- 6
    y[x=="pink"] <- 7
  } else if (n.cat == 5){
    y[x=="white"] <- 1
    y[x=="less_than_ten_pink"] <- 2
    y[x=="less_than_25_pink"] <- 3
    y[x=="equal"] <- 3
    y[x=="less_than_25_white"] <- 3
    y[x=="less_than_ten_white"] <- 4
    y[x=="pink"] <- 5
  }
  y
}

inv.logit <- function(phi) {
  exp(phi)/(1.0 + exp(phi))
}

predict.pink <- function(species, repens, long, map, elev, fecundity,
                         fl.per.head, seed.mass, infest,
                         beta.0, beta.long, beta.map, beta.elev, beta.fecundity,
                         beta.fl.per.head, beta.seed.mass, beta.infest,
                         sample.size)
{
  n.obs <- length(species)
  pink <- numeric(n.obs)
  for (i in 1:n.obs) {
    phi <- beta.0[species[i]] +
           beta.long[repens[i]]*long[i] +
           beta.map[repens[i]]*map[i] +
           beta.elev[repens[i]]*elev[i] +
           beta.fecundity[repens[i]]*fecundity[i] +
           beta.fl.per.head[repens[i]]*fl.per.head[i] +
           beta.seed.mass[repens[i]]*seed.mass[i] +
           beta.infest[repens[i]]*infest[i]
    pi <- inv.logit(phi)
    pink[i] <- rbinom(1, 20, pi)
  }
  list(pink=pink, sample.size=sample.size)
}

get.pink.cat <- function(pink, sample.size) {
  category <- numeric(length(pink))
  for (i in 1:length(pink)) {
    if (pink[i] == 0) {
      category[i] <- 1
    } else if (pink[i] < 0.10*sample.size) {
      category[i] <- 2
    } else if (pink[i] < 0.90*sample.size) {
      category[i] <- 3
    } else if (pink[i] < sample.size) {
      category[i] <- 4
    } else {
      category[i] <- 5
    }
  }
  category
}

rescale <- function(x) {
  for (i in 1:ncol(x)) {
    if (is.numeric(x[,i]) &&
        !(colnames(x[i]) %in% c("poly", "pink", "repens")))
    {
      x[,i] <- scale(x[,i])
    }
  }
  x
}

printInterval <- function(int, name) {
  if ((int[1,1] > 0) || (int[1,2] < 0)) {
    out <- sprintf("%33s  %8.3f %8.3f", name, int[1,1], int[1,2])
    cat(out)
    cat("*")
    cat("\n")
  }
}

SYMinterval <- function(x, prob) {
  lo <- (1.0 - prob)/2.0
  hi <- 1.0 - lo
  sym <- matrix(nrow=1, ncol=2)
  sym[1,] <- quantile(x, c(lo, hi))
  sym
}

intervalvector <- function(x, name, prob, type) {
  var <- x[[name]]
  n.cat <- dim(var)[2]
  for (i in 1:n.cat) {
    new.name <- paste(name, "[", i, "]", sep="")
    if (type == "symmetric") {
      sym <- SYMinterval(mcmcUpgrade(as.mcmc(var[,i])), prob=prob)
      printInterval(sym, new.name)
    } else if (type == "HPD") {
      hpd <- HPDinterval(mcmcUpgrade(as.mcmc(var[,i])), prob=prob)
      printInterval(hpd, new.name)
    }
  }
}

intervals <- function(x, prob, type) {
  x <- x$BUGSoutput$sims.list
  out <- sprintf("%2.0f%% %s intervals\n", prob*100, type)
  cat(out)
  out <- sprintf("%33s  %8s %8s\n", "Coefficient", "lo", "hi")
  cat(out)
  out <- sprintf("%33s  %8s %8s\n",
                 "---------------------------------",
                 "-------",
                 "-------")
  cat(out)
  for (name in names(x)) {
    if (name %in% c("alpha.poly.0", "alpha.pink.0",
                    "eps.fecundity.species", "eps.fl.per.head.species",
                    "eps.infest.species", "eps.seed.mass.species"))
    {
      intervalvector(x, name, prob, type)
    } else if (name %in% c("beta.elev", "beta.fecundity",
                           "beta.fl.per.head", "beta.infest",
                           "beta.long", "beta.map",
                           "beta.seed.mass", "beta.0"))
    {
      var <- x[[name]]
      n.cat <- dim(var)[2:3]
      for (i in 1:n.cat[1]) {
        for (j in 1:n.cat[2]) {
           new.name <- paste(name, "[", i, ",", j, "]", sep="")
           if (type == "symmetric") {
             sym <- SYMinterval(mcmcUpgrade(as.mcmc(var[,i,j])), prob=prob)
             printInterval(sym, new.name)
           } else if (type == "HPD") {
             hpd <- HPDinterval(mcmcUpgrade(as.mcmc(var[,i,j])), prob=prob)
             printInterval(hpd, new.name)
           }
        }
      }
    } else if (name %in% c("beta.fecundity.elev", "beta.fecundity.long",
                           "beta.fecundity.map", "beta.fl.per.head.elev",
                           "beta.fl.per.head.long", "beta.fl.per.head.map",
                           "beta.infest.elev", "beta.infest.fecundity",
                           "beta.infest.fl.per.head", "beta.infest.long",
                           "beta.infest.map", "beta.infest.seed.mass",
                           "beta.seed.mass.elev", "beta.seed.mass.long",
                           "beta.seed.mass.map"))
    {
      var <- x[[name]]
      for (i in 1:2) {
        new.name <- paste(name, "[", i, "]", sep="")
        if (type == "symmetric") {
          sym <- SYMinterval(mcmcUpgrade(as.mcmc(var[,i])), prob=prob)
          printInterval(sym, new.name)
        } else if (type == "HPD") {
          hpd <- HPDinterval(mcmcUpgrade(as.mcmc(var[,i])), prob=prob)
          printInterval(hpd, new.name)
        }
      }
    } else {
      if (name != "deviance") {
        if (type == "symmetric") {
          sym <- SYMinterval(mcmcUpgrade(as.mcmc(x[[name]])), prob=prob)
          printInterval(sym, name)
        } else if (type == "HPD") {
          hpd <- HPDinterval(mcmcUpgrade(as.mcmc(x[[name]])), prob=prob)
          printInterval(hpd, name)
        }
      }
    }
  }
}

logistic.intervals <- function(x, prob, type) {
  x <- x$BUGSoutput$sims.list
  out <- sprintf("%2.0f%% %s intervals\n", prob*100, type)
  cat(out)
  out <- sprintf("%33s  %8s %8s\n", "Coefficient", "lo", "hi")
  cat(out)
  out <- sprintf("%33s  %8s %8s\n",
                 "---------------------------------",
                 "-------",
                 "-------")
  cat(out)
  for (name in names(x)) {
    var <- x[[name]]
    if (name != "deviance") {
      for (i in 1:2) {
        new.name <- paste(name, "[", i, "]", sep="")
        if (type == "symmetric") {
          sym <- SYMinterval(mcmcUpgrade(as.mcmc(var[,i])), prob=prob)
          printInterval(sym, new.name)
        } else if (type == "HPD") {
          hpd <- HPDinterval(mcmcUpgrade(as.mcmc(var[,i])), prob=prob)
          printInterval(hpd, new.name)
        }
      }
    }
  }
}

compare.coefficients <- function(x, n.sample, prob=0.95, digits=3) {
  cov <- c("elev", "fecundity", "fl.per.head", "infest", "long",
           "map", "seed.mass")
  for (trait in cov) {
    name <- paste("beta.", trait, sep="")
    y <- x$BUGSoutput$sims.list[[name]]
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


## read in data
##
color.csv <- read.csv("color.csv", header=TRUE, na.strings=".")
color.csv <- subset(color.csv, species!="subvestita")
color.csv <- drop.levels(color.csv)

species <- color.csv$species
site <- color.csv$site
fecundity <- color.csv$seeds_per_SH
fl.per.head <- color.csv$flowers_per_head
seed.mass <- color.csv$seed_mass
infest <- color.csv$prop_infested_per_plant
long <- color.csv$long
map <- color.csv$gmap
elev <- color.csv$altitud
repens <- as.numeric(color.csv$species=="repens") + 1
poly <- diversity.to.poly(color.csv$diversity)
pink <- get.pink.scale(color.csv$ratio, n.cat)

color <- data.frame(species=species,
                    repens=repens,
                    site=site,
                    fecundity=fecundity,
                    fl.per.head=fl.per.head,
                    seed.mass=seed.mass,
                    infest=infest,
                    long=long,
                    map=map,
                    elev=elev,
                    poly=poly,
                    pink=pink)
rm(species, repens, site, fecundity, fl.per.head, seed.mass,
   infest, long, map, elev, poly, pink)

## exclude sites lacking fecundity, fl.per.head, seed.mass, or infest
##
traits <- c("fecundity", "fl.per.head", "seed.mass", "infest")
color <- color[!is.na(apply(color[,traits], 1, sum)),]
## exclude sites lacking poly
##
color <- color[!is.na(color$poly),]
color$site <- color$site[,drop=TRUE]
color$species <- color$species[,drop=TRUE]
## exclude sites lacking pink
##
color <- color[!is.na(color$pink),]
color$site <- color$site[,drop=TRUE]
color$species <- color$species[,drop=TRUE]

## recover categorical variables
##
species <- as.numeric(color$species)
repens <- color$repens
poly <- color$poly
pink <- color$pink

## save unscaled variates
##
unscaled <- data.frame(species=color$species,
                       elev=color$elev,
                       fecundity=color$fecundity,
                       fl.per.head=color$fl.per.head,
                       infest=color$infest,
                       long=color$long,
                       seed.mass=color$seed.mass,
                       map=color$map)

## rescale covariates
##
color <- rescale(color)

## get rid of scaled attributes from rescaling
##
fecundity <- as.numeric(color$fecundity)
fl.per.head <- as.numeric(color$fl.per.head)
seed.mass <- as.numeric(color$seed.mass)
infest <- as.numeric(color$infest)
long <- as.numeric(color$long)
map <- as.numeric(color$map)
elev <- as.numeric(color$elev)

species.table <- unique(data.frame(species=color$species,
                                   number=as.numeric(color$species)))
n.obs <- nrow(color)
n.species <- nrow(species.table)

tau.beta <- 0.1
tau.species <- 0.1

if (betas.at.zero) {
  ## checking false positive rate
  ##
  beta.0 <- numeric(n.species)
  for (i in 1:n.species) {
    beta.0[i] <- 0.0
  }
  beta.long <- numeric(2)
  beta.map <- numeric(2)
  beta.elev <- numeric(2)
  beta.fecundity <- numeric(2)
  beta.fl.per.head <- numeric(2)
  beta.seed.mass <- numeric(2)
  beta.infest <- numeric(2)
  for (i in 1:2) {
    beta.long[i] <- 0.0
    beta.map[i] <- 0.0
    beta.elev[i] <- 0.0
    beta.fecundity[i] <- 0.0
    beta.fl.per.head[i] <- 0.0
    beta.seed.mass[i] <- 0.0
    beta.infest[i] <- 0.0
  }
} else {
  ## checking power
  ##
  ## replace observed pink with predicted
  ##
  load("results-full.Rsave")
  means <- fit$BUGSoutput$mean
  ## minus because log odds for adjacent categories are
  ##
  ## exp(-L.pink)
  ##
  beta.0 <- -colMeans(means$beta.pink.0)
  beta.long <- -colMeans(means$beta.pink.long)
  beta.map <- -colMeans(means$beta.pink.map)
  beta.elev <- -colMeans(means$beta.pink.elev)
  beta.fecundity <- -colMeans(means$beta.pink.fecundity)
  beta.fl.per.head <- -colMeans(means$beta.pink.fl.per.head)
  beta.seed.mass <- -colMeans(means$beta.pink.seed.mass)
  beta.infest <- -colMeans(means$beta.pink.infest)
}

## split results between terminal and results file
##
results <- sub("power-analysis", "results", model.file)
sink(file=results, split=TRUE, append=TRUE)
for (i in 1:n.species) {
  cat("beta.0[", i, "]:           ", round(beta.0[i], 3), "\n", sep="")
}
for (i in 1:2) {
  cat("beta.long[", i, "]:        ", round(beta.long[i], 3), "\n", sep="")
  cat("beta.map[", i, "]:         ", round(beta.map[i], 3), "\n", sep="")
  cat("beta.elev[", i, "]:        ", round(beta.elev[i], 3), "\n", sep="")
  cat("beta.fecundity[", i, "]:   ", round(beta.fecundity[i], 3), "\n", sep="")
  cat("beta.fl.per.head[", i, "]: ", round(beta.fl.per.head[i], 3), "\n", sep="")
  cat("beta.seed.mass[", i, "]:   ", round(beta.seed.mass[i], 3), "\n", sep="")
  cat("beta.infest[", i, "]:      ", round(beta.infest[i], 3), "\n", sep="")
}
cat("\n")
sink()

for (i in 1:n.simulations) {
  freq <- predict.pink(species, repens, long, map, elev, fecundity, fl.per.head,
                       seed.mass, infest,
                       beta.0, beta.long, beta.map, beta.elev, beta.fecundity,
                       beta.fl.per.head, beta.seed.mass, beta.infest,
                       sample.size)
  pink <- get.pink.cat(freq$pink, freq$sample.size)

  jags.data <- c("fecundity",
                 "fl.per.head",
                 "seed.mass",
                 "infest",
                 "long",
                 "map",
                 "elev",
                 "pink",
                 "species",
                 "repens",
                 "n.obs",
                 "n.species",
                 "n.cat",
                 "tau.beta",
                 "tau.species")
  jags.parameters <- c("beta.fecundity",
                       "beta.fl.per.head",
                       "beta.seed.mass",
                       "beta.infest",
                       "beta.long",
                       "beta.map",
                       "beta.elev")

  fit <- jags(jags.data,
              inits=NULL,
              parameters=jags.parameters,
              model.file=model.file,
              n.chains=n.chains,
              n.burnin=n.burnin,
              n.iter=n.iter,
              n.thin=n.thin)

  ## split results between terminal and results file
  ##
  sink(file=results, split=TRUE, append=TRUE)

  model <- sub(".txt", "", model.file)
  cat(model, "\n")

  intervals(fit, prob=0.95, "HPD")

  ## close output file
  ##
  sink()

  pink <- freq$pink
  n.sample <- freq$sample.size
  jags.data <- c("fecundity",
                 "fl.per.head",
                 "seed.mass",
                 "infest",
                 "long",
                 "map",
                 "elev",
                 "pink",
                 "species",
                 "repens",
                 "n.obs",
                 "n.species",
                 "n.sample",
                 "tau.beta",
                 "tau.species")
  fit <- jags(jags.data,
              inits=NULL,
              parameters=jags.parameters,
              model.file=model.file.logistic,
              n.chains=n.chains,
              n.burnin=n.burnin,
              n.iter=n.iter,
              n.thin=n.thin)

  ## split results between terminal and results file
  ##
  sink(file=results, split=TRUE, append=TRUE)

  model <- sub(".txt", "", model.file.logistic)
  cat(model, "\n")

  logistic.intervals(fit, prob=0.95, "HPD")

  ## close output file
  ##
  sink()
}

elapsed <- proc.time() - ptm
print(elapsed)

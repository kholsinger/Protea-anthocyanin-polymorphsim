require(R2jags)

rm(list=ls())
options("width"=120)

debug <- FALSE
exclude.jkr.and.gp2 <- TRUE
jkr.and.gp2.in.data <- FALSE

## to seed RNG
##
runif(1)

n.poly.cat <- 3
n.pink.cat <- 5

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
  out <- sprintf("%33s  %8.3f %8.3f", name, int[1,1], int[1,2])
  cat(out)
  if ((int[1,1] > 0) || (int[1,2] < 0)) {
    cat("*")
  }
  cat("\n")
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
    } else if (name %in% c("beta.poly.elev", "beta.poly.fecundity",
                           "beta.poly.fl.per.head", "beta.poly.infest",
                           "beta.poly.long", "beta.poly.map",
                           "beta.poly.seed.mass", "beta.poly.0",
                           "beta.pink.elev", "beta.pink.fecundity",
                           "beta.pink.fl.per.head", "beta.pink.infest",
                           "beta.pink.long", "beta.pink.map",
                           "beta.pink.seed.mass", "beta.pink.0",
                           "Sigma"))
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


## read in data
##
color.csv <- read.csv("color.csv", header=TRUE, na.strings=".")
color.csv <- subset(color.csv, species!="subvestita")
if (exclude.jkr.and.gp2) {
  color.csv <- subset(color.csv, site!="JKR")
  color.csv <- subset(color.csv, site!="GP2")
}
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
poly <- diversity.to.poly(color.csv$diversity)
pink <- get.pink.scale(color.csv$ratio, n.pink.cat)
repens <- as.numeric(color.csv$species=="repens") + 1

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

y <- as.matrix(data.frame(fecundity,
                          fl.per.head,
                          seed.mass))

Omega <- diag(x=1.0, nrow=3, ncol=3)
wish.nu <- nrow(Omega) + 2

jags.data <- c("fecundity",
               "fl.per.head",
               "seed.mass",
               "infest",
               "long",
               "map",
               "elev",
               "poly",
               "pink",
               "species",
               "repens",
               "n.obs",
               "n.species",
               "n.poly.cat",
               "n.pink.cat",
               "tau.beta",
               "tau.species",
               "y",
               "Omega",
               "wish.nu")
jags.parameters <- c("alpha.poly.fecundity",
                     "alpha.poly.fl.per.head",
                     "alpha.poly.seed.mass",
                     "alpha.poly.infest",
                     "alpha.poly.long",
                     "alpha.poly.map",
                     "alpha.poly.elev",
                     "alpha.poly.0",
                     "alpha.pink.fecundity",
                     "alpha.pink.fl.per.head",
                     "alpha.pink.seed.mass",
                     "alpha.pink.infest",
                     "alpha.pink.long",
                     "alpha.pink.map",
                     "alpha.pink.elev",
                     "alpha.pink.0",
                     "beta.poly.fecundity",
                     "beta.poly.fl.per.head",
                     "beta.poly.seed.mass",
                     "beta.poly.infest",
                     "beta.poly.long",
                     "beta.poly.map",
                     "beta.poly.elev",
                     "beta.poly.0",
                     "beta.pink.fecundity",
                     "beta.pink.fl.per.head",
                     "beta.pink.seed.mass",
                     "beta.pink.infest",
                     "beta.pink.long",
                     "beta.pink.map",
                     "beta.pink.elev",
                     "beta.pink.0",
                     "beta.infest.fecundity",
                     "beta.infest.fl.per.head",
                     "beta.infest.seed.mass",
                     "beta.infest.long",
                     "beta.infest.map",
                     "beta.infest.elev",
                     "eps.infest.species",
                     "beta.fecundity.long",
                     "beta.fecundity.map",
                     "beta.fecundity.elev",
                     "eps.fecundity.species",
                     "tau.fecundity",
                     "beta.fl.per.head.long",
                     "beta.fl.per.head.map",
                     "beta.fl.per.head.elev",
                     "eps.fl.per.head.species",
                     "tau.fl.per.head",
                     "beta.seed.mass.long",
                     "beta.seed.mass.map",
                     "beta.seed.mass.elev",
                     "eps.seed.mass.species",
                     "tau.seed.mass",
                     "Sigma")

csv <- data.frame(model=NULL,
                  Dbar=NULL,
                  Dhat=NULL,
                  pD=NULL,
                  DIC=NULL)
for (model.file in c("analysis-full.txt")) {
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
  results <- sub("analysis", "results", model.file)
  sink(file=results, split=TRUE)

  model <- sub("analysis-", "", model.file)
  model <- sub(".txt", "", model.file)
  cat(model, "\n")
  if (jkr.and.gp2.in.data) {
    if (exclude.jkr.and.gp2) {
      cat("GP2 and JKR excluded from analysis\n")
    } else {
      cat("GP2 and JKR included in analysis\n")
    }
  }
  DIC <- fit$BUGSoutput$DIC
  pD <- fit$BUGSoutput$pD
  Dbar <- DIC - pD
  Dhat <- Dbar - pD
  cat("  Dbar: ", Dbar, "\n", sep="")
  cat("  Dhat: ", Dhat, "\n", sep="")
  cat("    pD: ", pD, "\n", sep="")
  cat("   DIC: ", DIC, "\n", sep="")
  tmp <- data.frame(model=model,
                    Dbar=Dbar,
                    Dhat=Dhat,
                    pD=pD,
                    DIC=DIC)
  csv <- rbind(csv, tmp)

  cat("\n\n\n", model.file, "results...")
  print(fit, digits.summary=3)
  cat("\n\n\n")
  intervals(fit, prob=0.95, "symmetric")
  intervals(fit, prob=0.95, "HPD")
  intervals(fit, prob=0.90, "symmetric")
  intervals(fit, prob=0.90, "HPD")

  cat("\n\n\n")
  source("compare-betas.R")
  compare(fit, prob=0.95)

  rsave <- paste(sub(".txt", "", model.file), ".Rsave", sep="")
  rsave <- sub("analysis", "results", rsave)
  save(fit, species.table, color, unscaled, file=rsave)

  ## close output file
  ##
  sink()
}

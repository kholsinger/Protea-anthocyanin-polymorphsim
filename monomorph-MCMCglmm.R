library(MCMCglmm)
library(car)
library(ape)

rm(list=ls())

na.strings <- c("NA", ".")
test <- FALSE
valente <- TRUE

n.reps <- 5

if (test) {
  n.sample <- 10000
  n.burnin <- 5000
  n.thin <- 1
  n.reps <- 2
} else {
  n.sample <- 10000000
  n.burnin <- 1000000
  n.thin <- 2000
}
n.iter <- n.burnin + n.sample

columns <- c("elevation",
             "long",
             "rainfall",
             "map",
             "tmax",
             "tmin",
             "mat",
             "rainy",
             "sp.range")


## Drop unused factor levels from all factors in a data.frame
## Author: Kevin Wright.  Idea by Brian Ripley.
##
drop.levels <- function(dat) {
  dat[] <- lapply(dat, function(x) x[,drop=TRUE])
  return(dat)
}

HPDintervals <- function(x, prob) {
  out <- sprintf("\n\n\n\n%2.0f%% HPD intervals\n", prob*100)
  cat(out)
  out <- sprintf("%33s  %8s %8s\n", "Coefficient", "lo", "hi")
  cat(out)
  out <- sprintf("%33s  %8s %8s\n",
                 "---------------------------------",
                 "-------",
                 "-------")
  cat(out)
  names <- colnames(x[[1]])
  for (name in names) {
    hpd <- HPDinterval(as.mcmc(unlist(x[,name])), prob=prob)
    out <- sprintf("%33s  %8.3f %8.3f\n", name, hpd[1,1], hpd[1,2])
    cat(out)
  }
}

standardize <- function(x) {
  if (is.numeric(x)) {
    y <- (x - mean(x, na.rm=TRUE))/sd(x, na.rm=TRUE)
  } else {
    y <- x
  }
  y
}

## read in data for traits, climate, merge them.
sp.color <- read.csv("pink_white_monomorphic_species.csv",
                    na.strings=".",
                    header=TRUE)
sp.color$animal <- sp.color$species
species <- sp.color$species
animal <- sp.color$animal
color <- sp.color$mono_color
elevation <- sp.color$elevation
rainfall <- sp.color$rainfall
map <- sp.color$map
tmax <- sp.color$tmin
tmin <- sp.color$tmin
long <- sp.color$long
mat <- sp.color$mat
rainy <- sp.color$rainy
sp.range <- sp.color$range

## put variables into new data frame
tmp <- data.frame(species=species,
                  animal=animal,
                  color=color,
                  elevation=elevation,
                  rainfall=rainfall,
                  map=map,
                  tmax=tmax,
                  tmin=tmin,
                  long=long,
                  mat=mat,
                  rainy=rainy,
                  sp.range=sp.range)
## remove polymorphic species
##
tmp <- subset(tmp, color != "poly")
## get rid of incomplete cases
ok <- complete.cases(tmp)
tmp <- tmp[ok,]
tmp <- drop.levels(tmp)

## 3. pull the vectors back out and standardize
raw.data <- tmp
species <- tmp$species
animal <- tmp$animal
color <- tmp$color
elevation <- standardize(tmp$elevation)
rainfall <- standardize(tmp$rainfall)
map <- standardize(tmp$map)
tmax <- standardize(tmp$tmax)
tmin <- standardize(tmp$tmin)
long <- standardize(tmp$long)
mat <- standardize(tmp$mat)
rainy <- standardize(tmp$rainy)
sp.range <- standardize(tmp$sp.range)

tmp <- data.frame(species=species,
                  animal=animal,
                  color=color,
                  elevation=elevation,
                  rainfall=rainfall,
                  map=map,
                  tmax=tmax,
                  tmin=tmin,
                  long=long,
                  mat=mat,
                  rainy=rainy,
                  sp.range=sp.range)
fitted.data <- tmp

## cladogram
##
clades <- read.tree("valentetree.tre")
## drop taxa that are not included in data set
##
phylo <- drop.tip(clades, setdiff(clades$tip.label, tmp$animal))
phylo <- makeNodeLabel(phylo)

n.fixed <- 1 + 3 + 1

prior <- list(R=list(V=1, nu=1, fix=1),
              G=list(G1=list(V=1, nu=3)))
results.mcmc <- list()
results.VCV <- list()
results.store <- list()
results.DIC <- numeric(n.reps)
results.deviance <- matrix(nrow=n.reps, ncol=n.sample/n.thin)
for (i in 1:n.reps) {
  cat("rep", i, "\n")
  flush.console()
  results <- MCMCglmm(color ~ elevation,
                      random = ~ animal,
                      data=tmp,
                      prior=prior,
                      pedigree=phylo,
                      nitt=n.iter,
                      burnin=n.burnin,
                      thin=n.thin,
                      family="categorical",
                      pr=TRUE,
                      verbose=FALSE,
                      slice=TRUE)
  results.mcmc[[i]] <- as.mcmc(results$Sol)
  results.VCV[[i]] <- as.mcmc(results$VCV)
  results.store[[i]] <- results
  results.DIC[i] <- results$DIC
  results.deviance[i,] <- results$Deviance
}
results.mcmc <- as.mcmc.list(results.mcmc)
results.VCV <- as.mcmc.list(results.VCV)

print(summary(results.mcmc, quantiles=c(0.025, 0.05, 0.5, 0.95, 0.975)))
HPDintervals(results.mcmc, prob=0.95)
HPDintervals(results.mcmc, prob=0.90)
if (n.reps > 1) {
  cat("\n\n\n\n")
  print(gelman.diag(results.mcmc))
}

save(results.mcmc, raw.data, fitted.data, file="poly-results.Rsave")

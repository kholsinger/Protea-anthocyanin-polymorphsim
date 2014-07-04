library(MCMCglmm)
library(car)
library(ape)

rm(list=ls())

na.strings <- c("NA", ".")
test <- FALSE
valente <- TRUE

drop.venusta <- FALSE

n.reps <- 5

if (test) {
  n.sample <- 10000
  n.burnin <- 5000
  n.thin <- 1
} else {
  n.sample <- 100000000
  n.burnin <- 10000000
  n.thin <- 100000
}
n.iter <- n.burnin + n.sample

columns <- c("ele",
             "ele.range",
             "long",
             "long.range",
             "map",
             "map.range",
             "area")


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
means <- read.csv("Carlson_colorDataAll.csv",
                    na.strings=".",
                    header=TRUE)
means <- subset(means, means$include_1 == 1)
if (valente) {
  means$animal <- means$Valente_name
} else {
  means$animal <- means$Schnitzler_name
}

## pull out desired variables
if (valente) {
  species <-means$Valente_name
} else {
  species <- means$Schnitzler_name
}
poly <- means$polymorph_yes
ele <- means$elevation
ele.range <- means$MAX_ELE - means$MIN_ELE
map <- means$MAP
map.range <- means$MAP_MAX - means$MAP_MIN
area <- means$range_area_albers
long <- means$long
long.range <- means$Max_long - means$min_long
animal <-means$animal

## put variables into new data frame
tmp <- data.frame(species=species,
                  poly=poly,
                  ele=ele,
                  ele.range=ele.range,
                  map=map,
                  map.range=map.range,
                  area=area,
                  long=long,
                  long.range=long.range,
                  animal=animal)
if (drop.venusta) {
  if (valente) {
    tmp <- subset(tmp, species != "Prvenu")
  } else {
    tmp <- subset(tmp, species != "Protea_venusta")
  }
}
## get rid of incomplete cases
ok <- complete.cases(tmp)
tmp <- tmp[ok,]
tmp <- drop.levels(tmp)

## 3. pull the vectors back out and standardize
raw.data <- tmp
species <- tmp$species
poly <- tmp$poly
animal <-tmp$animal
ele <- standardize(tmp$ele)
ele.range <- standardize(tmp$ele.range)
map <- standardize(tmp$map)
map.range <- standardize(tmp$map.range)
area <- standardize(tmp$area)
long <- standardize(tmp$long)
long.range <- standardize(tmp$long.range)

tmp <- data.frame(species=species,
                  poly=poly,
                  ele=ele,
                  ele.range=ele.range,
                  map=map,
                  map.range=map.range,
                  area=area,
                  long=long,
                  long.range=long.range,
                  animal=animal)
fitted.data <- tmp

## cladogram
##
if (valente) {
  clades <- read.tree("valentetree.tre")
  clades$tip.label <- tolower(clades$tip.label)
  substr(clades$tip.label, 1, 1) <- "P"
} else {
  clades <- read.nexus("schnitzler.nex")
}
## drop taxa that are not included in means data set
##
phylo <- drop.tip(clades, setdiff(clades$tip.label, tmp$animal))
phylo <- makeNodeLabel(phylo)

n.fixed <- length(columns) + 3 + 1

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
  results <- MCMCglmm(poly ~ ele + ele.range + long + long.range +
                             map + map.range + area,
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

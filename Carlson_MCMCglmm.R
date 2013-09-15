library(MCMCglmm)
library(car)
library(ape)

rm(list=ls())


phylo <- read.nexus("schnitzler.nex")
phylonew<-drop.tip(phylo, c("Faurea_galpinii",	"Faurea_macnaughtonii",	"Faurea_rochetiana",	"Faurea_rubriflora",	"Faurea_saligna",	"Leucadendron_gandogeri",	"Protea_acaulos",	"Protea_acuminata",	"Protea_amplexicaulis",	"Protea_angolensis",	"Protea_angustata",	"Protea_aspera",	"Protea_caespitosa",	"Protea_caffra",	"Protea_canaliculata",	"Protea_comptonii",	"Protea_convexa",	"Protea_cordata",	"Protea_coronata",	"Protea_curvata",	"Protea_decurrens",	"Protea_dracomontana",	"Protea_effusa",	"Protea_enervis",	"Protea_ericifolia",	"Protea_foliosa",	"Protea_gaguedi",	"Protea_glabra",	"Protea_heckmanniana",	"Protea_holosericea",	"Protea_humiflora",	"Protea_inopina",	"Protea_laetans",	"Protea_laevis",	"Protea_lanceolata",	"Protea_lorea",	"Protea_mucronifolia",	"Protea_mundiieast",	"Protea_namaquana",	"Protea_nana",	"Protea_nitida",	"Protea_nubigena",	"Protea_odorata",	"Protea_parvula",	"Protea_pendula",	"Protea_pityphylla",	"Protea_pruinosa",	"Protea_pudens",	"Protea_recondita",	"Protea_restionifolia",	"Protea_revoluta",	"Protea_rubropilosa",	"Protea_rupicola",	"Protea_scabra",	"Protea_scabriuscula",	"Protea_scolymocephala",	"Protea_scorzonerifolia",	"Protea_simplex",	"Protea_stokoei",	"Protea_subulifolia",	"Protea_sulphurea",	"Protea_welwitschii",	"Protea_wentzeliana",	"Protea_witzenbergiana",	"Protea_woeskaensis",	"Serruria_barbigera",	"Serruria_scoparia",	"Serruria_villosa"))

phylo2 <- read.nexus("valentetree.nex")
phylo2<-drop.tip(phylo2, c("Pracau",	"Pracum",	"Prampl",	"Prango",	"Prangu",		"Praspe",			"Prcab",	"Prcaes",	"Prcaff",	"Prcana",	"Prcomp",	"Prconv",	"Prcord",	"Prcoro",			"Prcurv",		"Prdecu",		"Prdrac",	"Preffu",	"Prener",		"Prfoli",	"Prgagu",	"Prglab",		"Prheck",	"Prholo",	"Prhumi",	"Prinop",			"Prlaet",	"Prlaev",	"Prlanc",				"Prlore",				"Prmucr",		"Prnama",	"Prnana",		"Prniti",	"Prnubi",		"Prodor",	"Prparv",	"Prpend",		"Prpity",	"Prprui",	"Prpude",		"Prreco",		"Prrest",	"Prrevo",		"Prrubr",	"Prrupi",	"Prscbr",		"Prscol",	"Prscor",	"Prsimp",		"Prstok",	"Prsubu",		"Prsulp",					"Prwelw",	"Prwent",	"Prwitz"))



na.strings <- c("NA", ".")
plot <- FALSE
test <- FALSE
timtree <- TRUE
nprs_trim <- FALSE


n.reps <- 5

if (test) {
  n.sample <- 5000
  n.burnin <- 10000
  n.thin <- 1
} else {
  n.sample <- 10000000
  n.burnin <- 1000000
  n.thin <- 10000
}
n.iter <- n.burnin + n.sample

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












## MCMCglmm ##

standardize <- function(x) {
  if (is.numeric(x)) {
    y <- (x - mean(x, na.rm=TRUE))/sd(x, na.rm=TRUE)
  } else {
    y <- x
  }
  y
}

drop.levels <- function(dat) {
  dat[] <- lapply(dat, function(x) x[,drop=TRUE])
  return(dat)
}


#polymorph_yes	elevation	rainfall_dry_season	MAP	long	
#Schnitzler_name	Valente_name

## read in data for traits, climate, merge them.
means <- read.csv("Carlson_colorDataAll.csv",
                    na.strings=".",
                    header=TRUE)

means$animal <- means$Schnitzler_name
means$animal

## pull out desired variables
species <-means$Schnitzler_name
poly <- means$polymorph_yes
alt <- means$elevation
map <- means$MAP
long <- means$long
dry <- means$rainfall_dry_season
animal <-means$animal

## put variables into new data frame
tmp <- data.frame(species=species,
                  poly=poly,
                  alt=alt,
                  map=map,
                  long=long,
                  dry=dry,
			animal=animal)
## get rid of incomplete cases
ok <- complete.cases(tmp)
tmp <- tmp[ok,]
tmp <- drop.levels(tmp)


columns <- c("ele",
	"long",	
	"map",	
	"dry")

if (plot) {
  plot(means[,columns])
}




## 3. pull the vectors back out and standardize
species <- tmp$species
poly <- tmp$poly
animal <-tmp$animal
alt <- standardize(tmp$alt)
map <- standardize(tmp$map)
long <- standardize(tmp$long)
dry <- standardize(tmp$dry)

tmp <- data.frame(species=species,
                  poly=poly,
                  alt=alt,
                  map=map,
                  long=long,
                  dry=dry,
			animal=animal)
phylo <- read.nexus("schnitzler.nex")


## cladogram
##
if (timtree) {
  clades <- read.nexus("schnitzler.nex")
  clades <- chronopl(clades, 0.1)
} else if (nprs_trim) {
  clades <- read.nexus("schnitzler.nex")
}
clades <- makeNodeLabel(clades)

## drop taxa that are not included in means data set
##
phylo <- drop.tip(phylo, setdiff(phylo$tip.label, means$animal))


## drop taxa not included in cladogram
##
#x <- means[match(phylo$tip.label, means$animal),]
#x$animal <- x$animal[,drop=TRUE]
#not.included <- setdiff(levels(means$animal), levels(x$animal))
#cat("Taxa not included in phylogeny...\n")
#for (i in not.included) {
#  cat("  ", i, "\n", sep="")
#}
#means <- x
#rm(x)


n.fixed <- length(columns) + 3 + 1

prior <- list(R=list(V=1, nu=1, fix=1),
              G=list(G1=list(V=1, nu=3)))
results.mcmc <- list()
results.VCV <- list()
results.store <- list()
results.DIC <- numeric(n.reps)
results.deviance <- numeric(n.reps)
for (i in 1:n.reps) {
  cat("rep", i, "\n")
  flush.console()
  results <- MCMCglmm(poly ~ alt
						+ long,
                      random = ~animal,
                      data=tmp,
                      prior=prior,
                      pedigree=phylo,
                      nitt=n.iter,
                      burnin=n.burnin,
                      thin=n.thin,
                      family="categorical",
                      verbose=FALSE,
                      slice=TRUE)
  results.mcmc[[i]] <- as.mcmc(results$Sol)
  results.VCV[[i]] <- as.mcmc(results$VCV)
  results.store[[i]] <- results
  results.DIC[i] <- results$DIC
  results.deviance[i] <- results$Deviance
}
results.mcmc <- as.mcmc.list(results.mcmc)
results.VCV <- as.mcmc.list(results.VCV)

print(summary(results.mcmc, quantiles=c(0.025, 0.05, 0.5, 0.95, 0.975)))
HPDintervals(results.mcmc, prob=0.95)
HPDintervals(results.mcmc, prob=0.90)
cat("\n\n\n\n")
print(gelman.diag(results.mcmc))



#### Actually using mcmcglmm
#G <- vcv(phylo, corr=TRUE)
#Ginv <- solve(G)

## convert Ginv to a dgCMatrix to use the MCMCglmm function
#Ginv<-as(Ginv, Class="dgCMatrix")
#prior = list(R = list(V = diag(6), nu = 0, fix = 4), G = list(G1 = list(V = diag(6), nu = 1)))
#model<-MCMCglmm(cbind(lma, area, lwr, fwc)~map+mat+mtmax+mtmin+ratio+cdd, random=~species, family=c("gaussian","gaussian","gaussian","gaussian"), nodes="TIPS", ginverse=list(species=Ginv),
#              rcov=~idh(trait):units,  data=tmp, pr=TRUE, pl=TRUE, verbose=TRUE)
#plot(model2$VCV)
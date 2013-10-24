require(ape)

clades <- read.nexus("valentetree.nex")
clades <- chronos(clades)
write.tree(clades, file="valentetree.tre")

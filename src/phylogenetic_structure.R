# Community Ecology and Phylogenetics

# I'm following this tutorial: http://ib.berkeley.edu/courses/ib200/2016/labs/14/lab14.R
# Also Picante manual: https://cran.r-project.org/web/packages/picante/vignettes/picante-intro.pdf

###############################################
# Analyzing NRI and NTI for Epidemic analysis #
###############################################

library(picante)
library(ggtree)
library(treeio)
library(phangorn)

# Load tree
#chrTree <- read.tree('/Users/arafat/GDrive/Sachs/Chapter3_Epidemic_Genotype/Sequences/Bradyrhizobium/Dataset/Complete/Dataset/tanglegram/chr.meta.unique.treefile')
#chrTree <- midpoint(chrTree)
#symTree <- read.tree('/Users/arafat/GDrive/Sachs/Chapter3_Epidemic_Genotype/Sequences/Bradyrhizobium/Dataset/Complete/Dataset/tanglegram/sym.meta.unique.treefile')
#symTree <- midpoint(symTree)

chrTree <- read.iqtree('/Users/arafat/GDrive/Sachs/Chapter3_Epidemic_Genotype/Sequences/Bradyrhizobium/Dataset/Complete/Dataset/chr/concatenate/aln/chr.concatenated.partitions.rxml.treefile')
chrTree <- midpoint(chrTree@phylo)

symTree <- read.iqtree('/Users/arafat/GDrive/Sachs/Chapter3_Epidemic_Genotype/Sequences/Bradyrhizobium/Dataset/Complete/Dataset/sym/concatenate/aln/sym.concatenated.partitions.rxml.treefile')
symTree <- midpoint(symTree@phylo)


par(mfrow = c(1, 2))
plot(chrTree, show.tip.label = FALSE, main='CHR tree')
plot(symTree, show.tip.label = FALSE, main='SYM tree')


# Prepare community
meta <- read.csv('/Users/arafat/GDrive/Sachs/Chapter3_Epidemic_Genotype/Sequences/Bradyrhizobium/Dataset/Complete/Dataset/meta_all.csv')
dim(meta)


# Update both tree Phylo object's tip name from isolate to haplotype name
meta.sub <- meta[which(meta$Strain %in% chrTree$tip.label), c('CHR_haplotype', 'Strain')]
meta.sub <- meta.sub[order(match(meta.sub$Strain, chrTree$tip.label)), ]
meta.sub$treetip <- chrTree$tip.label
chrTree$tip.label <- meta.sub$CHR_haplotype

meta.sub <- meta[which(meta$Strain %in% symTree$tip.label), c('SYM_haplotype', 'Strain')]
meta.sub <- meta.sub[order(match(meta.sub$Strain, symTree$tip.label)), ]
meta.sub$treetip <- symTree$tip.label
symTree$tip.label <- meta.sub$SYM_haplotype


# Distance matrix
phydistChr <- cophenetic(chrTree)
phydistSym <- cophenetic(symTree)

## Community by host
commHostCHR <- as.matrix(table(meta$Host, meta$CHR_haplotype))
commHostSYM <- as.matrix(table(meta$Host, meta$SYM_haplotype))
commHostCHR <- commHostCHR[,which(colnames(commHostCHR) %in% chrTree$tip.label)]
commHostSYM <- commHostSYM[,which(colnames(commHostSYM) %in% symTree$tip.label)]

## Community by site
commSiteCHR <- as.matrix(table(meta$Site, meta$CHR_haplotype))
commSiteSYM <- as.matrix(table(meta$Site, meta$SYM_haplotype))
commSiteCHR <- commSiteCHR[,which(colnames(commSiteCHR) %in% chrTree$tip.label)]
commSiteSYM <- commSiteSYM[,which(colnames(commSiteSYM) %in% symTree$tip.label)]


# Order community by host
commOrderChr <- commHostCHR[, chrTree$tip.label]
commOrderSym <- commHostSYM[, symTree$tip.label]

# Viz CHR
par(mfrow = c(3, 3))
for (i in row.names(commOrderChr)) {
  plot(chrTree, show.tip.label = FALSE, main = i)
  tiplabels(tip = which(commOrderChr[i, ] > 0), pch = 19, cex = 2)
}
mtext('CHR haplotype distribution by host', side=1, line = 1, adj=3)

# Viz SYM
par(mfrow = c(3, 3))
for (i in row.names(commOrderSym)) {
  plot(symTree, show.tip.label = FALSE, main = i)
  tiplabels(tip = which(commOrderSym[i, ] > 0), pch = 19, cex = 2)
}
mtext('SYM haplotype distribution by host', side=1, line = 1, adj=3)




# Analyze CHR x Host
mpdChrHost <- ses.mpd(commOrderChr, phydistChr,null.model="taxa.labels", abundance.weighted = T)
mntdChrHost <- ses.mntd(commOrderChr, phydistChr,null.model="taxa.labels", abundance.weighted = T)

# Analyze SYM x Host
mpdSymHost <- ses.mpd(commOrderSym, phydistSym,null.model="taxa.labels", abundance.weighted = T)
mntdSymHost <- ses.mntd(commOrderSym, phydistSym,null.model="taxa.labels", abundance.weighted = T)


###########
# By Site #
###########

# Order 
commOrderChr <- commSiteCHR[, chrTree$tip.label]
commOrderSym <- commSiteSYM[, symTree$tip.label]

# Viz CHR
par(mfrow = c(5, 4))
for (i in row.names(commOrderChr)) {
  plot(chrTree, show.tip.label = FALSE, main = i)
  tiplabels(tip = which(commOrderChr[i, ] > 0), pch = 19, cex = 2)
}
mtext('CHR haplotype distribution by sampling site', side=1, line = 1, adj=-1)

# Viz SYM
par(mfrow = c(5, 4))
for (i in row.names(commOrderSym)) {
  plot(symTree, show.tip.label = FALSE, main = i)
  tiplabels(tip = which(commOrderSym[i, ] > 0), pch = 19, cex = 2)
}
mtext('SYM haplotype distribution by sampling site', side=1, line = 1, adj=-1)


# Analyze CHR x Site
mpdChrSite <- ses.mpd(commOrderChr, phydistChr,null.model="taxa.labels", abundance.weighted = T)
mntdChrSite <- ses.mntd(commOrderChr, phydistChr,null.model="taxa.labels", abundance.weighted = T)

# Analyze SYM x Site
mpdSymSite <- ses.mpd(commOrderSym, phydistSym,null.model="taxa.labels", abundance.weighted = T)
mntdSymSite <- ses.mntd(commOrderSym, phydistSym,null.model="taxa.labels", abundance.weighted = T)

###############################
# Phylogenetic beta diversity #
###############################
library(cluster)
par(mfrow = c(1, 1))

# Beta diversity by Host
comdistHostChr <- comdist(commHostCHR, phydistChr)
comdistHosstChr.cluster <- hclust(comdistHostChr)
plot(comdistHosstChr.cluster)


comdistHostSym <- comdist(commHostSYM, phydistSym)
comdistHostSym.cluster <- hclust(comdistHostSym)
plot(comdistHostSym.cluster)



# Beta diversity by Sampling Site
comdistSiteChr <- comdist(commSiteCHR, phydistChr)
comdistSiteChr.cluster <- hclust(comdistSiteChr)
plot(comdistSiteChr.cluster)


comdistSiteSym <- comdist(commSiteSYM, phydistSym)
comdistSiteSym.cluster <- hclust(comdistSiteSym)
plot(comdistSiteSym.cluster)




###############
# Doing AMOVA #
###############

library(BiodiversityR) # also loads vegan
library(poppr) 
library(ggplot2)
library(ggsci)
library(ggforce)
library(dplyr)
library(ggrepel)
library(ape)

#######
# CHR #
#######

# Genclone object 
dna <- read.dna("~/GDrive/Sachs/Chapter3_Epidemic_Genotype/Epidemic Genotype/Sequences/Bradyrhizobium/Dataset/Complete/Dataset/chr/concatenate/aln/chr.meta.concatenated.fasta", format = "fasta")
pop <- DNAbin2genind(dna)

# Meta file
meta <- read.csv('~/GDrive/Sachs/Chapter3_Epidemic_Genotype/Epidemic Genotype/Sequences/Bradyrhizobium/Dataset/Complete/Dataset/meta_all.csv', header=T)
meta <- meta[,1:28]

# Matching strain names from meta to  the fasta file
meta$Sample<- meta$Strain
meta <- meta[meta$Strain %in% indNames(pop), ]
mat <- match(indNames(pop), meta$Sample)
meta <- meta[mat, ]

# Assigning population structure from meta
strata(pop) <- meta[, c('Site',   'Host', 'CHR_haplotype', 'SYM_haplotype')]

# Define the population
setPop(pop) <- ~ Site

# Convert to GenClone object
pop <- as.genclone(pop)

# AMOVA
amova <- poppr.amova(pop, ~Site/Host, clonecorrect = T)
amovaSig <-randtest(amova, nrepet = 999) 
plot(amovaSig)

# Randomize population structure
pop.rand <- pop
strata(pop.rand) <- strata(pop)[sample(nInd(pop)),]

amova.rand <- poppr.amova(pop.rand, ~Site/Host, clonecorrect=T)
amovaSig.rand <-randtest(amova.rand, nrepet = 999) 
plot(amovaSig.rand)

#######
# SYM #
#######

# Genclone object 
dna <- read.dna("~/GDrive/Sachs/Chapter3_Epidemic_Genotype/Epidemic Genotype/Sequences/Bradyrhizobium/Dataset/Complete/Dataset/sym/concatenate/aln/sym.meta.concatenated.fasta", format = "fasta")
pop <- DNAbin2genind(dna)

# Meta file
meta <- read.csv('~/GDrive/Sachs/Chapter3_Epidemic_Genotype/Epidemic Genotype/Sequences/Bradyrhizobium/Dataset/Complete/Dataset/meta_all.csv', header=T)
meta <- meta[,1:28]

# Matching strain names from meta to  the fasta file
meta$Sample<- meta$Strain
meta <- meta[meta$Strain %in% indNames(pop), ]
mat <- match(indNames(pop), meta$Sample)
meta <- meta[mat, ]

# Assigning population structure from meta
strata(pop) <- meta[, c('Site',   'Host', 'CHR_haplotype', 'SYM_haplotype')]

# Define the population
setPop(pop) <- ~ Site

# Convert to GenClone object
pop <- as.genclone(pop)

# AMOVA
amova <- poppr.amova(pop, ~Site/Host, clonecorrect=T)
amovaSig <-randtest(amova, nrepet = 999) 
plot(amovaSig)

# Randomize population structure
pop.rand <- pop
strata(pop.rand) <- strata(pop)[sample(nInd(pop)),]

amova.rand <- poppr.amova(pop.rand, ~Site/Host, clonecorrect=T)
amovaSig.rand <-randtest(amova.rand, nrepet = 999) 
plot(amovaSig.rand)


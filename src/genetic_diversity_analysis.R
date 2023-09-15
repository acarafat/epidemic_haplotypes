#####################
# Epidemic Genotype #
#####################

# 20 May 2022
# Final set of analysis required to complete the project:
# - Fully resolved phylogenetic tree, Show other meta feature by side [DONE]
# - Calculate overall nucleotide diversity [DONE]
# - SYM clad based nucleotide diversity [DONE]
# - Selection analysis on the SYM loci [DONE]
# - Summarizing selection analysis data
# - Calculate Tajima's D [DONE]
# - PCA like analysis to visualize clustering [DONE]
# - Clustering pattern by site/host (PCA/discriminant analysis) [DONE]
# - Plotting POCP/ANI data
# - Characterizing symICE types

###################################################################################
# Phylogenetic tree visualization of all isolates for both CHR and SYM haplotypes #
###################################################################################
library(ggtree)
library(phangorn)
library(ggnewscale)
library(treeio)
library(RColorBrewer)
library(ggpubr)

chr.tree <- read.iqtree('/Users/arafat/GDrive/Sachs/Chapter3_Epidemic_Genotype/Sequences/Bradyrhizobium/Dataset/Complete/Dataset/chr/concatenate/aln/chr.concatenated.partitions.rxml.treefile')
chr.tree@phylo <- midpoint(chr.tree@phylo)

sym.tree <- read.iqtree('/Users/arafat/GDrive/Sachs/Chapter3_Epidemic_Genotype/Sequences/Bradyrhizobium/Dataset/Complete/Dataset/sym/concatenate/aln/sym.concatenated.partitions.rxml.treefile')
sym.tree@phylo <- midpoint(sym.tree@phylo)
#sym.tree <- midpoint(sym.tree)

# anotate tree figures to visualize the other meta information
meta <- read.table('/Users/arafat/GDrive/Sachs/Chapter3_Epidemic_Genotype/Sequences/Bradyrhizobium/Dataset/Complete/Dataset/meta_clade.csv', sep=',', header=T)

# associated matrix
meta_matrix <- meta[,c('Strain', 'Site', 'Host', 'CHR.clade')]
rownames(meta_matrix) <- meta_matrix$Strain
meta_matrix <- meta_matrix[,-1]

# - Show bootstrap value, need a way to incorporate sym.tree$node.label
# - Add Site+Host info

# Ornamentation of the tree
t1 <- ggtree(chr.tree) %<+% meta
t2 <- ggtree(sym.tree) %<+% meta  #+  
  #geom_tippoint(aes(color = Site)) #+
  #geom_tiplab(aes(label=SYM_haplotype), size=1, check.overlap=T) 

# Show the internal nodes (which is necessary to find out the clades)    
t2 + geom_text(aes(label=node), hjust=-.3)

# Visualize the clades of interest
t2 + geom_hilight(node=563) + geom_hilight(node=569) +
  geom_hilight(node=631) + geom_hilight(node=671) + geom_hilight(node=750) +
  geom_cladelabel(node=563, label="SYM clade 1", align=TRUE,  color='red')

# Annotate the clades
t2 + geom_cladelabel(node=563, label="SYM 1", align=TRUE, offset=0.001,  color='red') +
  geom_cladelabel(node=570, label="SYM 2", align=TRUE, offset=-0.02,  color='blue') +
  geom_cladelabel(node=613, label="SYM 3", align=TRUE, offset=-0.02,  color='brown') +
  geom_cladelabel(node=632, label="SYM 4", align=TRUE, offset=-0.02,  color='green') +
  geom_cladelabel(node=671, label="SYM 5", align=TRUE, offset=-0.005,  color='orange')

# Get the isolate ids from each SYM clades 
cladeSYM1 <- get_taxa_name(t2, 563)
cladeSYM2 <- get_taxa_name(t2, 570)
cladeSYM3 <- get_taxa_name(t2, 613)
cladeSYM4 <- get_taxa_name(t2, 632)
cladeSYM5 <- get_taxa_name(t2, 671)

# Plot tree with associated matrix
meta_matrix.site <- as.data.frame(meta_matrix[,'Site'])
meta_matrix.host <- as.data.frame(meta_matrix[,'Host'])
meta_matrix.clade <- as.data.frame(meta_matrix[,'CHR.clade'])

colnames(meta_matrix.site) <- 'Sample\nSite'
colnames(meta_matrix.host) <- 'Host'
colnames(meta_matrix.clade) <- 'Clade'

rownames(meta_matrix.site) <- rownames(meta_matrix)
rownames(meta_matrix.host) <- rownames(meta_matrix)
rownames(meta_matrix.clade) <- rownames(meta_matrix)

# Plot tree for CHR 
t1 <- ggtree(chr.tree) %<+% meta
t1$data$bootstrap <- '0'
t1$data[which(t1$data$SH_aLRT >= 70 & t1$data$UFboot  >= 70),]$bootstrap <- '1'

t1_tips <- t1 + geom_tree(aes(color=bootstrap == '1')) +
  scale_color_manual(name='Bootstrap', values=setNames( c('black', 'grey'), c(T,F)), guide = "none") +
  #geom_tippoint(data=t1$data[which(t1$data$Host != 'A. strigosus'),] , aes(fill=Host), shape=21) +
  geom_tippoint(data=t1$data, aes(fill=CHR.clade), shape=21)  + new_scale_fill() 

t1_mat1 <- gheatmap(t1_tips, meta_matrix.site, width=0.1,
                    colnames_offset_y=5, colnames_position='top', custom_column_labels='Population') +
  scale_fill_viridis_d(option="D", name="Population") +
  coord_cartesian(clip = "off") +  new_scale_fill()

gheatmap(t1_mat1, meta_matrix.host, width=0.1, offset=0.0125,
         colnames_offset_y=5, colnames_position='top', custom_column_labels='Host') +
  scale_fill_viridis_d(option="B", name="Host") +
  scale_y_continuous(expand=c(0.03, 0.3)) +
  coord_cartesian(clip = "off")


#t1_mat2 <- gheatmap(t1_mat1, meta_matrix.host, width=0.1, offset=0.0125,
#         colnames_angle=0, legend_title="Host") +
#  scale_fill_viridis_d(option="C", name="Host") +
#  new_scale_fill()

# Plot tree for SYM 
t2$data$bootstrap <- '0'
t2$data[which(t2$data$SH_aLRT >= 70 & t2$data$UFboot  >= 70),]$bootstrap <- '1'


t2_tips <- t2 + geom_tree(aes(color=bootstrap == '1')) +
  scale_color_manual(name='Bootstrap', values=setNames( c('black', 'grey'), c(T,F)), guide = "none") +
  #geom_tippoint(data=t2$data[which(t2$data$Host != 'A. strigosus'),] , aes(fill=Host), shape=21) +
  geom_tippoint(data=t2$data, aes(fill=CHR.clade), shape=21) +
  new_scale_fill()


t2_mat1 <- gheatmap(t2_tips, meta_matrix.site, width=0.1, 
                    colnames_offset_y=8, legend_title="Population", colnames_position='top') +
  scale_fill_viridis_d(option="D", name="Population") +
  new_scale_fill()

gheatmap(t2_mat1, meta_matrix.host, width=0.1, offset=0.006,
         colnames_offset_y=2.5, legend_title="Host", colnames_position='top') +
  scale_fill_viridis_d(option="B", name="Host") +
  scale_y_continuous(expand=c(0.03, 0.3))

######################################################
# Same visualization scheme, but flip clade and host #
######################################################
# Plot tree for CHR 
t1_tips <- t1 + geom_tree(aes(color=bootstrap == '1')) +
  scale_color_manual(name='Bootstrap', values=setNames( c('black', 'grey'), c(T,F)), guide = "none") +
  geom_tippoint(data=t1$data, aes(fill=CHR.clade), shape=21) +
  labs(fill='Bradyrhizobium\nspecies') + new_scale_fill() +
  geom_treescale( y = -5, offset=1)  

t1_mat1 <- gheatmap(t1_tips, meta_matrix.site, width=0.1, 
                    colnames_offset_y=10, legend_title="Sample Site", colnames_position='top', font.size=3) +
  scale_fill_viridis_d(option="D", name="Sample Site") +
  new_scale_fill()

t1_mat1 <- gheatmap(t1_mat1, meta_matrix.host, width=0.1, offset=0.0125,
         colnames_offset_y=5, legend_title="Host", colnames_position='top', font.size = 3) +
  scale_fill_viridis_d(option="B", name="Host") +
  scale_y_continuous(expand=c(0.03, 0.3)) +
  ggtitle('A. Chromosomal haplotype phylogeny') 

# Plot tree for SYM 
t2_tips <- t2 + geom_tree(aes(color=bootstrap == '1')) +
  scale_color_manual(name='Bootstrap', values=setNames( c('black', 'grey'), c(T,F)), guide = "none") +
  geom_tippoint(data=t2$data, aes(fill=CHR.clade), shape=21) +
  labs(fill='Bradyrhizobium\nspecies') + new_scale_fill() +
  geom_treescale( y = -7, offset=1)  


t2_mat1 <- gheatmap(t2_tips, meta_matrix.site, width=0.1, 
                    colnames_offset_y=10, legend_title="Sample Site", colnames_position='top', font.size = 3) +
  scale_fill_viridis_d(option="D", name="Sample Site") +
  new_scale_fill()

t2_mat1 <- gheatmap(t2_mat1, meta_matrix.host, width=0.1, offset=0.006,
         colnames_offset_y=5, legend_title="Host", colnames_position='top', font.size = 3) +
  scale_fill_viridis_d(option="B", name="Host") +
  scale_y_continuous(expand=c(0.03, 0.3)) +
  ggtitle('B. symICE haplotype phylogeny')



t2_mat1 <- t2_mat1 + geom_cladelabel(node=563, label="SYM 1", align=TRUE, offset=-0.012,  color='red') +
  geom_cladelabel(node=570, label="SYM 2", align=TRUE, offset=-0.02,  color='blue') +
  geom_cladelabel(node=613, label="SYM 3", align=TRUE, offset=-0.02,  color='brown') +
  geom_cladelabel(node=632, label="SYM 4", align=TRUE, offset=-0.02,  color='green') +
  geom_cladelabel(node=671, label="SYM 5", align=TRUE, offset=-0.01,  color='orange')

ggarrange(t1_mat1, t2_mat1, ncol=2, common.legend = TRUE, legend="right")

########################
# Calculate Tajima's D #
########################
library(pegas)

chr.seq <- fasta2DNAbin('/Users/arafat/GDrive/Sachs/Chapter3_Epidemic_Genotype/Sequences/Bradyrhizobium/Dataset/Complete/Dataset/chr/concatenate/aln/chr.meta.concatenated.fasta')
sym.seq <- fasta2DNAbin('/Users/arafat/GDrive/Sachs/Chapter3_Epidemic_Genotype/Sequences/Bradyrhizobium/Dataset/Complete/Dataset/sym/concatenate/aln/sym.meta.concatenated.fasta')

tajima.test(chr.seq)
tajima.test(sym.seq)

# CHR genes
# dnaK, ITS, recA, glnII
dnaK.seq <- fasta2DNAbin('/Users/arafat/GDrive/Sachs/Chapter3_Epidemic_Genotype/Sequences/Bradyrhizobium/Dataset/Complete/Dataset/chr/concatenate/aln/dnaK.567.fasta.filter.fasta.aln')
its.seq <- fasta2DNAbin('/Users/arafat/GDrive/Sachs/Chapter3_Epidemic_Genotype/Sequences/Bradyrhizobium/Dataset/Complete/Dataset/chr/concatenate/aln/ITS.567.fasta.filter.fasta.aln')
recA.seq <- fasta2DNAbin('/Users/arafat/GDrive/Sachs/Chapter3_Epidemic_Genotype/Sequences/Bradyrhizobium/Dataset/Complete/Dataset/chr/concatenate/aln/recA.567.fasta.filter.fasta.aln')
glnII.seq <- fasta2DNAbin('/Users/arafat/GDrive/Sachs/Chapter3_Epidemic_Genotype/Sequences/Bradyrhizobium/Dataset/Complete/Dataset/chr/concatenate/aln/glnII.567.fasta.filter.fasta.aln')

tajima.test(dnaK.seq)
tajima.test(its.seq)
tajima.test(recA.seq)
tajima.test(glnII.seq)

nuc.div(dnaK.seq)
nuc.div(its.seq)
nuc.div(recA.seq)
nuc.div(glnII.seq)


# SYM genes
# nodZ, nodL, nodA, nifD
nodZ.seq <- fasta2DNAbin('/Users/arafat/GDrive/Sachs/Chapter3_Epidemic_Genotype/Sequences/Bradyrhizobium/Dataset/Complete/Dataset/sym/concatenate/aln/nodZ.brady.fasta.filter.fasta.aln')
nodL.seq <- fasta2DNAbin('/Users/arafat/GDrive/Sachs/Chapter3_Epidemic_Genotype/Sequences/Bradyrhizobium/Dataset/Complete/Dataset/sym/concatenate/aln/nodL.brady.fasta.filter.fasta.aln')
nodA.seq <- fasta2DNAbin('/Users/arafat/GDrive/Sachs/Chapter3_Epidemic_Genotype/Sequences/Bradyrhizobium/Dataset/Complete/Dataset/sym/concatenate/aln/nodA.brady.fasta.filter.fasta.aln')
nifD.seq <- fasta2DNAbin('/Users/arafat/GDrive/Sachs/Chapter3_Epidemic_Genotype/Sequences/Bradyrhizobium/Dataset/Complete/Dataset/sym/concatenate/aln/nifD.brady.fasta.filter.fasta.aln')

tajima.test(nodZ.seq)
tajima.test(nodL.seq)
tajima.test(nodA.seq)
tajima.test(nifD.seq)

nuc.div(nodZ.seq)
nuc.div(nodL.seq)
nuc.div(nodA.seq)
nuc.div(nifD.seq)

##########################################
# Calculate overall nucleotide diversity #
##########################################
library(pegas)
library(adegenet)
library(reshape2)
library(ggplot2)

chr.seq <- fasta2DNAbin('/Users/arafat/GDrive/Sachs/Chapter3_Epidemic_Genotype/Sequences/Bradyrhizobium/Dataset/Complete/Dataset/chr/concatenate/aln/chr.meta.concatenated.fasta')
sym.seq <- fasta2DNAbin('/Users/arafat/GDrive/Sachs/Chapter3_Epidemic_Genotype/Sequences/Bradyrhizobium/Dataset/Complete/Dataset/sym/concatenate/aln/sym.meta.concatenated.fasta')

# Es

# convert DNAbin to R matrix, subset, and convert back to DNAbin object
subsetCHR <- as.character(chr.seq)
subsetSYM <- as.character(sym.seq)
subsetCHR[1:5, 1:5]

# Subset both CHR and SYM sequences based on SYM clade membership
sym1.chr <- as.DNAbin(subsetCHR[which(rownames(subsetCHR) %in% cladeSYM1), ])
sym2.chr <- as.DNAbin(subsetCHR[which(rownames(subsetCHR) %in% cladeSYM2), ])
sym3.chr <- as.DNAbin(subsetCHR[which(rownames(subsetCHR) %in% cladeSYM3), ])
sym4.chr <- as.DNAbin(subsetCHR[which(rownames(subsetCHR) %in% cladeSYM4), ])
sym5.chr <- as.DNAbin(subsetCHR[which(rownames(subsetCHR) %in% cladeSYM5), ])

sym1.sym <- as.DNAbin(subsetSYM[which(rownames(subsetSYM) %in% cladeSYM1), ])
sym2.sym <- as.DNAbin(subsetSYM[which(rownames(subsetSYM) %in% cladeSYM2), ])
sym3.sym <- as.DNAbin(subsetSYM[which(rownames(subsetSYM) %in% cladeSYM3), ])
sym4.sym <- as.DNAbin(subsetSYM[which(rownames(subsetSYM) %in% cladeSYM4), ])
sym5.sym <- as.DNAbin(subsetSYM[which(rownames(subsetSYM) %in% cladeSYM5), ])

# Total genetic diversity
c.all <- nuc.div(chr.seq)
s.all <- nuc.div(sym.seq)


c1 <- nuc.div(sym1.chr)
s1 <- nuc.div(sym1.sym)

c2 <- nuc.div(sym2.chr)
s2 <- nuc.div(sym2.sym)

c3 <- nuc.div(sym3.chr)
s3 <- nuc.div(sym3.sym)

c4 <- nuc.div(sym4.chr)
s4 <- nuc.div(sym4.sym)

c5 <- nuc.div(sym5.chr)
s5 <- nuc.div(sym5.sym)

compare_nucdiv <- data.frame(clade=c('Whole Tree','SYM 1 Clade', 'SYM 2 Clade', 'SYM 3 Clade', 'SYM 4 Clade', 'SYM 5 Clade'),
           CHR=c(c.all, c1, c2, c3, c4, c5),
           SYM=c(s.all, s1, s2, s3, s4, c5))ß

compare_nucdiv <- melt(compare_nucdiv, id='clade')

ggplot(compare_nucdiv, aes(x=clade, y=value, fill=variable)) + 
  geom_bar(stat='identity', position='dodge') +ß
  xlab('Comparison level') + ylab('Nucleotide Diversity') +
  scale_fill_brewer(palette = "Paired", label=c('Chromosomal haplotypes', 'symICEs')) + 
  theme(legend.position = 'top', legend.title = element_blank(),
        axis.text.x=element_text(angle=30, vjust=0.5, hjust=0.5)) +
  ggtitle('Nucleotide diversity comparison based on symICE phylogeney clades')

# Why SYM 2 clade has opposite trend? i.e. higher diversity in CHR tree to SYM ?

t1 + geom_tippoint(aes(color = label %in% cladeSYM5)) +
  scale_color_manual(name='Member of SYM clade', values=setNames( c('steelblue', 'grey'), c(T, F)))
t2 + geom_tippoint(aes(color = label %in% cladeSYM5)) +
  scale_color_manual(name='Member of SYM clade', values=setNames( c('steelblue', 'grey'), c(T, F)))

# Estimate and compare genetic diversity of the clades



####################################################
# Nucleotide diversity by syn/nonsynonymous subset #
####################################################
library(PopGenome)

# Popgenome requires individual sequences in a folder. 
# So use a Python script to subset each sequence from master fasto into individual fasta file
GENOME.class <- readData('/Users/arafat/GDrive/Sachs/Chapter3_Epidemic_Genotype/Sequences/Bradyrhizobium/Dataset/Complete/Dataset/individual_seqs/CHR/')
#GENOME.class.split <- splitting.data(GENOME.class, subsites="coding")
GENOME.class.split <- neutrality.stats(GENOME.class, subsites="nonsyn")
GENOME.class.split@region.stats@nucleotide.diversity


#######
# PCA #
#######
library("adegenet")
library("ggplot2")

chr.snp <- fasta2genlight('/Users/arafat/GDrive/Sachs/Chapter3_Epidemic_Genotype/Sequences/Bradyrhizobium/Dataset/Complete/Dataset/chr/concatenate/aln/chr.meta.concatenated.fasta', snpOnly=T)
sym.snp <- fasta2genlight('/Users/arafat/GDrive/Sachs/Chapter3_Epidemic_Genotype/Sequences/Bradyrhizobium/Dataset/Complete/Dataset/sym/concatenate/aln/sym.meta.concatenated.fasta', snpOnly=T)
meta <- read.table('/Users/arafat/GDrive/Sachs/Chapter3_Epidemic_Genotype/Sequences/Bradyrhizobium/Dataset/Complete/Dataset/meta_clade.csv', sep=',', header=T)

epidemic <- c('K1_R1_I1_G1', 'K1_R1_I3_G1', 'K1_R3_I1_G1')
dominant <- c('K1_R1_I1_G1', 'K1_R1_I3_G1', 'K1_R3_I1_G1', 'K1_R1_I6_G1', 'K18_R1_I1_G1', 'K21_R23_I55_G4')

# Update SYM lineage
meta$SYM.lineage <- 'undefined'

meta[which(meta$Strain %in% cladeSYM1),]$SYM.lineage <- '1'
meta[which(meta$Strain %in% cladeSYM2),]$SYM.lineage <- '2'
meta[which(meta$Strain %in% cladeSYM3),]$SYM.lineage <- '3'
meta[which(meta$Strain %in% cladeSYM4),]$SYM.lineage <- '4'
meta[which(meta$Strain %in% cladeSYM5),]$SYM.lineage <- '5'

chr.snp$pop = as.factor(meta[match(chr.snp$ind.names, meta$Strain),]$Site)
sym.snp$pop = as.factor(meta[match(sym.snp$ind.names, meta$Strain),]$Site)

# PCA
chr.pca <- glPca(chr.snp, nf=10)
sym.pca <- glPca(sym.snp, nf=10)

# Plot the individuals in SNP-PCA space, with locality labels:
plot(chr.pca$scores[,1], chr.pca$scores[,2], 
     cex=2, pch=20, col=chr.snp$pop, 
     xlab="Principal Component 1", 
     ylab="Principal Component 2", 
     main="PCA on CHR data")
# Plot the individuals in SNP-PCA space, with locality labels:
plot(sym.pca$scores[,1], sym.pca$scores[,2], 
     cex=2, pch=20, col=sym.snp$pop, 
     xlab="Principal Component 1", 
     ylab="Principal Component 2", 
     main="PCA on SYM data")

# Which SNPs load most strongly on the 1st PC axis?
loadingplot(abs(chr.pca$loadings[,1]),
            threshold=quantile(abs(chr.pca$loadings), 0.999))

loadingplot(abs(sym.pca$loadings[,1]),
            threshold=quantile(abs(sym.pca$loadings), 0.999))

# Get their locus names
chr.snp$loc.names[which(quantile(abs(chr.pca$loadings))>0.999)]
sym.snp$loc.names[which(quantile(abs(sym.pca$loadings))>0.999)]

# Decorate+Visualize CHR pca
chr.pca.dataset = as.data.frame(chr.pca$scores)
chr.pca.dataset$isolates = rownames(chr.pca.dataset)
chr.pca.dataset$pop = as.factor(meta[match(chr.snp$ind.names, meta$Strain),]$Site)
chr.pca.dataset$host = as.factor(meta[match(chr.snp$ind.names, meta$Strain),]$Host)
chr.pca.dataset$clade = as.factor(meta[match(chr.snp$ind.names, meta$Strain),]$CHR.clade)
chr.pca.dataset$sym = as.factor(meta[match(chr.snp$ind.names, meta$Strain),]$SYM.lineage)
chr.pca.dataset$CHR_haplotype = as.factor(meta[match(chr.snp$ind.names, meta$Strain),]$CHR_haplotype)
chr.pca.dataset$SYM_haplotype = as.factor(meta[match(chr.snp$ind.names, meta$Strain),]$SYM_haplotype)

ggplot(chr.pca.dataset, aes(PC1, PC2, fill=host)) + geom_point(shape=21, size=3, alpha=0.7)
ggplot(chr.pca.dataset, aes(PC1, PC2, fill=pop)) + geom_point(shape=21, size=3, alpha=0.7)
ggplot(chr.pca.dataset, aes(PC1, PC2, fill=clade)) + geom_point(shape=21, size=3, alpha=0.7)
ggplot(chr.pca.dataset, aes(PC1, PC2, fill=sym)) + geom_point(shape=21, size=3, alpha=0.7)
ggplot(chr.pca.dataset, aes(PC1, PC2, fill=clade)) + geom_point(shape=21, size=3, alpha=0.7) +
  stat_ellipse(level = 0.95, size = 1) +
  geom_point(data=chr.pca.dataset[which(chr.pca.dataset$CHR_haplotype == 'K1_R1_I1_G1'),], aes(x=PC1, y=PC2), size=5, shape=24) 

ggplot(chr.pca.dataset, aes(PC1, PC2, fill=clade)) + geom_point(shape=21, size=3, alpha=0.7) + 
  stat_ellipse(level = 0.95, size = 1) +
  geom_point(data=chr.pca.dataset[which(chr.pca.dataset$CHR_haplotype %in% epidemic),], aes(x=PC1, y=PC2, fill=clade), size=5, shape=24, alpha=0.7)

ggplot(chr.pca.dataset, aes(PC1, PC2, fill=clade)) + 
  geom_point(data=chr.pca.dataset[which(chr.pca.dataset$CHR_haplotype %in% dominant),], aes(x=PC1, y=PC2), color='red', size=5, alpha=7, show.legend = F) +
  geom_point(shape=21, size=3, alpha=0.7) +
  stat_ellipse(level = 0.95, size = 0.7, show.legend = F, aes(color=clade)) +
  labs(fill='Bradyrhizobium Spp.')



# Decorate+Visualize SYM pca
sym.pca.dataset = as.data.frame(sym.pca$scores)
sym.pca.dataset$isolates = rownames(sym.pca.dataset)
sym.pca.dataset$pop = as.factor(meta[match(sym.snp$ind.names, meta$Strain),]$Site)
sym.pca.dataset$clade = as.factor(meta[match(sym.snp$ind.names, meta$Strain),]$CHR.clade)
sym.pca.dataset$sym = as.factor(meta[match(sym.snp$ind.names, meta$Strain),]$SYM.lineage)
sym.pca.dataset$CHR_haplotype = as.factor(meta[match(sym.snp$ind.names, meta$Strain),]$CHR_haplotype)
sym.pca.dataset$SYM_haplotype = as.factor(meta[match(sym.snp$ind.names, meta$Strain),]$SYM_haplotype)
sym.pca.dataset$host = as.factor(meta[match(sym.snp$ind.names, meta$Strain),]$Host)

ggplot(sym.pca.dataset, aes(PC1, PC2, fill=host)) + geom_point(shape=21, size=3, alpha=0.7)
ggplot(sym.pca.dataset, aes(PC1, PC2, fill=pop)) + geom_point(shape=21, size=3, alpha=0.7)
ggplot(sym.pca.dataset, aes(PC1, PC2, fill=clade)) + geom_point(shape=21, size=3, alpha=0.7)
ggplot(sym.pca.dataset, aes(PC1, PC2, fill=sym)) + geom_point(shape=21, size=3, alpha=0.7)

ggplot(sym.pca.dataset, aes(PC1, PC2, fill=clade)) + geom_point(shape=21, size=3, alpha=0.7) + 
  geom_point(data=sym.pca.dataset[which(sym.pca.dataset$CHR_haplotype == 'K1_R1_I1_G1'),], aes(x=PC1, y=PC2), size=5, shape=24, alpha=0.7)

ggplot(sym.pca.dataset, aes(PC1, PC2, fill=clade)) + geom_point(shape=21, size=3, alpha=0.7) + 
  geom_point(data=sym.pca.dataset[which(sym.pca.dataset$CHR_haplotype %in% epidemic),], aes(x=PC1, y=PC2, fill=clade), size=5, shape=24, alpha=0.7)

ggplot(sym.pca.dataset, aes(PC1, PC2, fill=clade)) + geom_point(shape=21, size=3, alpha=0.7) + 
  geom_point(data=sym.pca.dataset[which(sym.pca.dataset$CHR_haplotype %in% dominant),], aes(x=PC1, y=PC2,  shape=CHR_haplotype), size=3, alpha=0.5)

ggplot(sym.pca.dataset, aes(PC1, PC2, fill=clade)) + 
  geom_point(data=sym.pca.dataset[which(sym.pca.dataset$CHR_haplotype %in% dominant),], aes(x=PC1, y=PC2), color='red', size=5, alpha=7, show.legend = F) +
  geom_point(shape=21, size=3, alpha=0.7) +
  labs(fill='Bradyrhizobium Spp.')
  #stat_ellipse(level = 0.95, size = 0.7)

#######################################
# PD of CHR associated SYM haplotypes #
#######################################
library(geiger)

meta <- read.table('~/GDrive/Sachs/Chapter3_Epidemic_Genotype/Sequences/Bradyrhizobium/Dataset/Complete/Dataset/meta_all.csv', sep=',', header = T)
sym.tree <- read.tree('~/GDrive/Sachs/Chapter3_Epidemic_Genotype/Sequences/Bradyrhizobium/Dataset/Complete/Dataset/sym/concatenate/aln/sym.concatenated.partitions.rxml.treefile')
sym.tree <- midpoint(sym.tree)

dom.haplotypes <- c('K1_R1_I1_G1', 'K1_R1_I3_G1', 'K1_R3_I1_G1',
                    'K1_R1_I6_G1', 'K18_R1_I1_G1', 'K21_R23_I55_G4')

# PD for K1_R1_I1_G1 
rownames(meta) <- meta$Strain

for (i in dom.haplotypes) {
  meta.h1 <- meta[which(meta$CHR_haplotype == i), ]
  pruned.h1 <- treedata(sym.tree, meta.h1, warnings=F)$phy
  print(sum(pruned.h1$edge.length))
}

# test correlation
cor <- data.frame(association = c(70, 11, 5, 3, 3, 2),
                  pd = c(0.1953965, 0.1079797, 0.07129972, 0.001160326,0.01119363, 0.0003878465))

p <- ggplot(cor, aes(x=association, y=pd)) +
  geom_point() + geom_smooth(method='lm') +
  theme_bw() +
  xlab('#Association') + ylab('SYM diversity')

model = lm(pd~association, data=cor)
summary(model)  

#########################
# Genome stat: ANI+POCP #
#########################
library(ggplot2)
library(gridExtra)

ani.popspp = read.csv('~/GDrive/Sachs/Chapter3_Epidemic_Genotype/Genome_stat/popsppANI.csv')
ani.haplo = read.csv('~/GDrive/Sachs/Chapter3_Epidemic_Genotype/Genome_stat/withinHaploANI.csv')

colnames(ani.haplo)[5:6] <- c('population', 'iteration')
colnames(ani.popspp)[5:6] <- c('population', 'iteration')

ani_merged <- rbind(ani.haplo, ani.popspp)
ani.summary <- aggregate(ani~population, data=ani_merged, FUN=mean)
ani.se <- aggregate(ani~population, data=ani_merged, FUN=se)
colnames(ani.se)[2] <- 'se'

ani.summary <- merge(ani.summary, ani.se, by='population')

ani.summary$population <- factor(ani.summary$population, levels = c('K1_R1_I1_G1', 'K1_R1_I3_G1', 'K1_R3_I1_G1', 'species', 'metapop'))
p1 <- ggplot(ani.summary, aes(x=population, y=ani)) + geom_point(stat='identity', color='red') +
  geom_errorbar(aes(x=population, ymin= as.numeric(ani) - as.numeric(se), ymax = as.numeric(ani)+as.numeric(se) ), width=0.4) +
  ggtitle('A. ANI for core genes') +
  xlab('Population' ) + ylab('ANI') + scale_y_continuous(limits=c(90, 100)) +
  theme_bw() + theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))



pocp.chr = read.csv('~/GDrive/Sachs/Chapter3_Epidemic_Genotype/Genome_stat/pangenome_chr_mat_bootstrap.pocp.csv')
pocp.sym = read.csv('~/GDrive/Sachs/Chapter3_Epidemic_Genotype/Genome_stat/pangenome_sym_mat_bootstrap.pocp.csv')

summary.chr = as.data.frame(as.matrix(aggregate(pocp~population, data=pocp.chr, FUN= function(x) c(mean(x), se(x)) )))
colnames(summary.chr) <- c('population', 'pocp_mean', 'se')
summary.chr$population <- c('K1_R1_I1_G1', 'K1_R1_I3_G1', 'metapop', 'K1_R3_I1_G1', 'species')

summary.sym = as.data.frame(as.matrix(aggregate(pocp~population, data=pocp.sym, FUN= function(x) c(mean(x), se(x)) )))
colnames(summary.sym) <- c('population', 'pocp_mean', 'se')
summary.sym$population <- c('K1_R1_I1_G1', 'K1_R1_I3_G1', 'metapop', 'K1_R3_I1_G1', 'species')

summary.chr$population <- factor(summary.chr$population, levels = c('K1_R1_I1_G1', 'K1_R1_I3_G1', 'K1_R3_I1_G1', 'species', 'metapop'))
p2 <- ggplot(summary.chr, aes(x=population, y=as.numeric(pocp_mean))) + geom_point(stat='identity', color='red') +
  geom_errorbar(aes(x=population, ymin= as.numeric(pocp_mean) - as.numeric(se), ymax = as.numeric(pocp_mean)+as.numeric(se) ), width=0.4) +
  ggtitle('B. CHR genes POCP') +
  xlab('Population' ) + ylab('POCP') + scale_y_continuous(limits=c(0.5, 1)) +
  theme_bw() + theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))

summary.sym$population <- factor(summary.chr$population, levels = c('K1_R1_I1_G1', 'K1_R1_I3_G1', 'K1_R3_I1_G1', 'species', 'metapop'))
p3 <- ggplot(summary.sym, aes(x=population, y=as.numeric(pocp_mean))) + geom_point(stat='identity', color='red') +
  geom_errorbar(aes(x=population, ymin= as.numeric(pocp_mean) - as.numeric(se), ymax = as.numeric(pocp_mean)+as.numeric(se) ), width=0.4) +
  ggtitle('C. SYM genes POCP') +
  xlab('Population' ) + ylab('POCP') + scale_y_continuous(limits=c(0.5, 1)) +
  theme_bw() + theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))

grid.arrange(p1, p2, p3, nrow=1)


##########################
# Characterizing symICEs #
##########################
meta <- read.table('/Users/arafat/GDrive/Sachs/Chapter3_Epidemic_Genotype/Sequences/Bradyrhizobium/Dataset/Complete/Dataset/meta_clade.csv', sep=',', header=T)

file <- read.table('/Users/arafat/GDrive/Sachs/Chapter3_Epidemic_Genotype/symICEs/brady_symICE_orthologgroups.txt' ,sep='\t', header=F)
file <- separate(file,col=V2, into=c('isolate','identifier'), sep="_")

head(file)
unique(file$V1)

file <- file[!grepl("B_", file$V1), ]

meta.ices <- merge(meta, file, by.x='Strain', by.y='isolate')

epidemic <- c('K1_R1_I1_G1', 'K1_R1_I3_G1', 'K1_R3_I1_G1')
dominant <- c('K1_R1_I1_G1', 'K1_R1_I3_G1', 'K1_R3_I1_G1', 'K1_R1_I6_G1', 'K18_R1_I1_G1', 'K21_R23_I55_G4')
                
unique(meta.ices[which(meta.ices$CHR_haplotype == 'K1_R1_I1_G1'), ]$V1)
unique(meta.ices[which(meta.ices$CHR_haplotype == 'K1_R1_I3_G1'), ]$V1)
unique(meta.ices[which(meta.ices$CHR_haplotype == 'K1_R3_I1_G1'), ]$V1)

unique(meta.ices[which(meta.ices$CHR_haplotype == 'K1_R1_I6_G1'), ]$V1)
unique(meta.ices[which(meta.ices$CHR_haplotype == 'K18_R1_I1_G1'), ]$V1)
unique(meta.ices[which(meta.ices$CHR_haplotype == 'K21_R23_I55_G4'), ]$V1)

unique(meta.ices[which(!meta.ices$CHR_haplotype %in% dominant), ]$V1)

# Dominant/Epidemic haplotype symICE

############################################################
# T3SS effectors presence absence clustering by haplotypes #
############################################################
library(gplots)

epidemic <- c('K1_R1_I1_G1', 'K1_R1_I3_G1', 'K1_R3_I1_G1')
dominant <- c('K1_R1_I1_G1', 'K1_R1_I3_G1', 'K1_R3_I1_G1', 'K1_R1_I6_G1', 'K18_R1_I1_G1', 'K21_R23_I55_G4')

matrix <- read.csv('/Users/arafat/GDrive/Sachs/Chapter3_Epidemic_Genotype/symICEs/blast_final_output.csv')
rownames(matrix) <- matrix$X

# Convert the coverage values to present/absent scores
matrix[is.na(matrix)] <- 0
matrix[matrix < 70] <- 0
matrix[matrix > 69] <- 1
matrix$X <- rownames(matrix)

meta <- read.table('/Users/arafat/GDrive/Sachs/Chapter3_Epidemic_Genotype/Sequences/Bradyrhizobium/Dataset/Complete/Dataset/meta_clade.csv', sep=',', header=T)

meta.part <- meta[,c('Strain', 'CHR_haplotype', 'CHR.clade', 'Host')]

t3ss <- merge(matrix, meta.part, by.x='X', by.y='Strain')
rownames(t3ss) <- t3ss$X


# PCA
library(ggbiplot)

data <- t3ss[,2:39]
df_f <- data[,apply(data, 2, var, na.rm=TRUE) != 0]

pca <- prcomp(df_f, scale=F, center=T)

pca.coords <- as.data.frame(pca$x)
pca.coords$isolate <- rownames(pca$x)
pca.coords <- merge(pca.coords, meta.part, by.x='isolate', by.y='Strain')
pca.coords$dominant <- 'non-dominant'
pca.coords$epidemic <- 'non-epidemic'
pca.coords[which(pca.coords$CHR_haplotype %in% dominant),]$dominant <- 'dominant'
pca.coords[which(pca.coords$CHR_haplotype %in% epidemic),]$epidemic <- 'epidemic'



ggplot(pca.coords, aes(x=PC1, y=PC2, color=Host)) + geom_point()
ggplot(pca.coords, aes(x=PC1, y=PC2, color=dominant)) + geom_point()
ggplot(pca.coords, aes(x=PC1, y=PC2, color=epidemic)) + geom_point()


# Heatmap
heatmap( as.matrix(t3ss[,2:39]), col = c("grey", "salmon"), Rowv = F)

heatmap( as.matrix(t3ss[which(t3ss$CHR_haplotype %in% epidemic),2:39]), col = c("grey", "salmon"), Rowv = F)

heatmap( as.matrix(t3ss[which(t3ss$CHR_haplotype %in% dominant),2:39]), col = c("grey", "salmon"), Rowv = F)
heatmap( as.matrix(t3ss[which(!t3ss$CHR_haplotype %in% dominant),2:39]), col = c("grey", "salmon"), Rowv = F)

heatmap( as.matrix(t3ss[which(t3ss$CHR_haplotype == 'K1_R1_I1_G1'),2:39]), col = c("grey", "salmon"), Rowv = F)
heatmap( as.matrix(t3ss[which(t3ss$CHR_haplotype == 'K1_R1_I3_G1'),2:39]), col = c("grey", "salmon"), Rowv = F)
heatmap( as.matrix(t3ss[which(t3ss$CHR_haplotype == 'K1_R3_I1_G1'),2:39]), col = c("grey", "salmon"), Rowv = F)




########################################################
# Discriminant analysis of principal components (DAPC) #
########################################################

# Find the optimum number of PC
prmx <- xvalDapc(tab(sym.snp, NA.method = "mean"), pop(chr.snp), n.pca = 1:30, n.rep = 100, parallel = "multicore", ncpus = 4L)

# DAPC and plot Plot 
dapc.chr.snp <- dapc(chr.snp, var.contrib = TRUE, scale = FALSE, n.pca = 28, n.da = nPop(chr.snp) - 1)
scatter(dapc.chr.snp, cell = 0, pch = 18:23, cstar = 0, mstree = TRUE, lwd = 2, lty = 2)

######################################################
# Selection analysis: prepare tree without bootstrap #
######################################################
library(ape)

tree <- read.tree("/Users/arafat/GDrive/Sachs/Chapter3_Epidemic_Genotype/Selection_Analysis/ete3_pipeline/SI_selection_analysis/tree/rxml.txt.treefile")
tree$node.label <- NULL
write.tree(tree, "/Users/arafat/GDrive/Sachs/Chapter3_Epidemic_Genotype/Selection_Analysis/ete3_pipeline/SI_selection_analysis/tree/coreSItree.nwk")

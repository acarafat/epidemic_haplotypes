#####################
# Epidemic Genotype #
#####################

#####################################
# Plotting final tree with outgroup #
#####################################

# Dec 19, 2022

library(ggtree)
library(ggplot2)
library(phangorn)
library(ggnewscale)
library(treeio)
library(viridis)

chr.tree <- read.iqtree('/Users/arafat/GDrive/Sachs/Chapter3_Epidemic_Genotype/Sequences/Bradyrhizobium/Dataset/Complete/Dataset/CHR_tree_with_outgroup/partition.raxml.treefile')
chr.tree@phylo <- midpoint(chr.tree@phylo)

# anotate tree figures to visualize the other meta information
meta <- read.table('/Users/arafat/GDrive/Sachs/Chapter3_Epidemic_Genotype/Sequences/Bradyrhizobium/Dataset/Complete/Dataset/CHR_tree_with_outgroup/meta_clade_ncbi.csv', sep=',', header=T)

# check if tips from phylogeny are absent/mismatch from the meta file
for (i in chr.tree@phylo$tip.label){
  if (!(i %in% meta$Strain)){
    print(i)
  }
}

# associated matrix
meta_matrix <- meta[,c('Strain', 'Site', 'Host', 'CHR.clade', 'Source')]
rownames(meta_matrix) <- meta_matrix$Strain
meta_matrix <- meta_matrix[,-1]
#meta_matrix[which(meta_matrix$Source != 'NCBI'),]$Source 
head(meta_matrix)
# - Show bootstrap value, need a way to incorporate sym.tree$node.label
# - Add Site+Host info

# Ornamentation of the tree
t1 <- ggtree(chr.tree) %<+% meta

#geom_tippoint(aes(color = Site)) #+
#geom_tiplab(aes(label=SYM_haplotype), size=1, check.overlap=T) 

# Plot tree with associated matrix
meta_matrix.site <- as.data.frame(meta_matrix[,'Site'])
meta_matrix.host <- as.data.frame(meta_matrix[,'Host'])
meta_matrix.clade <- as.data.frame(meta_matrix[,'CHR.clade'])
meta_matrix.source <- as.data.frame(meta_matrix[,'Source'])

colnames(meta_matrix.site) <- 'Sample\nSite'
colnames(meta_matrix.host) <- 'Host'
colnames(meta_matrix.clade) <- 'Clade'
colnames(meta_matrix.source) <- 'Source'

rownames(meta_matrix.site) <- rownames(meta_matrix)
rownames(meta_matrix.host) <- rownames(meta_matrix)
rownames(meta_matrix.clade) <- rownames(meta_matrix)
rownames(meta_matrix.source) <- rownames(meta_matrix)

# Long tree for saving into a PDF

t1 <- ggtree(chr.tree, size = 0.8) %<+% meta
t1$data$bootstrap <- '0'
t1$data[which(t1$data$SH_aLRT >= 70 & t1$data$UFboot  >= 70),]$bootstrap <- '1'

t1_tips <- t1 + geom_tree(aes(color=bootstrap == '1')) +
  scale_color_manual(name='Bootstrap', values=setNames( c('black', 'grey'), c(T,F)), guide = "none") +
  geom_tippoint(aes(fill=Source), shape=21) + 
  #geom_tiplab(align = T, linetype = "dotted", linesize = 0.1, size=2.5, hjust = -0.1) +
  new_scale_fill()


t1_mat1 <- gheatmap(t1_tips, meta_matrix.site, width=0.05, offset=0.015, colnames_angle=90, 
                    colnames_offset_y=4, colnames_position='top', custom_column_labels='Site') +
  scale_fill_viridis_d(option="A", name="Sample\nSite (population)") +
  new_scale_fill()

t1_mat1 <- gheatmap(t1_mat1, meta_matrix.clade, width=0.05, offset=0.025, colnames_angle=90,
         colnames_offset_y=10, colnames_position='top', custom_column_labels='Clade') +
  scale_fill_viridis_d(option="B", name="Clade") +
  scale_y_continuous(expand=c(0.03, 0.3)) + 
  new_scale_fill()

t1_mat1 <- gheatmap(t1_mat1, meta_matrix.host, width=0.05, offset=0.035,colnames_angle=90,
                    colnames_offset_y=10, colnames_position='top', custom_column_labels='Host') +
  scale_fill_viridis_d(option="C", name="Host") +
  scale_y_continuous(expand=c(0.03, 0.3))

t_mat1 <- t1_mat1 + theme(legend.key.size = unit(1,"line"))

# Plot tree for CHR 
pdf(file = "/Users/arafat/GDrive/Sachs/Chapter3_Epidemic_Genotype/Sequences/Bradyrhizobium/Dataset/Complete/Dataset/CHR_tree_with_outgroup/phylogeny_brady_epidemic.pdf",   # The directory you want to save the file in
    width = 15, # The width of the plot in inches
    height = 60) # The height of the plot in inches

t_mat1
dev.off()


# Tree without the tiplabs
t1 <- ggtree(chr.tree, size = 0.8) %<+% meta
t1$data$bootstrap <- '0'
t1$data[which(t1$data$SH_aLRT >= 70 & t1$data$UFboot  >= 70),]$bootstrap <- '1'

t1_tips <- t1 + geom_tree(aes(color=bootstrap == '1')) +
  scale_color_manual(name='Bootstrap', values=setNames( c('black', 'grey'), c(T,F)), guide = "none") +
  geom_tippoint(aes(fill=Source), shape=21) + 
  #geom_tiplab(align = T, linetype = "dotted", linesize = 0.1, size=2.5, hjust = -0.1) +
  new_scale_fill()


t1_mat1 <- gheatmap(t1_tips, meta_matrix.site, width=0.05, offset=0.015, colnames_angle=90, 
                    colnames_offset_y=4, colnames_position='top', custom_column_labels='Site') +
  scale_fill_viridis_d(option="A", name="Sample\nSite (population)") +
  new_scale_fill()

t1_mat1 <- gheatmap(t1_mat1, meta_matrix.clade, width=0.05, offset=0.025, colnames_angle=90,
                    colnames_offset_y=10, colnames_position='top', custom_column_labels='Clade') +
  scale_fill_viridis_d(option="B", name="Clade") +
  scale_y_continuous(expand=c(0.03, 0.3)) + 
  new_scale_fill()

t1_mat1 <- gheatmap(t1_mat1, meta_matrix.host, width=0.05, offset=0.035,colnames_angle=90,
                    colnames_offset_y=10, colnames_position='top', custom_column_labels='Host') +
  scale_fill_viridis_d(option="C", name="Host") +
  scale_y_continuous(expand=c(0.03, 0.3))

t_mat1 <- t1_mat1 + theme(legend.key.size = unit(1,"line"))

t_mat1

# Circular tree
t1 <- ggtree(chr.tree, layout="fan") %<+% meta
t1$data$bootstrap <- '0'
t1$data[which(t1$data$SH_aLRT >= 70 & t1$data$UFboot  >= 70),]$bootstrap <- '1'

t1_tips <- t1 + geom_tree(aes(color=bootstrap == '1')) +
  scale_color_manual(name='Bootstrap', values=setNames( c('black', 'grey'), c(T,F)), guide = "none") +
  geom_tippoint(aes(fill=Source), shape=21) + 
  new_scale_fill() + geom_treescale(x=0.32, y=310, offset=5)  


t1_mat1 <- gheatmap(t1_tips, meta_matrix.site, width=0.1, offset=0.0, colnames_angle=80, 
                    colnames_offset_y=10, colnames_position='top', custom_column_labels='Site') +
  scale_fill_viridis_d(option="A", name="Sample\nSite (population)") +
  new_scale_fill()

t1_mat1 <- gheatmap(t1_mat1, meta_matrix.clade, width=0.1, offset=0.02, colnames_angle=80,
                    colnames_offset_y=10, colnames_position='top', custom_column_labels='Clade') +
  scale_fill_viridis_d(option="H", name="Clade") +
  scale_y_continuous(expand=c(0.03, 0.3)) + 
  new_scale_fill()

t1_mat1 <- gheatmap(t1_mat1, meta_matrix.host, width=0.1, offset=0.04,colnames_angle=80,
                    colnames_offset_y=10, colnames_position='top', custom_column_labels='Host') +
  scale_fill_viridis_d(option="G", name="Host") +
  scale_y_continuous(expand=c(0.03, 0.3))

t_mat1 <- t1_mat1 + theme(legend.key.size = unit(1,"line"))
t_mat1 


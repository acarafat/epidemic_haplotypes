library(ggplot2)
library(ggtree)
library(phangorn)
library(dplyr)


# Meta file
meta <- read.table('~/GDrive/Sachs/Chapter3_Epidemic_Genotype/Sequences/Bradyrhizobium/Dataset/Complete/Dataset/meta_all.csv', sep=',', header = T)

# Concatenation of Chr loci
chr.tree <- read.tree('~/GDrive/Sachs/Chapter3_Epidemic_Genotype/Sequences/Bradyrhizobium/Dataset/Complete/Dataset/tanglegram/chr.meta.unique.treefile')
chr.tree <- midpoint(chr.tree)

# Concatenation of Sym loci
sym.tree <- read.tree('~/GDrive/Sachs/Chapter3_Epidemic_Genotype/Sequences/Bradyrhizobium/Dataset/Complete/Dataset/tanglegram/sym.meta.unique.treefile')
sym.tree <- midpoint(sym.tree)

# Another Meta for assisting the funtion
symHapList = unique(meta$SYM_haplotype)


# Function to make the tanglegrams
# example chrHaplotype = 'K1_R1_I1_G1'

makeTanglegram <- function(chrHaplotype, title){   
  meta2 <- meta
  # Need to update the meta file so that no dominant CHR missed
  for (symhap in symHapList){
    if (chrHaplotype %in% meta[which(meta$SYM_haplotype == symhap), 'CHR_haplotype']){
      meta[which(meta$SYM_haplotype == symhap),]$CHR_haplotype = chrHaplotype
      }
  }
  # Attach meta with tree dataframe
  t1 <-ggtree(chr.tree)  %<+%  meta2[,2:28]
  t2 <- ggtree(sym.tree)  %<+%  meta[,2:28] 
  
  d1 <- t1$data
  d2 <- t2$data
  d1$tree <-'t1'
  d2$tree <-'t2'
  
  d2$x <- max(d2$x) - d2$x + max(d1$x) + 0.05
  pp <- t1 + geom_tree(data=d2)
  
  dd <- bind_rows(d1, d2) %>% 
    filter(isTip == TRUE)
  
  
  dd1 <- as.data.frame(dd) 
  
  dd1$dom <- 0
  dd1[which(dd1$CHR_haplotype == chrHaplotype),]$dom <- chrHaplotype
  
  # Tangleram for Haplotype "K1_R1_I1_G1"
  hap1_tree <- dd1[which(dd1$CHR_haplotype == chrHaplotype), c('label', 'x', 'y', 'tree')]
  hap1_tree.repeat <- hap1_tree[2:dim(hap1_tree)[1], ]
  hap1_tree.repeat$x <- hap1_tree[which(hap1_tree$tree == 't1'), 'x']
  hap1_tree.repeat$y <- hap1_tree[which(hap1_tree$tree == 't1'), 'y']
  hap1_tree <- rbind(hap1_tree, hap1_tree.repeat)
  pp + geom_line(aes(x, y, group=label), data=hap1_tree, color='grey') + 
    ggtitle(title)
}

p1 <- makeTanglegram('K1_R1_I1_G1', 'A. K1_R1_I1_G1') # 70 lines?
p2 <- makeTanglegram('K1_R1_I3_G1', 'B. K1_R1_I3_G1') # 11 lines ok
p3 <- makeTanglegram('K1_R3_I1_G1', 'C. K1_R3_I1_G1') # 5 lines ok
p4 <- makeTanglegram('K1_R1_I6_G1', 'D. K1_R1_I6_G1') # 3 lines ok
p5 <- makeTanglegram('K18_R1_I1_G1', 'E. K18_R1_I1_G1') # 3 lines ok
p6 <- makeTanglegram('K21_R23_I55_G4', 'F. K21_R23_I55_G4') # 2 lines ok

grid.arrange(p1, p2, p3, p4, p5, p6, nrow=2)

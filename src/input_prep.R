###################################
# Preparing input files for Roary #
###################################

matrix = read.table('/Users/arafat/GDrive/Sachs/Chapter3_Epidemic_Genotype/GWAS_analysis/pangenome_matrix_genes_t0.tr.tab', header=T, sep='\t')

# Following isolates id changes
# 48-1
# 61-1
# inoc15-2
# inoc4-2

# Edit first column as Gene header
colnames(matrix)[1] = 'Gene'

# Remove the isolates with B_ initials
matrix = matrix[, !startsWith(colnames(matrix), 'B_')]

# Remove .gbk 
colnames(matrix) = gsub('.gbk', '', colnames(matrix))

# Remove last X column
matrix = subset(matrix, select= -c(X))

# Remove X from the beginning of strain id
colnames(matrix) = gsub('X', '', colnames(matrix))


# Give the matrix gene_presence_absence format
matrix$`Non-unique Gene Name` <- NA
matrix$Annotation <- NA
matrix$`No. isolates` <- NA
matrix$`No. sequences` <- NA
matrix$`Avg sequences per isolate` <- NA
matrix$`Genome Framgment` <- NA
matrix$`Order within Fragment` <- NA
matrix$`Accessory Fragment` <- NA
matrix$`Accessory Order with Fragment` <- NA
matrix$QC <- NA
matrix$`Min group size nuc` <- NA
matrix$`Max group size nuc` <- NA
matrix$`Avg group size nuc` <- NA


# Rearrange columns
matrix =  matrix[, c(1, 267:279, 2:266)]

# Update some values
matrix$`No. isolates` <- rowSums(matrix[, 15:179] != '-')
#matrix$`No. sequences` <- rowSums(matrix[, 15:179])
#matrix$`Avg sequences per isolate`<- matrix$`No. sequences` / matrix$`No. isolates`


# lapply(colnames(matrix[, 15:179]), function(x) gsub('-', '', matrix[[x]]))
test = data.frame(lapply(matrix, function(x) gsub('-', '', x)))

# Save the matrix and try SCOARY on it
write.csv(matrix, file='/Users/arafat/GDrive/Sachs/Chapter3_Epidemic_Genotype/Epidemic Genotype/GWAS_analysis/pan_gene_presence2.csv', sep=',')

# Matrix with only gene presence absence
matrix.trunc <- matrix[, c(1, 15:279)]
colnames(matrix.trunc)
write.csv(matrix.trunc, file='/Users/arafat/GDrive/Sachs/Chapter3_Epidemic_Genotype/Epidemic Genotype/GWAS_analysis/pan_gene_presence2.trunc.csv', sep=',', row.names=F)

#######################
# Make the trait file #
#######################

# Subset Meta file according to traits to be tested
meta <- read.table('/Users/arafat/GDrive/Sachs/Chapter3_Epidemic_Genotype/Epidemic Genotype/Sequences/Bradyrhizobium/Dataset/Complete/Dataset/meta_clade.csv')

# B. canariense species
meta.acs <- meta[which(meta$CHR.clade == "B. canariense"), 'Strain']

# Epidemic Haplotypes
meta.all1 <- meta[which(meta$CHR_haplotype == "K1_R1_I1_G1"), 'Strain']
meta.I3 <- meta[which(meta$CHR_haplotype == "K1_R1_I3_G1"), 'Strain']
meta.R3 <- meta[which(meta$CHR_haplotype == "K1_R3_I1_G1"), 'Strain']
meta.epi <- meta[which(meta$CHR_haplotype == "K1_R1_I1_G1" | meta$CHR_haplotype == "K1_R1_I3_G1" | 
                         meta$CHR_haplotype == "K1_R3_I1_G1"), 'Strain']

# Only dominant haplotypes
meta.I6 <- meta[which(meta$CHR_haplotype == "K1_R1_I6_G1"), ]
meta.K18 <- meta[which(meta$CHR_haplotype == "K18_R1_I1_G1"), ]
meta.K21 <- meta[which(meta$CHR_haplotype == "K21_R23_I55_G4"), ]
meta.dom <- meta[which(meta$CHR_haplotype == "K1_R1_I1_G1" | meta$CHR_haplotype == "K1_R1_I3_G1" | 
                         meta$CHR_haplotype == "K1_R3_I1_G1" | meta$CHR_haplotype == "K1_R1_I6_G1" |
                         meta$CHR_haplotype == "K18_R1_I1_G1" | meta$CHR_haplotype == "K21_R23_I55_G4"), 'Strain']

## Trait file
trait <- as.data.frame(colnames(matrix)[15:279])
colnames(trait) <- 'isolate'

# Initiate
trait$dominant <- 0 
trait$epidemic <- 0
trait$all1 <- 0
trait$I3 <- 0
trait$R3 <- 0
trait$spp <- 0

trait[which(trait$isolate %in% meta.dom), 'dominant'] <- 1
trait[which(trait$isolate %in% meta.epi), 'epidemic'] <- 1
trait[which(trait$isolate %in% meta.acs), 'spp'] <- 1

trait[which(trait$isolate %in% meta.all1), 'all1'] <- 1
trait[which(trait$isolate %in% meta.I3), 'I3'] <- 1
trait[which(trait$isolate %in% meta.R3), 'R3'] <- 1

write.table(trait, file='/Users/arafat/GDrive/Sachs/Chapter3_Epidemic_Genotype/Epidemic Genotype/GWAS_analysis/trait.csv', row.names=F, sep=',')




##############################
# Playing with Scoary output #
##############################
# - venn diagram
# - enrichment analysis

dominant <- read.table('~/GDrive/Sachs/Chapter3_Epidemic_Genotype/Epidemic Genotype/GWAS_analysis/dominant_17_11_2021_2222.results.csv', 
           sep=',', header=T)
epidemic <- read.table('~/GDrive/Sachs/Chapter3_Epidemic_Genotype/Epidemic Genotype/GWAS_analysis/epidemic_17_11_2021_2222.results.csv', 
                       sep=',', header=T)
spp <- read.table('~/GDrive/Sachs/Chapter3_Epidemic_Genotype/Epidemic Genotype/GWAS_analysis/spp_17_11_2021_2222.results.csv', 
                       sep=',', header=T)

all1 <- read.table('~/GDrive/Sachs/Chapter3_Epidemic_Genotype/Epidemic Genotype/GWAS_analysis/all1_17_11_2021_2222.results.csv', 
                       sep=',', header=T)

I3 <- read.table('~/GDrive/Sachs/Chapter3_Epidemic_Genotype/Epidemic Genotype/GWAS_analysis/I3_17_11_2021_2222.results.csv', 
                       sep=',', header=T)

R3 <- read.table('~/GDrive/Sachs/Chapter3_Epidemic_Genotype/Epidemic Genotype/GWAS_analysis/R3_17_11_2021_2222.results.csv', 
                       sep=',', header=T)


dominant.bonf <- dominant[which(dominant$Bonferroni_p <= 0.05), ]
epidemic.bonf <- epidemic[which(epidemic$Bonferroni_p <= 0.05), ]
spp.bonf <- spp[which(spp$Bonferroni_p <= 0.05), ]
all1.bonf <- all1[which(all1$Bonferroni_p <= 0.05), ]
I3.bonf <- I3[which(I3$Bonferroni_p <= 0.05), ]
R3.bonf <-R3[which(R3$Bonferroni_p <= 0.05), ]

dim(dominant.bonf)

library(ggvenn)

x <- list(
  Dominant = dominant.bonf$Gene,
  Epidemic = epidemic.bonf$Gene,
  B_canariense = spp.bonf$Gene,
  K1_R1_I1_G1 = all1.bonf$Gene
)

ggvenn(
  x, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 4
)

####################################
# Converting gene names to product #
####################################

hash = read.csv('~/GDrive/Sachs/Chapter3_Epidemic_Genotype/Epidemic Genotype/GWAS_analysis/hash_table.csv', header=F)
colnames(hash) = c('Gene', 'Product')

# based on Feb 1 analysis

dominant.contrast <- read.csv('/Users/arafat/GDrive/Sachs/Chapter3_Epidemic_Genotype/Epidemic Genotype/GWAS_analysis/Feb1_allTrait_allPanGenes/Results/contrast/contrast_dominant.csv', header=T)
epidemic.contrast <- read.csv('/Users/arafat/GDrive/Sachs/Chapter3_Epidemic_Genotype/Epidemic Genotype/GWAS_analysis/Feb1_allTrait_allPanGenes/Results/contrast/contrast_epidemic.csv', header=T)
spp.contrast <- read.csv('/Users/arafat/GDrive/Sachs/Chapter3_Epidemic_Genotype/Epidemic Genotype/GWAS_analysis/Feb1_allTrait_allPanGenes/Results/contrast/contrast_spp.csv', header=T)
all1.contrast <- read.csv('/Users/arafat/GDrive/Sachs/Chapter3_Epidemic_Genotype/Epidemic Genotype/GWAS_analysis/Feb1_allTrait_allPanGenes/Results/contrast/contrast_all1haplotype.csv', header=T)
I3.contrast <- read.csv('/Users/arafat/GDrive/Sachs/Chapter3_Epidemic_Genotype/Epidemic Genotype/GWAS_analysis/Feb1_allTrait_allPanGenes/Results/contrast/contrast_I3haplotype.csv', header=T)
R3.contrast <- read.csv('/Users/arafat/GDrive/Sachs/Chapter3_Epidemic_Genotype/Epidemic Genotype/GWAS_analysis/Feb1_allTrait_allPanGenes/Results/contrast/contrast_R3haplotype.csv', header=T)

dominant.product <- unique(hash[which(hash$Gene %in% dominant.contrast[which(dominant.contrast$Association == 'Positive'),]$Gene), ]$Product)
epidemic.product <- unique(hash[which(hash$Gene %in% epidemic.contrast[which(epidemic.contrast$Association == 'Positive'),]$Gene), ]$Product)
spp.product <- unique(hash[which(hash$Gene %in% spp.contrast[which(spp.contrast$Association == 'Positive'),]$Gene), ]$Product)
all1.product <- unique(hash[which(hash$Gene %in% all1.contrast[which(all1.contrast$Association == 'Positive'),]$Gene), ]$Product)
I3.product <- unique(hash[which(hash$Gene %in% I3.contrast[which(I3.contrast$Association == 'Positive'),]$Gene), ]$Product)
R3.product <- unique(hash[which(hash$Gene %in% R3.contrast[which(R3.contrast$Association == 'Positive'),]$Gene), ]$Product)

# Prepare input for PANTHER
paste(R3.product, sep=',', collapse = ',')

# Gene set enrichment in PANTHER
dominant.panther <- read.table('~/GDrive/Sachs/Chapter3_Epidemic_Genotype/Epidemic Genotype/GWAS_analysis/dominant.bonf.product_pantherGeneList.txt', sep='\t')
dominant.panther$V1



gsub('.fna', '', dominant.bonf$Gene)[499]

dominant.bonf$Gene[grepl('..', dominant.bonf$Gene, fixed=T)]

##########################################################################
# Converting gene names to product for pan-gene presence absense dataset #
##########################################################################
hash = read.csv('~/GDrive/Sachs/Chapter3_Epidemic_Genotype/Epidemic Genotype/GWAS_analysis/hash_table.csv', header=F)
colnames(hash) = c('Gene', 'Product')
head(hash)

pan_gene_pa = read.csv(file='/Users/arafat/GDrive/Sachs/Chapter3_Epidemic_Genotype/Epidemic Genotype/GWAS_analysis/pan_gene_presence.csv', sep=',')
dim(pan_gene_pa)
hash <- hash[order(match(hash$Gene, pan_gene_pa$Gene)), ]

##################################################### 
# Visualize Gene Set Enrichment Output from Panther #
##################################################### 
pan.data <- read.table('~/GDrive/Sachs/Chapter3_Epidemic_Genotype/Epidemic Genotype/GWAS_analysis/PANTHER_GO_Term_list.csv', sep=',', header=T)
head(pan.data)


library('GOsummaries')
dom = data.frame(Term= pan.data$dom_name[1:48], Score=0.9/pan.data$dom_count[1:48])
gs = gosummaries(wc_data = list(Results1 = dom))
plot(gs, filename = 'Dominant_haplotype.pdf')
dom


# Visualize panther over representation test output using customized code

pan.data <- read.table('~/GDrive/Sachs/Chapter3_Epidemic_Genotype/Epidemic Genotype/GWAS_analysis/nov17_scoary_analysis/panther_overrepresentation_test.csv', sep=',', header=T)
head(pan.data)

ggplot(pan.data, aes(y=reorder(spp.go.terms, spp.fold.enrich), x=spp.fold.enrich)) + 
  geom_point(aes(size=spp.gene.count, color=spp.p.adj)) + xlab('Fold Enrichment') + ylab('Biological processes') + ggtitle('B. canariense')

ggplot(pan.data, aes(y=reorder(all1.go.terms, all1.fold.enrichment), x=all1.fold.enrichment)) + 
  geom_point(aes(size=all1.gene.count, color=spp.p.adj)) + xlab('Fold Enrichment') + ylab('Biological processes') + ggtitle('K1_R1_I1_G1')

ggplot(pan.data, aes(y=reorder(epi.go.terms, epi.fold.enrichment), x=epi.fold.enrichment)) + 
  geom_point(aes(size=ep.gene.count, color=epi.p.adj)) + xlab('Fold Enrichment') + ylab('Biological processes') + ggtitle('Epidemic')

ggplot(pan.data, aes(y=reorder(dom.go.terms, dom.gene.enrichment), x=dom.gene.enrichment)) + 
  geom_point(aes(size=dom.genes.count, color=dom.p.adj)) + xlab('Fold Enrichment') + ylab('Biological processes') + ggtitle('Dominant')

# 2 Feb 2022
dom <- read.csv('~/GDrive/Sachs/Chapter3_Epidemic_Genotype/GWAS_analysis/Feb1_allTrait_allPanGenes/Results/contrast/GO/dominant.tsv', sep='\t', header = T)
epi <- read.csv('~/GDrive/Sachs/Chapter3_Epidemic_Genotype/GWAS_analysis/Feb1_allTrait_allPanGenes/Results/contrast/GO/epidemic.tsv', sep='\t', header = T)
spp <- dom <- read.csv('~/GDrive/Sachs/Chapter3_Epidemic_Genotype/GWAS_analysis/Feb1_allTrait_allPanGenes/Results/contrast/GO/species.tsv', sep='\t', header = T)
all1 <- read.csv('~/GDrive/Sachs/Chapter3_Epidemic_Genotype/GWAS_analysis/Feb1_allTrait_allPanGenes/Results/contrast/GO/all1.tsv', sep='\t', header = T)
I3 <- read.csv('~/GDrive/Sachs/Chapter3_Epidemic_Genotype/GWAS_analysis/Feb1_allTrait_allPanGenes/Results/contrast/GO/I3.tsv', sep='\t', header = T)
R3 <- read.csv('~/GDrive/Sachs/Chapter3_Epidemic_Genotype/GWAS_analysis/Feb1_allTrait_allPanGenes/Results/contrast/GO/R3.tsv', sep='\t', header = T)


p1 <- ggplot(dom, aes(y=reorder(GO_biological_processes, fold_enrichment), x=fold_enrichment)) + 
  geom_point(aes(size=X.Genes, color=p.value)) + xlab('Fold Enrichment') + ylab('Biological processes') + 
  ggtitle('Dominant') + theme(legend.position = "none")

p2 <- ggplot(epi, aes(y=reorder(GO_biological_processes, fold_enrichment), x=fold_enrichment)) + 
  geom_point(aes(size=X.Genes, color=p.value)) + xlab('Fold Enrichment') + ylab('Biological processes') + 
  ggtitle('Epidemic') + theme(legend.position = "none")

p3 <- ggplot(spp, aes(y=reorder(GO_biological_processes, fold_enrichment), x=fold_enrichment)) + 
  geom_point(aes(size=X.Genes, color=p.value)) + xlab('Fold Enrichment') + ylab('Biological processes') + 
  ggtitle('B. canariense') + theme(legend.position = "none")

p4 <- ggplot(all1, aes(y=reorder(GO_biological_processes, fold_enrichment), x=fold_enrichment)) + 
  geom_point(aes(size=X.Genes, color=p.value)) + xlab('Fold Enrichment') + ylab('Biological processes') + 
  ggtitle('K1_R1_I1_G1') + theme(legend.position = "none")

p5 <- ggplot(I3, aes(y=reorder(GO_biological_processes, fold_enrichment), x=fold_enrichment)) + 
  geom_point(aes(size=X.Genes, color=p.value)) + xlab('Fold Enrichment') + ylab('Biological processes') + 
  ggtitle('K1_R1_I3_G1') + theme(legend.position = "none")

p6 <- ggplot(R3, aes(y=reorder(GO_biological_processes, fold_enrichment), x=fold_enrichment)) + 
  geom_point(aes(size=X.Genes, color=p.value)) + xlab('Fold Enrichment') + ylab('Biological processes') + 
  ggtitle('K1_R3_I1_G1') 

grid.arrange(p1, p2, p3, p4, p5, p6, nrow=3)

#####################################################
# Prepare traits file for trait by trait comparison #
#####################################################
meta <- read.csv('/Users/arafat/GDrive/Sachs/Chapter3_Epidemic_Genotype/Epidemic Genotype/Sequences/Bradyrhizobium/Dataset/Complete/Dataset/meta_clade.csv')
meta$trait <- 0
trait <- read.csv('/Users/arafat/GDrive/Sachs/Chapter3_Epidemic_Genotype/Epidemic Genotype/GWAS_analysis/trait.csv', sep=',')
pan_gene_matrix <- read.table('/Users/arafat/GDrive/Sachs/Chapter3_Epidemic_Genotype/Epidemic Genotype/GWAS_analysis/pan_gene_presence2.csv', sep=',')

dominant <- c('K1_R1_I1_G1', 'K1_R1_I3_G1', 'K1_R3_I1_G1', 'K1_R1_I6_G1', 'K18_R1_I1_G1', 'K21_R23_I55_G14')
epidemic <-  c('K1_R1_I1_G1', 'K1_R1_I3_G1', 'K1_R3_I1_G1')

dom_vs_spp <- meta[which(meta$CHR.clade == 'B. canariense'), c('Strain', 'CHR_haplotype', 'trait')]
epi_vs_dom <- meta[which(meta$CHR.clade %in% dominant), c('Strain', 'CHR_haplotype', 'trait')]
all1_vs_r3 <- meta[which(meta$CHR_haplotype %in% c('K1_R1_I1_G1', 'K1_R1_I3_G1')), c('Strain', 'CHR_haplotype', 'trait')]
r3_vs_i3 <- meta[which(meta$CHR_haplotype %in% c('K1_R3_I1_G1', 'K1_R1_I3_G1')), c('Strain', 'CHR_haplotype', 'trait')]
i3_vs_all1 <- meta[which(meta$CHR_haplotype %in% c('K1_R1_I1_G1', 'K1_R1_I3_G1')), c('Strain', 'CHR_haplotype', 'trait')]



###################################################
# Scoary analysis with truncated dataset          #
# i.e. not all isolates needed in all comparison  #
###################################################
meta <- read.csv('/Users/arafat/GDrive/Sachs/Chapter3_Epidemic_Genotype/Epidemic Genotype/Sequences/Bradyrhizobium/Dataset/Complete/Dataset/meta_clade.csv')
trait <- read.csv('/Users/arafat/GDrive/Sachs/Chapter3_Epidemic_Genotype/Epidemic Genotype/GWAS_analysis/trait.csv', sep=',')

pan_gene_matrix <- read.table('/Users/arafat/GDrive/Sachs/Chapter3_Epidemic_Genotype/Epidemic Genotype/GWAS_analysis/pan_gene_presence2.csv', sep=',')

colnames(pan_gene_matrix) <- pan_gene_matrix[1,]
pan_gene_matrix <- pan_gene_matrix[-1,-1]




# Epidemic vs B canariense
dummy_columns <- colnames(pan_gene_matrix)[1:14]
spp_set <- meta[which(meta$CHR.clade == 'B. canariense'), ]$Strain



spp_compare_matrix <- pan_gene_matrix[,which( colnames(pan_gene_matrix) %in% dummy_columns |
                         colnames(pan_gene_matrix) %in% spp_set)]
spp_compare_trait <- trait[which(trait$X %in% spp_set),c('X', 'epidemic')]

write.csv(spp_compare_matrix, file='/Users/arafat/GDrive/Sachs/Chapter3_Epidemic_Genotype/Epidemic Genotype/GWAS_analysis/spp_compare_matrix.csv', row.names=F)
write.csv(spp_compare_trait, file='/Users/arafat/GDrive/Sachs/Chapter3_Epidemic_Genotype/Epidemic Genotype/GWAS_analysis/spp_compare_trait.csv', row.names=F)



# Read Scoary output
spp_compare <- read.table('~/GDrive/Sachs/Chapter3_Epidemic_Genotype/Epidemic Genotype/GWAS_analysis/epidemic_vs_spp/epidemic_01_12_2021_1046.results.csv', 
                       sep=',', header=T)
spp_compare.bonf <- spp_compare[which(spp_compare$Bonferroni_p <= 0.05), ]
dim(spp_compare.bonf)

hash = read.csv('~/GDrive/Sachs/Chapter3_Epidemic_Genotype/Epidemic Genotype/GWAS_analysis/hash_table.csv', header=F)
colnames(hash) = c('Gene', 'Product')

spp_compare.bonf.product <- unique(hash[which(hash$Gene %in% spp_compare.bonf$Gene), ]$Product)

paste(spp_compare.bonf.product, sep=',', collapse = ',')

# Panther visualization
pan.data <- read.table('~/GDrive/Sachs/Chapter3_Epidemic_Genotype/Epidemic Genotype/GWAS_analysis/epidemic_vs_spp_compare.csv', sep='\t', header=T)
head(pan.data)

ggplot(pan.data, aes(y=reorder(GO_terms, fold_enrichment), x=fold_enrichment)) + 
  geom_point(aes(size=query_genes, color=p_adj)) + xlab('Fold Enrichment') + ylab('Biological processes') + ggtitle('B. canariense vs Epidemic')


# Dominant vs Epidemic

dom_set <- meta[which(meta$CHR_haplotype == "K1_R1_I1_G1" | meta$CHR_haplotype == "K1_R1_I3_G1" | 
                        meta$CHR_haplotype == "K1_R3_I1_G1" | meta$CHR_haplotype == "K1_R1_I6_G1" |
                        meta$CHR_haplotype == "K18_R1_I1_G1" | meta$CHR_haplotype == "K21_R23_I55_G4"), 'Strain']


dom_compare_matrix <- pan_gene_matrix[,which( colnames(pan_gene_matrix) %in% dummy_columns |
                                                colnames(pan_gene_matrix) %in% dom_set)]
dom_compare_trait <- trait[which(trait$X %in% dom_set),c('X', 'epidemic')]

write.csv(dom_compare_matrix, file='/Users/arafat/GDrive/Sachs/Chapter3_Epidemic_Genotype/Epidemic Genotype/GWAS_analysis/dom_compare_matrix.csv', row.names=F)
write.csv(dom_compare_trait, file='/Users/arafat/GDrive/Sachs/Chapter3_Epidemic_Genotype/Epidemic Genotype/GWAS_analysis/dom_compare_trait.csv', row.names=F)

# Read Scoary output
dom_compare <- read.table('~/GDrive/Sachs/Chapter3_Epidemic_Genotype/Epidemic Genotype/GWAS_analysis/epidemic_vs_spp/epidemic_01_12_2021_1046.results.csv', 
                          sep=',', header=T)
dom_compare.bonf <- dom_compare[which(dom_compare$Bonferroni_p <= 0.05 & 
                                        dom_compare$Number_pos_present_in != 0 & 
                                        dom_compare$Number_pos_not_present_in < 10), ]
dim(dom_compare.bonf)

hash = read.csv('~/GDrive/Sachs/Chapter3_Epidemic_Genotype/Epidemic Genotype/GWAS_analysis/hash_table.csv', header=F)
colnames(hash) = c('Gene', 'Product')

dom_compare.bonf.product <- unique(hash[which(hash$Gene %in% dom_compare.bonf$Gene), ]$Product)

paste(dom_compare.bonf.product, sep=',', collapse = ',')

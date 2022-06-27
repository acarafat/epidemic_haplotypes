# Bootstrapped epiGeno FreqAnalysis
library("ggplot2")
library("gridExtra")
library('dplyr')
library('geosphere')
library('reshape2')
library('ggmap')
library('scatterpie')
library('ggrepel')
library('ggtree')
library('ggrepel')
library('phangorn')


###############
# Master meta #
###############

master <- read.csv('~/GDrive/Sachs/Chapter3_Epidemic_Genotype/Sequences/Bradyrhizobium/Dataset/Complete/Dataset/haplotype_meta.csv')
master.tab <- as.data.frame(table(master$Site, master$Host))
master.tab <- master.tab[which(master.tab$Freq != 0), ]
sum(master.tab[order(master.tab$Var1),]$Freq)

###################################
# Bootstrap without replacement   #
# One nodule from each plant only #
###################################

setwd('~/GDrive/Sachs/Chapter3_Epidemic_Genotype/Sequences/Bradyrhizobium/Dataset/Complete/Dataset/bootstrap/')

##########################################
# Incorporate haplotype_meta file as csv #
##########################################

temp <- list.files(pattern='*csv')
bootstrap_csv = lapply(temp, read.delim, sep=',')


# CHR concatenated haplotype
for (i in 1:1000){
  chr.haplotypes <- bootstrap_csv[[i]][, c('dnaK_coded', 'recA_coded', 'ITS_coded', 'glnII_coded')]
  bootstrap_csv[[i]]$CHR_haplotype <- apply(chr.haplotypes, 1, paste, collapse='_')
}

# SYM concatenated haplotype
for (i in 1:1000){
  sym.haplotypes <- bootstrap_csv[[i]][, c('nifD_coded', 'nodA_coded', 'nodL_coded', 'nodZ_coded')]
  bootstrap_csv[[i]]$SYM_haplotype <- apply(sym.haplotypes, 1, paste, collapse='_')
}

# Genome haplotype

for (i in 1:1000){
  genome.haplotypes <- bootstrap_csv[[i]][, c('CHR_haplotype', 'SYM_haplotype')]
  bootstrap_csv[[i]]$Genome_haplotype <- apply(genome.haplotypes, 1, paste, collapse='_')
}

# Initiate count variable
for (i in 1:1000){
  bootstrap_csv[[i]]$count <- 1
}


###################
# Define dominant #
###################

## Which haplotypes are dominant in each site?

# Need to filter for same plant isolation of strain
# Count per site >= 5, frequency per site >= 10%

chr_doms = list()
sym_doms = list()
genome_doms = list()


for (i in 1:1000){
  # CHR x Site
  count.table <- table(bootstrap_csv[[i]]$Site, bootstrap_csv[[i]]$CHR_haplotype)
  site.count <- as.data.frame(count.table)
  site.count <- site.count[which(site.count$Freq >= 4), ]
  site.freq <- as.data.frame(count.table/rowSums(count.table))
  site.freq <- site.freq[which(site.freq$Freq >= 0.1), ]
  chr.dom <- merge(site.count, site.freq, by=c('Var2', 'Var1'))
  colnames(chr.dom) <- c('Haplotype', 'Site', 'Count', 'Freq')
  chr_doms[[i]] <- chr.dom
  
  # SYM x Site
  count.table <- table(bootstrap_csv[[i]]$Site, bootstrap_csv[[i]]$SYM_haplotype)
  site.count <- as.data.frame(count.table)
  site.count <- site.count[which(site.count$Freq >= 4), ]
  site.freq <- as.data.frame(count.table/rowSums(count.table))
  site.freq <- site.freq[which(site.freq$Freq >= 0.1), ]
  sym.dom <- merge(site.count, site.freq, by=c('Var2', 'Var1'))
  colnames(sym.dom) <- c('Haplotype', 'Site', 'Count', 'Freq')
  sym_doms[[i]] <- sym.dom
  
  # Genome x Site
  count.table <- table(bootstrap_csv[[i]]$Site, bootstrap_csv[[i]]$Genome_haplotype)
  site.count <- as.data.frame(count.table)
  site.count <- site.count[which(site.count$Freq >= 4), ]
  site.freq <- as.data.frame(count.table/rowSums(count.table))
  site.freq <- site.freq[which(site.freq$Freq >= 0.1), ]
  genome.dom <- merge(site.count, site.freq, by=c('Var2', 'Var1'))
  colnames(genome.dom) <- c('Haplotype', 'Site', 'Count', 'Freq')
  genome_doms[[i]] <- genome.dom
}

#######################
# Plotting statistics #
#######################

# CHR
all_chr_doms <- do.call("rbind" ,chr_doms)
se <- function(x) sd(x)/sqrt(length(x))

all_chr_doms.summary <- aggregate(data= all_chr_doms[, c(1,3)], Count~Haplotype, mean)
all_chr_doms.se <- aggregate(data = all_chr_doms[, c(1,3)], Count~Haplotype, se)
colnames(all_chr_doms.se) <- c('Haplotype', 'SE')
colnames(all_chr_doms.summary) <- c('Haplotype', 'Mean')
all_chr_dom <- merge(all_chr_doms.se, all_chr_doms.summary, by='Haplotype')

# SYM
all_sym_doms <- do.call("rbind" ,sym_doms)

all_sym_doms.summary <- aggregate(data= all_sym_doms[, c(1,3)], Count~Haplotype, mean)
all_sym_doms.se <- aggregate(data = all_sym_doms[, c(1,3)], Count~Haplotype, se)
colnames(all_sym_doms.se) <- c('Haplotype', 'SE')
colnames(all_sym_doms.summary) <- c('Haplotype', 'Mean')
all_sym_dom <- merge(all_sym_doms.se, all_sym_doms.summary, by='Haplotype')

# Genome
all_genome_doms <- do.call("rbind" ,genome_doms)

all_genome_doms.summary <- aggregate(data= all_genome_doms[, c(1,3)], Count~Haplotype, mean)
all_genome_doms.se <- aggregate(data = all_genome_doms[, c(1,3)], Count~Haplotype, se)
colnames(all_genome_doms.se) <- c('Haplotype', 'SE')
colnames(all_genome_doms.summary) <- c('Haplotype', 'Mean')
all_genome_dom <- merge(all_genome_doms.se, all_genome_doms.summary, by='Haplotype')

# Plotting count
p1 <- ggplot(all_chr_dom, aes(x=reorder(Haplotype, -Mean), y=Mean)) +
  #geom_bar(stat='identity', fill='lightgray', color='black') +
  geom_point(size=1.5, color='blue') +
  geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE, width = 0.2)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1), 
        panel.border = element_blank(), 
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle('A. CHR') + ylab('Frequency') + xlab('') +
  scale_y_continuous(limits = c(3.5, 5)) 

p2 <- ggplot(all_sym_dom, aes(x=reorder(Haplotype, -Mean), y=Mean)) +
  #geom_bar(stat='identity',fill='lightgray', color='black') +
  geom_point(size=1.5, color='blue') +
  geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE, width = 0.2)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1), 
        panel.border = element_blank(), 
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle('B. SYM') + ylab('Frequency') + xlab('') +
  scale_y_continuous(limits = c(3.5, 5))

p3 <- ggplot(all_genome_dom, aes(x=reorder(Haplotype, -Mean), y=Mean)) +
  #geom_bar(stat='identity', fill='lightgray', color='black') +
  geom_point(size=1.5, color='blue') +
  geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE, width = 0.2)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1), 
        panel.border = element_blank(), 
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle('C. Genome') + ylab('Frequency') + xlab('') +
  scale_y_continuous(limits = c(3.5, 5))


grid.arrange(p1, p2, nrow=1)


####################################
# Define dominant without Acmispon #
####################################
## Which haplotypes are dominant in each site?

# Need to filter for same plant isolation of strain
# Count per site >= 5, frequency per site >= 10%

chr_doms = list()
sym_doms = list()
genome_doms = list()


for (i in 1:1000){
  # CHR x Site
  wo.acs <- bootstrap_csv[[i]][which(bootstrap_csv[[i]]$Host != 'A. strigosus'), ]
  count.table <- table(wo.acs$Site, wo.acs$CHR_haplotype)
  site.count <- as.data.frame(count.table)
  site.count <- site.count[which(site.count$Freq >= 3), ]
  site.freq <- as.data.frame(count.table/rowSums(count.table))
  site.freq <- site.freq[which(site.freq$Freq >= 0.1), ]
  chr.dom <- merge(site.count, site.freq, by=c('Var2', 'Var1'))
  colnames(chr.dom) <- c('Haplotype', 'Site', 'Count', 'Freq')
  chr_doms[[i]] <- chr.dom
  
  # SYM x Site
  count.table <- table(wo.acs$Site, wo.acs$SYM_haplotype)
  site.count <- as.data.frame(count.table)
  site.count <- site.count[which(site.count$Freq >= 3), ]
  site.freq <- as.data.frame(count.table/rowSums(count.table))
  site.freq <- site.freq[which(site.freq$Freq >= 0.1), ]
  sym.dom <- merge(site.count, site.freq, by=c('Var2', 'Var1'))
  colnames(sym.dom) <- c('Haplotype', 'Site', 'Count', 'Freq')
  sym_doms[[i]] <- sym.dom
  
  # Genome x Site
  count.table <- table(wo.acs$Site, wo.acs$Genome_haplotype)
  site.count <- as.data.frame(count.table)
  site.count <- site.count[which(site.count$Freq >= 3), ]
  site.freq <- as.data.frame(count.table/rowSums(count.table))
  site.freq <- site.freq[which(site.freq$Freq >= 0.1), ]
  genome.dom <- merge(site.count, site.freq, by=c('Var2', 'Var1'))
  colnames(genome.dom) <- c('Haplotype', 'Site', 'Count', 'Freq')
  genome_doms[[i]] <- genome.dom
}

## Plotting statistics

# CHR
all_chr_doms <- do.call("rbind" ,chr_doms)
se <- function(x) sd(x)/sqrt(length(x))

all_chr_doms.summary <- aggregate(data= all_chr_doms[, c(1,3)], Count~Haplotype, mean)
all_chr_doms.se <- aggregate(data = all_chr_doms[, c(1,3)], Count~Haplotype, se)
colnames(all_chr_doms.se) <- c('Haplotype', 'SE')
colnames(all_chr_doms.summary) <- c('Haplotype', 'Mean')
all_chr_dom <- merge(all_chr_doms.se, all_chr_doms.summary, by='Haplotype')

# SYM
all_sym_doms <- do.call("rbind" ,sym_doms)

all_sym_doms.summary <- aggregate(data= all_sym_doms[, c(1,3)], Count~Haplotype, mean)
all_sym_doms.se <- aggregate(data = all_sym_doms[, c(1,3)], Count~Haplotype, se)
colnames(all_sym_doms.se) <- c('Haplotype', 'SE')
colnames(all_sym_doms.summary) <- c('Haplotype', 'Mean')
all_sym_dom <- merge(all_sym_doms.se, all_sym_doms.summary, by='Haplotype')

# Genome
all_genome_doms <- do.call("rbind" ,genome_doms)

all_genome_doms.summary <- aggregate(data= all_genome_doms[, c(1,3)], Count~Haplotype, mean)
all_genome_doms.se <- aggregate(data = all_genome_doms[, c(1,3)], Count~Haplotype, se)
colnames(all_genome_doms.se) <- c('Haplotype', 'SE')
colnames(all_genome_doms.summary) <- c('Haplotype', 'Mean')
all_genome_dom <- merge(all_genome_doms.se, all_genome_doms.summary, by='Haplotype')

# Ploting count
p1 <- ggplot(all_chr_dom, aes(x=reorder(Haplotype, -Mean), y=Mean)) +
  #geom_bar(stat='identity', fill='lightgray', color='black') +
  geom_point(size=1.5, color='blue') +
  geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE, width = 0.2)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1), 
        panel.border = element_blank(), 
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle('CHR') + ylab('Frequency') + xlab('') +
  scale_y_continuous(limits = c(2.5, 5))

p2 <- ggplot(all_sym_dom, aes(x=reorder(Haplotype, -Mean), y=Mean)) +
  #geom_bar(stat='identity',fill='lightgray', color='black') +
  geom_point(size=1.5, color='blue') +
  geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE, width = 0.2)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1), 
        panel.border = element_blank(), 
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle('SYM') + ylab('Frequency') + xlab('') +
  scale_y_continuous(limits = c(2.5, 5))

p3 <- ggplot(all_genome_dom, aes(x=reorder(Haplotype, -Mean), y=Mean)) +
  #geom_bar(stat='identity', fill='lightgray', color='black') +
  geom_point(size=1.5, color='blue') +
  geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE, width = 0.2)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1), 
        panel.border = element_blank(), 
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle('Genome') + ylab('Frequency') + xlab('') +
  scale_y_continuous(limits = c(2.5, 5))


grid.arrange(p1, p2, nrow=1)



####################################
# Site host haplotype distribution #
####################################
hostSite <- as.data.frame(table(meta$Host, meta$CHR_haplotype, meta$Site))
hostSite <- hostSite[which(hostSite$Freq != 0), ]
colnames(hostSite) <- c('Host', 'CHR_haplotype', 'Site', 'Count')
sampleSite <- aggregate(Count~Site, hostSite, sum)

hostSite.multiple <- hostSite[which(hostSite$Count >= 1),]
hostSite.dcast <- dcast(hostSite.multiple, Site+Host ~ CHR_haplotype, value.var='Count')
hostSite.dcast[is.na(hostSite.dcast)] <- 0
hostSite.dcast$sum <- rowSums(hostSite.dcast[3:243])
hostSite.dcast[3:243] <- hostSite.dcast[3:243]/hostSite.dcast$sum

#################################
# Dominant Distribution by Site #
# Based on the Haplotype data   #
#################################
chr_dom.by_site <- aggregate(data= all_chr_doms[, c(1,2,3)], Count~Haplotype+Site, mean)
sym_dom.by_site <- aggregate(data= all_sym_doms[, c(1,2,3)], Count~Haplotype+Site, mean)
genome_dom.by_site <- aggregate(data= all_genome_doms[, c(1,2,3)], Count~Haplotype+Site, mean)

ggplot(chr_dom.by_site, aes(x=Haplotype, y=Count)) +
  geom_bar(stat='identity') + 
  facet_wrap(~Site) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1 )) +
  ggtitle('CHR Dominant Haplotypes by site')

ggplot(sym_dom.by_site, aes(x=Haplotype, y=Count)) +
  geom_bar(stat='identity') + 
  facet_wrap(~Site) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  ggtitle('SYM Dominant Haplotypes by site')

ggplot(genome_dom.by_site, aes(x=Haplotype, y=Count)) +
  geom_bar(stat='identity') + 
  facet_wrap(~Site) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  ggtitle('Genome Dominant Haplotypes by site')

###################################################
# Save a new meta file for concatenated haplotype #
##################################################
#meta.all <- read.table('~/GDrive/Sachs/Chapter3_Epidemic_Genotype/Sequences/Bradyrhizobium/Dataset/Complete/Dataset/haplotype_meta.csv', sep=',', header=T)

#chr.haplotypes <- meta.all[, c('dnaK_coded', 'recA_coded', 'ITS_coded', 'glnII_coded')]
#meta.all$CHR_haplotype <- apply(chr.haplotypes, 1, paste, collapse='_')

#sym.haplotypes <- meta.all[, c('nifD_coded', 'nodA_coded', 'nodL_coded', 'nodZ_coded')]
#meta.all$SYM_haplotype <- apply(sym.haplotypes, 1, paste, collapse='_')

#genome.haplotypes <- meta.all[, c('CHR_haplotype', 'SYM_haplotype')]
#meta.all$Genome_haplotype <- apply(genome.haplotypes, 1, paste, collapse='_')

# Now, need to find dominants in one site, are also present in other sites?
# Do it in Python??

#write.csv(meta.all, '~/GDrive/Sachs/Chapter3_Epidemic_Genotype/Sequences/Bradyrhizobium/Dataset/Complete/Dataset/meta_all.csv')
#write.csv(chr_dom.by_site, '~/Sachs/Epidemic Genotype/Sequences/Bradyrhizobium/Dataset/Complete/Dataset/chr_dom.csv')
#write.csv(sym_dom.by_site, '~/Sachs/Epidemic Genotype/Sequences/Bradyrhizobium/Dataset/Complete/Dataset/sym_dom.csv')


#############################################
# Haplotype Diversity 
# Now let's bring back the original dataset 
#############################################
meta = read.csv('~/GDrive/Sachs/Chapter3_Epidemic_Genotype/Sequences/Bradyrhizobium/Dataset/Complete/Dataset/meta_all.csv')

# CHR loci
# Get count of haplotype per site
hap_count.chr.table <- table(meta$Site, meta$CHR_haplotype)
hap_count.chr <- as.data.frame(hap_count.chr.table) 
hap_freq.chr <- as.data.frame(hap_count.chr.table/rowSums(hap_count.chr.table))

hap_count.chr <- hap_count.chr[which(hap_count.chr$Freq > 0), ]
#hap_freq.site <- hap_freq.site[which(hap_freq.site$Freq > 0), ]

hap_chr <- merge(hap_count.chr, hap_freq.chr, by=c('Var1', 'Var2'))
colnames(hap_chr) <- c('Site', 'Haplotype', 'Count', 'Freq')

# Stacked summary of diversity by site
ggplot(hap_chr, aes(x=Site, y=Count, fill=Haplotype)) +
  geom_bar(position='fill', stat='identity', show.legend = FALSE) + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) + 
  ggtitle('CHR Haplotype Diversity by Site')

# Table summarizing haplotype diversity in different sites


# SYM loci
# Get count of haplotype per site
hap_count.sym.table <- table(meta$Site, meta$SYM_haplotype)
hap_count.sym <- as.data.frame(hap_count.sym.table) 
hap_freq.sym <- as.data.frame(hap_count.sym.table/rowSums(hap_count.sym.table))

hap_count.sym <- hap_count.sym[which(hap_count.sym$Freq > 0), ]
hap_freq.sym <- hap_freq.sym[which(hap_freq.sym$Freq > 0), ]

hap_sym <- merge(hap_count.sym, hap_freq.sym, by=c('Var1', 'Var2'))
colnames(hap_sym) <- c('Site', 'Haplotype', 'Count', 'Freq')

# Stacked summary of diversity by site
ggplot(hap_sym, aes(x=Site, y=Count, fill=Haplotype)) +
  geom_bar(position='fill', stat='identity', show.legend = FALSE) + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) + 
  ggtitle('SYM Haplotype Diversity by Site')

# Table summarizing haplotype diversity in different sites
hap_chr$n <- 1
hap_summary <- aggregate(n~Site, data=hap_chr, FUN=sum)

hap_sym$n <- 1
hap_summary <- merge(hap_summary, aggregate(n~Site, data=hap_sym, FUN=sum), by='Site')
colnames(hap_summary) <- c('Site', 'CHR', 'SYM')
head(hap_summary)

hap_summary.melted <- melt(hap_summary, id.vars = 'Site')
colnames(hap_summary.melted) <- c('Site','Haplotype', 'Count')

ggplot(hap_summary.melted, aes(x=reorder(Site, -Count), y=Count, fill=Haplotype)) + 
  geom_bar(stat='identity') +
  theme(axis.text.x = element_text(angle=60, vjust=1, hjust=1)) +
  ggtitle('How many different haplotypes in each site?') +
  xlab('Site')

# Doing Chi-square test on haplotype distribution
hap_summary %>%
  rowwise() %>%
  mutate(
    test_stat = chisq.test(c(CHR, SYM))$statistic,
    p_val = chisq.test(c(CHR, SYM))$p.value
  ) %>%
  filter(p_val <= 0.05)

# Genome
# Get count of haplotype per site
hap_count.genome.table <- table(meta$Site, meta$Genome_haplotype)
hap_count.genome <- as.data.frame(hap_count.genome.table) 
hap_freq.genome <- as.data.frame(hap_count.genome.table/rowSums(hap_count.genome.table))

hap_count.genome <- hap_count.genome[which(hap_count.genome$Freq > 0), ]
hap_freq.genome <- hap_freq.genome[which(hap_freq.genome$Freq > 0), ]

hap_genome <- merge(hap_count.genome, hap_freq.genome, by=c('Var1', 'Var2'))
colnames(hap_genome) <- c('Site', 'Haplotype', 'Count', 'Freq')

# Stacked summary of diversity by site
ggplot(hap_genome, aes(x=Site, y=Count, fill=Haplotype)) +
  geom_bar(position='fill', stat='identity', show.legend = FALSE) + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) + 
  ggtitle('Genome Haplotype Diversity by Site')


###########################################
# Local distribution of dominant haplotype
###########################################
meta = read.csv('~/GDrive/Sachs/Chapter3_Epidemic_Genotype/Sequences/Bradyrhizobium/Dataset/Complete/Dataset/meta_all.csv')

dom_chr.local_dis <-hap_chr[which(hap_chr$Haplotype %in% chr_dom.by_site$Haplotype), ]

ggplot(dom_chr.local_dis, aes(x=reorder(Site, -Freq), y=Freq)) +
  geom_bar(stat='identity') + 
  theme(axis.text.x= element_text(angle=45, vjust=1, hjust=1)) +
  xlab('Sites') +
  facet_wrap(~Haplotype) +
  ggtitle('Local Distribution of dominant haplotypes')

dom_sym.local_dis <-hap_sym[which(hap_sym$Haplotype %in% sym_dom.by_site$Haplotype), ]

ggplot(dom_sym.local_dis, aes(x=reorder(Site, -Freq), y=Freq)) +
  geom_bar(stat='identity') + 
  theme(axis.text.x= element_text(angle=45, vjust=1, hjust=1)) +
  xlab('Sites') +
  facet_wrap(~Haplotype) +
  ggtitle('Local Distribution of dominant haplotypes')

dom_genome.local_dis <-hap_genome[which(hap_genome$Haplotype %in% genome_dom.by_site$Haplotype), ]

ggplot(dom_genome.local_dis, aes(x=reorder(Site, -Freq), y=Freq)) +
  geom_bar(stat='identity') + 
  theme(axis.text.x= element_text(angle=45, vjust=1, hjust=1)) +
  xlab('Sites') +
  facet_wrap(~Haplotype) +
  ggtitle('Local Distribution of dominant Genome haplotypes')

#############################################
# Haplotype Diversity by host and site
# Now let's bring back the original data set 
#############################################

# CHR loci
# Get count of haplotype per site per host
hap_count.chr.host.table <- table(meta$Site, meta$Host, meta$CHR_haplotype)
hap_count.chr.host <- as.data.frame(hap_count.chr.host.table) 
hap_freq.chr.host <- as.data.frame(hap_count.chr.host.table/rowSums(hap_count.chr.host.table))

hap_count.chr.host <- hap_count.chr.host[which(hap_count.chr.host$Freq > 0), ]
hap_freq.chr.host <- hap_freq.chr.host[which(hap_freq.chr.host$Freq > 0), ]

hap_chr.host <- merge(hap_count.chr.host, hap_freq.chr.host, by=c('Var1', 'Var2', 'Var3'))
colnames(hap_chr.host) <- c('Site', 'Host', 'Haplotype', 'Count', 'Freq')

dom_chr.site.host <- hap_chr.host[which(hap_chr.host$Haplotype %in% unique(chr_dom.by_site$Haplotype)), ]
dom_chr.site.host$Haplotype <- factor(dom_chr.site.host$Haplotype, levels=c('K1_R1_I1_G1', 'K1_R1_I3_G1', 'K1_R3_I1_G1',
                                                                  'K21_R23_I55_G4', 'K18_R1_I1_G1', 'K1_R1_I6_G1'))

# Stacked summary of diversity by site
ggplot(dom_chr.site.host, aes(x=reorder(Site, -Freq), y=Freq, fill=Host)) +
  geom_bar(stat='identity', show.legend = T, color='black') + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) + 
  facet_wrap(~Haplotype) +
  ggtitle('Dominant CHR haplotype') +
  xlab('Sampling Site') + ylab('Frequency / Sampling site')
  
# SYM loci
# Get count of haplotype per site per host
hap_count.sym.host.table <- table(meta$Site, meta$Host, meta$SYM_haplotype)
hap_count.sym.host <- as.data.frame(hap_count.sym.host.table) 
hap_freq.sym.host <- as.data.frame(hap_count.sym.host.table/rowSums(hap_count.sym.host.table))
  
hap_count.sym.host <- hap_count.sym.host[which(hap_count.sym.host$Freq > 0), ]
hap_freq.sym.host <- hap_freq.sym.host[which(hap_freq.sym.host$Freq > 0), ]
  
hap_sym.host <- merge(hap_count.sym.host, hap_freq.sym.host, by=c('Var1', 'Var2', 'Var3'))
colnames(hap_sym.host) <- c('Site', 'Host', 'Haplotype', 'Count', 'Freq')
  
dom_sym.site.host <- hap_sym.host[which(hap_sym.host$Haplotype %in% unique(sym_dom.by_site$Haplotype)), ]
  
# Stacked summary of diversity by site
ggplot(dom_sym.site.host, aes(x=reorder(Site, -Freq), y=Freq, fill=Host)) +
  geom_bar(stat='identity', show.legend = T) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) + 
  facet_wrap(~Haplotype) +
  ggtitle('Dominant SYM haplotypes') +
  xlab('Sampling Site') + ylab('Frequency / site')

# Genome loci
# Get count of haplotype per site per host
hap_count.genome.host.table <- table(meta$Site, meta$Host, meta$Genome_haplotype)
hap_count.genome.host <- as.data.frame(hap_count.genome.host.table) 
hap_freq.genome.host <- as.data.frame(hap_count.genome.host.table/rowSums(hap_count.genome.host.table))
  
hap_count.genome.host <- hap_count.genome.host[which(hap_count.genome.host$Freq > 0), ]
hap_freq.genome.host <- hap_freq.genome.host[which(hap_freq.genome.host$Freq > 0), ]
  
hap_genome.host <- merge(hap_count.genome.host, hap_freq.genome.host, by=c('Var1', 'Var2', 'Var3'))
colnames(hap_genome.host) <- c('Site', 'Host', 'Haplotype', 'Count', 'Freq')
  
dom_genome.site.host <- hap_genome.host[which(hap_genome.host$Haplotype %in% unique(genome_dom.by_site$Haplotype)), ]
  
# Stacked summary of diversity by site
ggplot(dom_genome.site.host, aes(x=reorder(Site, -Freq), y=Freq, fill=Host)) +
  geom_bar(stat='identity', show.legend = T) + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) + 
  facet_wrap(~Haplotype) +
  ggtitle('Dominant Genome haplotypes') +
  xlab('Site') + ylab('Frequency / site')

######################
# Epidemic Analysis  #
######################

# Subset of meta dataset in which haplotypes belong to dominant group
meta.epidemic.chr <- meta[which(meta$CHR_haplotype %in% unique(chr_dom.by_site$Haplotype)), ]
meta.epidemic.sym <- meta[which(meta$SYM_haplotype %in% unique(sym_dom.by_site$Haplotype)), ]
meta.epidemic.genome <- meta[which(meta$Genome_haplotype %in% unique(genome_dom.by_site$Haplotype)), ]

# Get mean coordinates for each dominant site
dom.chr.mean.coord <- meta.epidemic.chr %>% 
  group_by(CHR_haplotype, Site) %>% 
  summarize_each(funs(mean)) %>% 
  select(CHR_haplotype, Site, Longitude, Latitude)  

dom.sym.mean.coord <- meta.epidemic.sym %>% 
  group_by(SYM_haplotype, Site) %>% 
  summarize_each(funs(mean)) %>% 
  select(SYM_haplotype, Site, Longitude, Latitude)  

dom.chr.mean.coord <- as.data.frame(dom.chr.mean.coord)
dom.sym.mean.coord <- as.data.frame(dom.sym.mean.coord)
  
dom.chr.mean.coord <- dom.chr.mean.coord[c(10, 13, 18, 19, 23, 24),]
dom.sym.mean.coord <- dom.sym.mean.coord[c(1, 3, 5, 6, 7),]

# Get subset of dominant genotype found in other sites as well
epi.chr.sites <- meta.all[which(meta.all$CHR_haplotype %in% dom.chr.mean.coord$CHR_haplotype),]
epi.sym.sites <- meta.all[which(meta.all$SYM_haplotype %in% dom.sym.mean.coord$SYM_haplotype),]

# Function to calculate Distance
calcDist <- function(row, domLong, domLat) {distm(c(domLong, domLat), 
                                                  c(as.numeric(row[[2]]), as.numeric(row[[3]])), 
                                                  fun=distHaversine)*0.001}

# Calculate distance CHR
epi.chr.sites$Dist <- NA
epi.chr.sites$domHapLong <- NA
epi.chr.sites$domHapLat <- NA

for (i in 1:dim(epi.chr.sites)[1]){
  row = epi.chr.sites[i,c(1,21,20,19,25)]
  domHapSite <- dom.chr.mean.coord[which(dom.chr.mean.coord$CHR_haplotype == row[,5]), 2][1]
  domLong = dom.chr.mean.coord[which(dom.chr.mean.coord$CHR_haplotype == row[,5]), 3][1]
  domLat = dom.chr.mean.coord[which(dom.chr.mean.coord$CHR_haplotype == row[,5]), 4][1]
  
  dist <- as.numeric(unlist(calcDist(row, domLong, domLat)))
  
  epi.chr.sites[i,]$domHapLong <- domLong
  epi.chr.sites[i,]$domHapLat <- domLat

  epi.chr.sites[i,]$domHapSite <- domHapSite
  epi.chr.sites[i,]$Dist <- dist
}

# Calculate distance
epi.sym.sites$Dist <- NA
epi.sym.sites$domHapLong <- NA
epi.sym.sites$domHapLat <- NA

for (i in 1:dim(epi.sym.sites)[1]){
  row = epi.sym.sites[i,c(1,21,20,19,26)]
  domHapSite <- dom.sym.mean.coord[which(dom.sym.mean.coord$SYM_haplotype == row[,5]), 2]
  domLong = dom.sym.mean.coord[which(dom.sym.mean.coord$SYM_haplotype == row[,5]), 3]
  domLat = dom.sym.mean.coord[which(dom.sym.mean.coord$SYM_haplotype == row[,5]), 4]
  
  dist <- as.numeric(unlist(calcDist(row, domLong, domLat)))
  
  epi.sym.sites[i,]$domHapLong <- domLong
  epi.sym.sites[i,]$domHapLat <- domLat
  
  epi.sym.sites[i,]$domHapSite <- domHapSite
  epi.sym.sites[i,]$Dist <- dist
}

# Filter epiSpread for distance, unique location, unique
epi.chr.sites <- epi.chr.sites[which(epi.chr.sites$Dist >= 10),]
epi.sym.sites <- epi.sym.sites[which(epi.sym.sites$Dist >= 10),]

################################
# Haplotype count distribution #
########3#######################
meta = read.csv('~/GDrive/Sachs/Chapter3_Epidemic_Genotype/Sequences/Bradyrhizobium/Dataset/Complete/Dataset/meta_all.csv')

haploHostCount <-  as.data.frame.matrix(table(meta$CHR_haplotype, meta$Host))
haploHostArray <- rowSums(haploHostCount != 0)
hist(haploHostArray, xlab='# Host Species', ylab='# Haplotypes', main='Some haplotype infects multiple hosts')

haploHostDF <- as.data.frame(haploHostArray)

p0 <- ggplot(haploHostDF, aes(x=haploHostArray, y=..density..)) + 
  geom_histogram(binwidth = 1, fill = "blue", colour=1, alpha=0.5) +
  geom_density(aes(y=..scaled..), linetypgeogeo) +
  xlab('#Host Species') + ylab('Haplotype Density') +
  theme_bw() +
  ggtitle('A. Host range distribution of haplotypes')

######################################################################################
# Does dominant haplotype associated with more species compared to non-dominant one? #
######################################################################################
library(ggsignif)
library(ggpubr)

meta = read.csv('~/GDrive/Sachs/Chapter3_Epidemic_Genotype/Sequences/Bradyrhizobium/Dataset/Complete/Dataset/meta_all.csv')

# Subset meta for sites where multiple species has been sampled
meta.multiple <- as.data.frame.matrix(table(meta$Site, meta$Host))
meta.multiple$hostCount <- rowSums(meta.multiple != 0)
meta.multiple$sites <- rownames(meta.multiple)
meta.multiple.host <- meta.multiple[which(meta.multiple$hostCount > 2), 'sites']

chr.dom.haplotypes <- c('K1_R1_I1_G1', 'K1_R1_I3_G1', 'K1_R3_I1_G1',
                        'K21_R23_I55_G4', 'K18_R1_I1_G1', 'K1_R1_I6_G1')
#
meta.sub <- meta[which(meta$Site %in% meta.multiple.host),]
chr.host.stat <- as.data.frame.matrix(table(meta.sub$CHR_haplotype, meta.sub$Host))
chr.host.stat$totalHost <- rowSums(chr.host.stat != 0)
chr.host.stat$haplotype <- rownames(chr.host.stat)
chr.host.stat$chrDom <- 'Non-dominant'
chr.host.stat[which(chr.host.stat$haplotype %in% chr.dom.haplotypes),]$chrDom <- 'Dominant'


t.test(chr.host.stat[which(chr.host.stat$chrDom == 'Non-dominant'), ]$totalHost, 
       chr.host.stat[which(chr.host.stat$chrDom == 'Dominant'), ]$totalHost,
       alternative='two.sided', var.equal=T)

chr.hostAffinity <- aggregate(totalHost ~ chrDom, data=chr.host.stat, mean)
colnames(chr.hostAffinity) <- c('dominant', 'mean')
chr.hostAffinity <- merge(chr.hostAffinity, aggregate(totalHost ~ chrDom, data=chr.host.stat, se), by.x='dominant', by.y='chrDom')
colnames(chr.hostAffinity) <- c('dominant', 'mean', 'se')

p1 <- ggplot(chr.hostAffinity, aes(dominant, mean)) + geom_bar(stat='identity', colour='black', fill='blue', alpha=0.5) + 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2) +
  ggtitle('B. Host range by dominance') + theme_bw() + ylab('Mean Host Range') +
  annotate('text', label="paste(italic(P), \" = \", 3.007^{-06})", x=2, y=3, parse=T) +
  xlab('Haplotype Type')

p2 <- ggplot(chr.hostAffinity, aes(dominant, mean)) + geom_bar(stat='identity', fill='blue', alpha=0.5) + 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2) +
  ggtitle('Three or more host sampled') + theme_classic() +
  xlab('Haplotype')

#grid.arrange(p1, p2, ncol=2)

grid.arrange(p0, p1, nrow=1, widths=c(2,1))



###############################
# Plot epidemic spread on map #
###############################

# Dominant haplotype abundance dataframe CHR
hap_abund.chr <- hap_freq.chr[which(hap_freq.chr$Var2 %in% unique(chr_dom.by_site$Haplotype)),]
hap_abund.chr <- dcast(hap_abund.chr, Var1~Var2)
hap_abund.chr$other <- 1 - rowSums(hap_abund.chr[,2:7])
#hap_abund.chr$Longitude <- c(-116.4194,-116.4447, -116.4751, -121.5420, -123.063,-118.3089, -120.6088, -121.5631, -122.4090, -117.2575, -120.6224, -117.7080, -117.7800, -117.7691, -116.4548, -119.6897, -120.0493, -117.323, -116.6547)
#hap_abund.chr$Latitude <- c(33.2713,33.2713, 34.15315, 36.05631, 38.31930,34.12196, 35.00812, 36.37866, 38.85985, 33.80517, 35.05448, 34.11039, 34.13436, 34.16382, 33.58357, 34.01139, 34.69150, 33.96593, 33.97972)
#hap_abund.chr <- hap_abund.chr[which(hap_abund.chr$other != 1),]
#hap_abund.chr <- merge(hap_abund.chr, dom.chr.mean.coord[,c(2,3,4)], by.x='Var1', by.y='Site')

# Dominant haplotype abundance dataframe SYM
hap_abund.sym <- hap_freq.sym[which(hap_freq.sym$Var2 %in% unique(sym_dom.by_site$Haplotype)),]
hap_abund.sym <- dcast(hap_abund.sym, Var1~Var2)
hap_abund.sym[is.na(hap_abund.sym)] <- 0
hap_abund.sym$other <- 1 - rowSums(hap_abund.sym[,2:6])
#hap_abund.chr$Longitude <- c(-116.4194,-116.4447, -116.4751, -121.5420, -123.063,-118.3089, -120.6088, -121.5631, -122.4090, -117.2575, -120.6224, -117.7080, -117.7800, -117.7691, -116.4548, -119.6897, -120.0493, -117.323, -116.6547)
#hap_abund.chr$Latitude <- c(33.2713,33.2713, 34.15315, 36.05631, 38.31930,34.12196, 35.00812, 36.37866, 38.85985, 33.80517, 35.05448, 34.11039, 34.13436, 34.16382, 33.58357, 34.01139, 34.69150, 33.96593, 33.97972)
#hap_abund.chr <- hap_abund.chr[which(hap_abund.chr$other != 1),]

# Pie-chart
hap_abund.chr.melt <- melt(hap_abund.chr, id.vars = 'Var1')
ggplot(hap_abund.chr.melt, aes(x='', y=value, fill=variable)) +
   geom_bar(stat='identity', width=1, color='black') +
   coord_polar("y", start=0) +
   facet_wrap(~Var1) +
   theme_void() +
   scale_color_manual(values=c("#e28743", "#E69F00", "#21130d", '#00fffd', '#ff00eb', "#56B4E9", "#999999"))

# SYM Pie-chart
hap_abund.sym.melt <- melt(hap_abund.sym, id.vars = 'Var1')
ggplot(hap_abund.sym.melt, aes(x='', y=value, fill=variable)) +
  geom_bar(stat='identity', width=1, color='black') +
  coord_polar("y", start=0) +
  facet_wrap(~Var1) +
  theme_void()

# Generate map of CA
# Used bboxfinder to get CA bbox coordinates
# http://bboxfinder.com
ca_bbox = c(left = -124.782715, bottom = 32.472695, right = -114.147949, top = 40.032974)
map <- get_map(ca_bbox, zoom = 7, maptype = "terrain")


m1 <- ggmap(map) + 
  geom_curve(data=epi.chr.sites, 
             aes(x=domHapLong, y=domHapLat, xend=Longitude, yend=Latitude, color=CHR_haplotype), 
               arrow=arrow(length= unit(0.03, "npc"))) +
  coord_cartesian()

m1 + geom_label_repel(data=hap_abund.chr, aes(x=Longitude, y=Latitude, label=Var1)) + 
  geom_scatterpie(aes(x=Longitude, y=Latitude, r=0.2), 
                           cols=colnames(hap_abund.chr)[2:8], 
                           data=hap_abund.chr) + coord_equal() 


m2 <- ggmap(map) + 
  geom_curve(data=epi.sym.sites, 
             aes(x=domHapLong, y=domHapLat, xend=Longitude, yend=Latitude, color=SYM_haplotype), 
             arrow=arrow(length= unit(0.03, "npc"))) +
  coord_cartesian()

m2
##################################################################################################

#############################
# Bootstrap with replacement #
##############################
setwd('~/Sachs/Epidemic Genotype/Sequences/Bradyrhizobium/Dataset/Complete/Dataset/bootstrap_with_replacement')

##########################################
# Incorporate haplotype_meta file as csv #
##########################################

temp <- list.files(pattern='*csv')
bootstrap_csv = lapply(temp, read.delim, sep=',')


# CHR concatenated haplotype
for (i in 1:1000){
  chr.haplotypes <- bootstrap_csv[[i]][, c('dnaK_coded', 'recA_coded', 'ITS_coded', 'glnII_coded')]
  bootstrap_csv[[i]]$CHR_haplotype <- apply(chr.haplotypes, 1, paste, collapse='_')
}

# SYM concatenated haplotype
for (i in 1:1000){
  sym.haplotypes <- bootstrap_csv[[i]][, c('nifD_coded', 'nodA_coded', 'nodL_coded', 'nodZ_coded')]
  bootstrap_csv[[i]]$SYM_haplotype <- apply(sym.haplotypes, 1, paste, collapse='_')
}

# Initiate count variable
for (i in 1:1000){
  bootstrap_csv[[i]]$count <- 1
}


###################
# Define dominant #
###################

## Which haplotypes are dominant in each site?

# Need to filter for same plant isolation of strain
# Count per site >= 5, frequency per site >= 10%

chr_doms = list()
sym_doms = list()

for (i in 1:1000){
  # CHR x Site
  count.table <- table(bootstrap_csv[[i]]$Site, bootstrap_csv[[i]]$CHR_haplotype)
  site.count <- as.data.frame(count.table)
  site.count <- site.count[which(site.count$Freq >= 10), ]
  site.freq <- as.data.frame(count.table/rowSums(count.table))
  site.freq <- site.freq[which(site.freq$Freq >= 0.1), ]
  chr.dom <- merge(site.count, site.freq, by=c('Var2', 'Var1'))
  colnames(chr.dom) <- c('Haplotype', 'Site', 'Count', 'Freq')
  chr_doms[[i]] <- chr.dom
  
  # SYM x Site
  count.table <- table(bootstrap_csv[[i]]$Site, bootstrap_csv[[i]]$SYM_haplotype)
  site.count <- as.data.frame(count.table)
  site.count <- site.count[which(site.count$Freq >= 10), ]
  site.freq <- as.data.frame(count.table/rowSums(count.table))
  site.freq <- site.freq[which(site.freq$Freq >= 0.1), ]
  sym.dom <- merge(site.count, site.freq, by=c('Var2', 'Var1'))
  colnames(sym.dom) <- c('Haplotype', 'Site', 'Count', 'Freq')
  sym_doms[[i]] <- sym.dom
}


# Plotting statistics
all_chr_doms <- do.call("rbind" ,chr_doms)
se <- function(x) sd(x)/sqrt(length(x))

all_chr_doms.summary <- aggregate(data= all_chr_doms[, c(1,3)], Count~Haplotype, mean)
all_chr_doms.se <- aggregate(data = all_chr_doms[, c(1,3)], Count~Haplotype, se)
colnames(all_chr_doms.se) <- c('Haplotype', 'SE')
colnames(all_chr_doms.summary) <- c('Haplotype', 'Mean')
all_chr_dom <- merge(all_chr_doms.se, all_chr_doms.summary, by='Haplotype')


# Ploting count
ggplot(all_chr_dom, aes(x=Haplotype, y=Mean)) +
  geom_bar(stat='identity') +
  geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ggtitle('Dominance threshold count = 10')



# Dominant groups

# CHR
list_K1R1I1G1 = meta[which(meta$CHR_haplotype == 'K1_R1_I1_G1'), c('Strain')]
list_I3 = meta[which(meta$CHR_haplotype == 'K1_R1_I3_G1'), c('Strain')]
list_R1 = meta[which(meta$CHR_haplotype == 'K1_R3_I1_G1'), c('Strain')]

list_I6 = meta[which(meta$CHR_haplotype == 'K1_R1_I6_G1'), c('Strain')]
list_K18 = meta[which(meta$CHR_haplotype == 'K18_R1_I1_G1'), c('Strain')]
list_K21 = meta[which(meta$CHR_haplotype == 'K21_R23_I55_G4'), c('Strain')]



##########################
# Subset genome for GWAS #
##########################
meta = read.csv('~/GDrive/Sachs/Chapter3_Epidemic_Genotype/Sequences/Bradyrhizobium/Dataset/Complete/Dataset/meta_all.csv')
genomes = meta[which(meta$Source == 'Genomes'),]
head(genomes)

# CHR Haplotypes per sites
d1 <- colSums(table(genomes$CHR_haplotype, genomes$Site) != 0)

# SYM Haplotypes per sites
d2 <- colSums(table(genomes$SYM_haplotype, genomes$Site) != 0)

# Hosts per sites
d3 <- colSums(table(genomes$Host, genomes$Site) != 0)

par(mfrow=c(3,1))
barplot(d1, main='#Unique CHR haplotype per sampling site')
barplot(d2, main='#Unique SYM haplotype per sampling site')
barplot(d3, main='#Host species per sampling site')

write.csv(genomes, file='~/GDrive/Sachs/Chapter3_Epidemic_Genotype/Sequences/Bradyrhizobium/Dataset/Complete/Dataset/meta_genomes.csv', row.names = F)

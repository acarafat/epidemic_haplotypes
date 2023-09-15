#  Epidemic spread of *Bradyrhizobium* haplotypes


## Manuscript
**Epidemic Bradyrhizobium haplotypes adapt to host metapopulations via acquisition of diverse, host-specific symbiosis ICEs.**  Arafat Rahman, Lorena Torres Martínez, Alexandra J. Weisberg, Jeff H. Chang, Joel L. Sachs. *Under Review*


## Directory Structure
`data` folder contains dataset to do different evolutionary and population genetics  analysis. You need to update the paths to the data files in the R scripts. 

`src` folder contains R scripts use to analyze the data.
- `dominant_haplotype_bootstrap.R`: Do bootstrap analysis, define dominant haplotypes, distribution and statistics of domiannt haplotypes. 
- `genetic_diversity_analysis.R`: Phylogenetic tree, nucleotide diversity, selection analysis, Tajima's D, ANI etc.
- `main_tree.R`: Plotting final trees with outgroup
- `input_prep.R`: Preparing input files for Roary
- `phylogenetic_structure.R`: Community ecology analyses
- `tanglegram.R`: Create co-phylogeny of chromosomal and symICE haplotypes

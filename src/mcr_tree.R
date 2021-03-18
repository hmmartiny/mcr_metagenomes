library(tidyverse)
library(ape)
library(igraph)
library(phytools)
library(ggtree)
library(genbankr)
library(gggenes)
library(ggnewscale)
library(RColorBrewer)
library(gtable)
library(grid)
library(kableExtra)
library(ggtreeExtra)

setwd('~/repos/mcr_metagenomes/')

# load data from blast+kraken2+plasmidfinder
assembly_res <- read_csv('data/assembly_results.csv') %>%
  filter(method != 'kraken2') %>%
  mutate(
    sample=run_accession,
    contig=qseqid
  ) %>%
  mutate(
    id = paste(sample, contig, sep='_'),
    gene_variant = str_split(sseqid, '_', simplify = T)[, 1]
  ) %>%
  separate(
    gene_variant, into=c('gene', NA), sep = '\\.', remove = F
  )


# load flanker clusters
flanker_clusters <- read_csv('data/flanker_all.clusters', col_names = F) %>%
  rename(
    guiid=X1,
    group=X2
  ) %>%
  separate(
    guiid, into=c('id', NA),sep = '_mcr', remove=F
  ) %>%
  separate(
    id, into=c('sample', 'contig'), sep=':', remove=F
  ) %>%
  mutate(
    # group = as.numeric(group),
    id = str_replace(id, ':', '_'),
    window = str_split(guiid, '_', simplify = T)[, 9]
  )

# cdhit clusters
cdhit_clusters.list <- read.cdhit('data/all_contigs.cdhit97.res.clstr.clstr')
cdhit_clusters <- tibble()
for (n in names(cdhit_clusters.list)) {
  s <- as_tibble(cdhit_clusters.list[[n]]) %>%
    separate(col = value, into = c('run_accession', 'qseqid'), sep = ':', remove = F) %>%
    mutate(
      cluster=n,
      group=str_split(n, ' ', simplify = T)[, 2],
      id=str_replace(value, ':', '_')
    )
  if(nrow(s) == 1){
    s$group = NA
  } else {
    s$group <- as.numeric(s$group)
  }
  cdhit_clusters <- rbind(cdhit_clusters, s)
}


# add cluster info
x_blast <- data.frame(tree$tip.label) %>%
  mutate(id=tree.tip.label) %>%
  left_join(filter(assembly_res, method=='blastn'), by=c("id" = "id"))

plasmids <- assembly_res %>%
  filter(method == 'plasmidFinder') %>%
  select(id, sseqid)

# read tree
tree <- read.tree('data/all_contigs.dnd')
tree <- midpoint.root(tree)
g <- ggtree(tree)

#colour the tiplabels by mcr gene
(p1 <- g %<+% x_blast + 
    geom_tippoint(aes(color=gene)) 
)
?facet_plot

facet_plot(p1,
           panel="Identified plasmid",
           data=plasmids,
           geom=geom_tile,
           aes(x=1,fill=plasmids$sseqid)
           )
  
p2 <- facet_plot(
  p1, 
  panel='CD-HIT clusters',
  data=select(cdhit_clusters, id, group),
  geom=geom_tile,
  aes(x=1, fill=as.factor(cdhit_clusters$group)) 
) +
  guides(fill=guide_legend(title="CD-HIT cluster"))

p3 <- facet_plot(
  p2, 
  panel='Flanker clusters',
  data = select(flanker_clusters, id, group),
  geom = geom_tile,
  aes(x=1, fill=as.factor(flanker_clusters$group))
) #+     #guides(fill=guide_legend(title="Flanker cluster"))
p3

read.kmadist <- function(distfile){
  num_cols <- as.numeric(read.table(file = kmafile,header = F,nrows = 1))
  mat <- read.table(kmafile, fill=T, skip=1, col.names=rep("", num_cols))
  mat <- as.data.frame(mat)
  rownames(mat) <- mat$X
  return(mat)
}
kmafile='data/all_contigs.KMA.phy'
mat <- read.kmadist(kmafile)
mat
head(mat)

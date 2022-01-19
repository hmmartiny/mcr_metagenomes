#!/usr/bin/env Rscript

library(tidyverse, quietly=TRUE)
library(ggtree, quietly=TRUE)
library(treeio, quietly=TRUE)
library(tidytree, quietly=TRUE)
library(ggnewscale, quietly=TRUE)
library(patchwork, quietly=TRUE)
library(ggtreeExtra, quietly=TRUE)
library(grDevices, quietly=TRUE)
library(argparser, quietly=TRUE)

pargs <- arg_parser("Make pretty plots of trees")
pargs <- add_argument(pargs, "path", help="source path of files needed", type="character")
pargs <- add_argument(pargs, "dest", help="destination path of output files", type="character")
pargs <- add_argument(pargs, '-height', help="height of produced plots", type="numeric", default=7)
pargs <- add_argument(pargs, '-width', help="width of produced plots", type="numeric", default=10)
argv <- parse_args(pargs)


path = argv$path #'~/repos/mcr_metagenomes/data/snp_data/chosen_cov98_depth0.2_AD5.0_AF0.75/'
parameters <- unlist(strex::str_extract_numbers(path, decimals = T))
cov_parameters <- paste("Template Coverage:", parameters[1], ", Depth of Coverage:", parameters[2])
snp_parameters <- paste("SNP AD:", parameters[3], ", SNP AF:", parameters[4])

dest <- argv$dest

caption <- paste("**Minimum thresholds**:", cov_parameters, snp_parameters, sep='<br>')

tree <- read.newick(paste0(path, 'all_consensus_found.aln.tree'), node.label='label')
allele.pattern <- "(mcr-\\d(_|\\.)\\d+)"
variant.pattern <- "(v\\d)$"

snpdists <- read.csv(paste0(path, 'all_consensus_found.aln.snp_dists.tsv'), sep = '\t', row.names = 1, check.names = F)

colnames(snpdists) <- paste(str_extract_all(colnames(snpdists), allele.pattern, simplify = T), str_extract_all(colnames(snpdists), variant.pattern, simplify = T), sep='.')

occurences <- read.csv(paste0(path, 'new_variant_overview.csv'), row.names = 1)

tree.data <- as.tibble(list(tip=tree$tip.label)) %>%
  mutate(
    SeqType=case_when(
      substr(tip, nchar(tip)-1, nchar(tip)-1) == 'v' ~ 'Variant',
      TRUE ~ 'Template',
    ),
    Sequence_no = 
      case_when(
        substr(tip, nchar(tip)-1, nchar(tip)-1) == 'v' ~ substr(tip, nchar(tip), nchar(tip)),
      ),
    gene=substr(tip, 1, 5),
    template=stringr::str_extract(tip, allele.pattern)
)

noplotcols <- grep(pattern = 'min', colnames(occurences))
occurences2 <- occurences %>%
  select(.dots=-grep(pattern = 'min', colnames(occurences)))
occurences.melted <- occurences2 %>%
  tibble::rownames_to_column(var = "row") %>%
   reshape2::melt(id.vars="row") 

gt <- ggtree(tree) %<+% tree.data +
  geom_tiplab(aes(label=template), size=2, hjust = -.15, align = T)+ 
  geom_tippoint(aes(color=gene, shape=SeqType), size=2) +
  new_scale_fill()

gt.occ <- gheatmap(
  p = gt,
  data=occurences2,
  width=2,
  offset = .5,
  colnames = F,
  low = 'white', high='darkblue',
  legend_title = 'sample occurence\n(yes: 1, no: 0)'
) + 
  labs(caption = caption) + 
  theme(plot.caption=ggtext::element_markdown())

ggsave(paste0(dest, 'consensus_tree_occurences.png'), plot=gt.occ, width = argv$width, height=argv$height)

gt.snpdists <- gt + 
  geom_fruit(
    data = reshape2::melt(rownames_to_column(snpdists), id.vars='rowname'),
    geom = geom_tile,
    mapping = aes(x=variable, y=rowname, fill=value),
    offset=0.1,
    pwidt=1,
    color='white'
  ) +
  scale_fill_continuous("SNP\ndistance") +
  labs(caption = caption) + 
  theme(plot.caption=ggtext::element_markdown())
ggsave(paste0(dest, 'consensus_tree_snpdist.png'), plot=gt.snpdists, width = argv$width, height=argv$height)


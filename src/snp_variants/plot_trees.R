#!/usr/bin/env Rscript
options(warn=-1)

suppressMessages(library(tidyverse, quietly=TRUE))
suppressMessages(library(ggtree, quietly=TRUE))
suppressMessages(library(treeio, quietly=TRUE))
suppressMessages(library(tidytree, quietly=TRUE))
suppressMessages(library(ggnewscale, quietly=TRUE))
suppressMessages(library(patchwork, quietly=TRUE))
suppressMessages(library(ggtreeExtra, quietly=TRUE))
suppressMessages(library(grDevices, quietly=TRUE))
suppressMessages(library(argparser, quietly=TRUE))
suppressMessages(library(aplot, quietly=TRUE))

pargs <- arg_parser("Make pretty plots of trees")
# pargs <- add_argument(pargs, "path", help="source path of files needed", type="character")
pargs <- add_argument(pargs, "tree", help="tree file", type="character")
pargs <- add_argument(pargs, "snpdist", help="snpdistance file", type="character")
pargs <- add_argument(pargs, "occ", help="snp-sample occurence file", type="character")
pargs <- add_argument(pargs, "dest", help="destination path of output files", type="character")
pargs <- add_argument(pargs, '-height', help="height of produced plots", type="numeric", default=7)
pargs <- add_argument(pargs, '-width', help="width of produced plots", type="numeric", default=10)
argv <- parse_args(pargs)

dest <- argv$dest

fulltreepath <- tools::file_path_as_absolute(argv$tree)
parameters <- unlist(strex::str_extract_numbers(fulltreepath, decimals = T))
cov_parameters <- paste("Template Coverage:", parameters[1], ", Depth of Coverage:", parameters[2], "Query Identity:", parameters[3], "P-value:", parameters[4])
snp_parameters <- paste("SNP AD:", parameters[5], ", SNP AF:", parameters[6])
caption <- paste("**Minimum thresholds**:", cov_parameters, snp_parameters, sep='<br>')

tree <- read.newick(argv$tree, node.label='label')
allele.pattern <- "(mcr-\\d(_|\\.)\\d+)"
variant.pattern <- "v\\d+"

snpdists <- read.csv(argv$snpdist, sep = '\t', row.names = 1, check.names = F)
snpdists.org <- snpdists
colnames(snpdists) <- paste(str_extract_all(colnames(snpdists), allele.pattern, simplify = T), str_extract_all(colnames(snpdists), variant.pattern, simplify = T), sep='.')

occurences <- read.csv(argv$occ, row.names = 1)

tree.data <- as.tibble(list(tip=tree$tip.label)) %>%
  mutate(
    SeqType=case_when(
      grepl('v\\d+', tip) ~ 'Variant',
      grepl('mcr-', tip) ~ 'Template',
      TRUE ~ 'Variant'
    ),
    Variant_no = str_extract(tip, variant.pattern),
    gene = str_extract(tip,'mcr-\\d+'),
    template=stringr::str_extract(tip, allele.pattern)
) %>%
  fill(gene) %>%
  fill(template) %>%
  mutate(variant=case_when(
    !is.na(Variant_no) ~ str_c(template, Variant_no, sep='.'),
    TRUE ~ template
  ))

noplotcols <- grep(pattern = 'min', colnames(occurences))
occurences2 <- occurences %>%
  select(.dots=-grep(pattern = 'min', colnames(occurences)))
occurences.melted <- occurences2 %>%
  tibble::rownames_to_column(var = "row") %>%
   reshape2::melt(id.vars="row") 

sumOccurences <- occurences2 %>% 
  mutate(sum=rowSums(.)) %>% 
  select(sum) %>% 
  rownames_to_column('tip') %>%
  mutate(gene = str_extract(tip,'mcr-\\d'))



gt <- ggtree(tree) %<+% tree.data +
  geom_tiplab(aes(label=variant), size=2, hjust = -.15, align = F)+ 
  geom_tippoint(aes(color=gene, shape=SeqType), size=2) +
  new_scale_fill() +
  geom_treescale(width=round(mean(tree$edge.length), 2), fontsize=2, linesize = .2, label = 'substitutions per site ', y=5)
ggsave(paste0(dest, 'consensus_tree.png'), plot=gt +labs(caption = caption) + theme(plot.caption=ggtext::element_markdown()), width = argv$width, height=argv$height)


g.occ <- ggplot(data=occurences.melted) + 
  geom_tile(aes(x=variable, y=row, fill=as.factor(value))) +
  scale_fill_manual(
    "sample occurence\n(yes: 1, no: 0)",
    values=c('white', 'darkblue')) +
  labs(caption = caption) + 
  theme(
    axis.text.x = element_text(angle=90, size = 4),
    plot.caption=ggtext::element_markdown()
  )
ggsave(paste0(dest, 'consensus_occurences.png'), plot=g.occ, width = argv$width, height=argv$height)

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

snpdist.melt <- reshape2::melt(rownames_to_column(snpdists.org), id.vars='rowname') %>%
  arrange(rowname, variable)
g.snpdists <- ggplot(data=snpdist.melt) +
  geom_tile(aes(x=variable, y=rowname, fill=value)) +
  geom_text(aes(x=variable, y=rowname, label=value), color='white', size=2) +
  labs(caption = caption) +
  theme(
    axis.text.x = element_text(angle=90, size = 4),
    plot.caption=ggtext::element_markdown()
  )
ggsave(paste0(dest, 'consensus_snpdist.png'), plot=g.snpdists, width = argv$width, height=argv$height)

gt.snpdists <- gt + 
  geom_fruit(
    data = reshape2::melt(rownames_to_column(snpdists), id.vars='rowname'),
    geom = geom_tile,
    mapping = aes(x=variable, y=rowname, fill=value),
    offset=0.5,
    pwidt=1,
    color='white'
  ) +
  scale_fill_continuous("SNP\ndistance") +
  labs(caption = caption) + 
  theme(plot.caption=ggtext::element_markdown())
ggsave(paste0(dest, 'consensus_tree_snpdist.png'), plot=gt.snpdists, width = argv$width, height=argv$height)


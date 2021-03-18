library(Biostrings)
library(karyoploteR)
library(circlize)
library(dplyr)
library(stringr)
library(tidyr)
library(stringr)
library(RColorBrewer)
library(cowplot)
library(ggplot2)

setwd('~/repos/mcr_metagenomes')
source('src/extra_funcs.R')
#clusters <- read.csv('data/all_contigs.cdhit.res.clstr', )


clusters <- read.cdhit('data/all_contigs.cdhit97.res.clstr.clstr')

# read assembly res file
assembly.res <- read.csv('data/assembly_results.csv')
head(assembly.res)

assembly.res <- assembly.res %>%
  mutate(group=case_when(
    str_starts(sseqid, 'mcr') ~ str_split(str_split(sseqid, '\\.', simplify = T)[1], '_', simplify=T)[1]
  )) %>% 
  mutate(group = coalesce(group, sseqid))

create_plots <- function(clusters, assembly.res, prefix='') {
  
  # 433 colors
  color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
  
  for (chosen.cluster in names(clusters)) {
    chosen.data <- as_tibble(clusters[[chosen.cluster]]) %>%
      separate(col = value, into = c('run_accession', 'qseqid'), sep = ':', remove = F) %>%
      mutate(cluster=chosen.cluster) %>%
      left_join(
        assembly.res, by=c('run_accession', 'qseqid')
      ) %>%
      rename(
        chromosome=qseqid,
        start=qstart,
        end=qend,
        header=value,
      )%>%
      mutate(method=str_replace(method, 'blastn', 'mcr')) %>%
      mutate(method=str_replace(method, 'plasmidFinder', 'plasmid'))
    
    # make GRange object
    GR <- GenomicRanges::makeGRangesFromDataFrame(
      filter(chosen.data, !is.na(start)),
      keep.extra.columns = T
    )
    
    pg <- ggbio::autoplot(GR, layout='karyogram', aes(fill=method)) + ggtitle(chosen.cluster)
    p.kar <- pg@ggplot
    p.run <- chosen.data %>% 
      select(chromosome, run_accession) %>%
      distinct() %>%
      ggplot() +
      geom_bar(aes(x=1, y=chromosome, fill=run_accession), stat='identity', width=1) +
      theme_void() +
      theme(panel.spacing.x = unit(1, "mm")) +
      scale_fill_manual(values = color)
    legends <- plot_grid(get_legend(p.run), get_legend(p.kar), ncol = 1)
    p.run <- p.run + theme(legend.position="none")
    p.kar <- p.kar + theme(legend.position="none")
    plot <- plot_grid(p.run, p.kar, align='h', axis='tblr', rel_widths = c(0.5, 15))
    
    g <- plot_grid(plot, legends, ncol=2, nrow=1, rel_widths = c(10, 2))
    
    filename=paste0(prefix, str_replace(chosen.cluster, ' ', '_'), '.pdf')
    ggsave(
      filename=filename, 
      plot=g
    )
  }
}

#create_plots(clusters, assembly.res, prefix='output/contigs_')

# try a cluster
chosen.cluster <- 'Cluster 8'

# combine
chosen.data <- as_tibble(clusters[[chosen.cluster]]) %>%
  separate(col = value, into = c('run_accession', 'qseqid'), sep = ':', remove = F) %>%
  mutate(cluster=chosen.cluster) %>%
  left_join(
    assembly.res, by=c('run_accession', 'qseqid')
  ) %>%
  dplyr::rename(
    chromosome=qseqid,
    start=qstart,
    end=qend,
    header=value,
  )%>%
  mutate(method=str_replace(method, 'blastn', 'mcr')) %>%
  mutate(method=str_replace(method, 'plasmidFinder', 'plasmid')) %>%
  mutate(header2=substr(str_replace(header, ':', '_'), start = 1, stop = 30))

GR <- GenomicRanges::makeGRangesFromDataFrame(
  filter(chosen.data, !is.na(start)),
  keep.extra.columns = T
)

## You need to expand palette size
color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]

pg <- ggbio::autoplot(GR, layout='karyogram', aes(fill=group)) + ggtitle(chosen.cluster)
p.kar <- pg@ggplot
p.run <- chosen.data %>% 
  select(chromosome, run_accession) %>%
  distinct() %>%
  ggplot() +
  geom_bar(aes(x=1, y=chromosome, fill=run_accession), stat='identity', width=1) +
  theme_void() +
  theme(panel.spacing.x = unit(1, "mm")) +
  scale_fill_manual(values = color)
#  scale_fill_brewer(palette="PuBuGn")
legends <- plot_grid(get_legend(p.run), get_legend(p.kar), ncol = 1)
p.run <- p.run + theme(legend.position="none")
p.kar <- p.kar + theme(legend.position="none")
plot <- plot_grid(p.run, p.kar, align='h', axis='tblr', rel_widths = c(0.5, 15))

g <- plot_grid(plot, legends, ncol=2, nrow=1, rel_widths = c(10, 2))

# ALIGN
contig.file <- paste0('data/', str_replace(chosen.cluster, ' ', '_'), '.contigs')
lapply(unique(GR$header), write, contig.file, append=TRUE, ncolumns=1000)

fasta.file <- paste0('data/', str_replace(chosen.cluster, ' ', '_'), '.fasta')
cluster8.aln <- readDNAMultipleAlignment('data/cluster8.aln')
system(
  paste(
    "samtools", "faidx", "data/all_contigs.fasta", '--region-file', contig.file, '-o', fasta.file
  )
)
library("alignfigR")
library(reshape2)
contig.aln <- read_alignment('data/cluster8.aln.fasta.txt')
plot_alignment(contig.aln)

# p <- ggplot() + 
#   geom_rect(
#     plot_frame, 
#     mapping = aes(xmin = x1 - 1, xmax = x2 - 1, ymin = y1 - 1, ymax = y2 - 1, fill = seq), 
#     linetype = 0) + 
#   scale_fill_manual(values = pal, name = legend_title)

x <- as_tibble(contig.aln) %>% t()

xm <- melt(data = x) %>% mutate(
  type=ifelse(value!='-', 'AA', 'GAP')
)
p <- ggplot(xm) + 
  geom_tile(
    mapping = aes(
      x = Var2,
      y = Var1,
      fill=type
    )
  ) +
  scale_fill_manual(values=c("gray", "white"))
#plot <- plot_grid(p.run, p.kar, align='h', axis='tblr', rel_widths = c(0.5, 15))

plot_grid(g, p, align='h', axis='tblr', rel_widths = c(2, 1))

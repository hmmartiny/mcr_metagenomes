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

setwd('~/repos/mcr_metagenomes/')

# load output from flanker clustering
flanker_output<-read_csv('./kpc/all')

#here we just parse some messy filenames - you might need to adjust or be able to omit this
flanker_output<-flanker_output %>% 
  mutate(window = map_chr(isolate, function(s) rev(strsplit(s, "_")[[1]])[3]))

flanker_output$window<-as.numeric(flanker_output$window)
x<-unique(flanker_output$window) %>% sort()
flanker_output<-flanker_output %>% 
  mutate(guuid = map_chr(isolate, function(s) rev(strsplit(s, "_")[[1]])[6]))

flanker_output<-flanker_output %>% 
  mutate(gene = map_chr(isolate, function(s) rev(strsplit(s, "_")[[1]])[5]))

gene<-select(flanker_output,guuid,gene)
flanker_output<-select(flanker_output,guuid,window,group)
#go from long to wide
flanker_output<-flanker_output %>% pivot_wider(id_cols = guuid,names_from=window,values_from=group,names_sort=TRUE)

#give an ID based on the pattern over all 72 windows (7200/100)
flanker_output$ID<-flanker_output %>% group_indices(flanker_output[,2:73])

#take the ids and add them back to the long format for plotting
test<-select(flanker_output,guuid,ID) %>% distinct()

flanker_output<-select(flanker_output,-ID)
flanker_output<-flanker_output %>% pivot_longer(-guuid,names_to = "window",values_to="group")

flanker_output<-left_join(flanker_output,test,by=c("guuid"="guuid"))


flanker_output$window<-as.numeric(flanker_output$window)
flanker_output<-left_join(flanker_output,gene,by=c("guuid"="guuid"))
kpc_store<-flanker_output # we'll use this later



# read tree
tree <- read.tree('data/kpc.tree')
tree <- midpoint.root(tree)
g <- ggtree(tree)

# FLANKER TREE
x_flanker <-data.frame(tree$tip.label) %>%
  mutate(id=tree.tip.label) %>%
  left_join(gene, by=c("id" = "guuid"))

(p1 <- g %<+% x_flanker + geom_tippoint(aes(color=gene)))

px<-p1 + geom_facet(panel="flankergram", data=flanker_output,aes(x=window,fill=group), geom=geom_tile) 
px

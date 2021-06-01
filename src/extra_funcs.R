read.cdhit <- function(file, verbose=TRUE) {
  cluster_lines <- readLines(file)
  clusters = list()
  cluster_seqs=c()
  cluster=NA
  for(i in 1:length(cluster_lines)){
    line = cluster_lines[i]
    
    if(str_starts(line, '>Cluster')){
      if (!is.na(cluster)) {
        clusters[[cluster]] = cluster_seqs
        if (verbose) {
          print(paste(cluster, 'sequences:', length(cluster_seqs)))
        }
        
      }
      
      cluster = str_replace(line, '>', '')
      cluster_seqs = c()
    } else {
      nodeid=str_split(str_split(line, '>', 2, simplify = T)[2],'\\.\\.\\.', 2, simplify = T)[1]
      cluster_seqs = c(cluster_seqs, nodeid)
    }
  }
  return(clusters)
}

read.kmadist <- function(distfile){
  num_cols <- as.numeric(read.table(file = kmafile,header = F,nrows = 1))
  mat <- read.table(kmafile, fill=T, skip=1, col.names=rep("", num_cols))
  mat <- as.data.frame(mat)
  rownames(mat) <- mat$X
  return(mat)
}

fill_kma.dist <- function(mat){
  cols <- colnames(mat)
  mat <- as.matrix(mat)
  for(i in 1:nrow(mat)){
    for(j in 1:ncol(mat)){
      if(i == j){
        mat[i, j] = 0
      } else {
        if(is.na(mat[i, j])){
         tryCatch({mat[i, j] = mat[j, i]}, error = function(e){print(paste(i, j))})
        }
        
      }
    }
  }
  mat <- as.data.frame(mat)
  colnames(mat) <- cols
  return(mat)
}

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

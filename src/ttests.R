library(gt)
library(dplyr)
library(tidyr)
library(plyr)
library(ALDEx2)

run_aldex2 <- function(data, parts, conditions, match) {
  
  # subselect parts that are actually in the data
  parts <- parts[parts %in% colnames(data)]
  
  # select composition and transpose to match aldex2 requirements for data input
  data.t <- round(t(dplyr::select(data, parts)))
  
  # create conditions list
  conditions <- mapvalues(
    conditions == match,
    from = c(0, 1),
    to = c('Other', match)
  )
  
  # run aldex2
  x.all <- aldex(data.t, conditions=conditions, verbose = F)
  
  # format data.frame
  x.all$parts <- rownames(x.all)
  
  conditionA <- str_split(colnames(x.all[2]), '\\.', simplify = T)[3]
  conditionB <- str_split(colnames(x.all[3]), '\\.', simplify = T)[3]
  
  x.all$conditionA <- conditionA
  x.all$conditionB <- conditionB
  
  colnames(x.all)[2] <- str_replace(colnames(x.all)[2], conditionA, "conditionA") 
  colnames(x.all)[3] <- str_replace(colnames(x.all)[3], conditionB, "conditionB") 
  
  
  return(x.all)
}

full_aldex2 <- function(data, parts, condition_column, verbose=F) {
  
  all_conditions <- dplyr::select(data, condition_column)
  condition_counts <- as.data.frame(table(all_conditions))
  unique_conditions <- condition_counts[condition_counts$Freq > 1, 'all_conditions']
  n_conditions <- length(unique_conditions)
  
  res <- data.frame()
  i <- 1
  while (i < n_conditions) {
    tryCatch({
      cond <- as.character(unique_conditions[i])
      if (verbose) {
        print(paste("Condition:", cond))
      }
      cond.res <- run_aldex2(data, parts, conditions = all_conditions, match = cond)
      res <- rbind(res, cond.res)
    }, error = function(e){
      print(paste("Failed to run ALDEx2 for", cond, e))
    }
    )
    
    i <- i + 1
  }
  
  if (nrow(res) > 0) {
    rownames(res) <- seq(1, dim(res)[1])
    res$conditions.A.vs.B <- paste(res$conditionA, "vs", res$conditionB)
  }
  return(res)
}

analyse_aldex <- function(data, parts, totcol, baccol, conds, alpha = 0.05) {
  
  # check argument types 
  if (!is.vector(conds)) {
    stop("'conds' must be a vector!")
  }

  res <- vector(mode="list", length=length(conds))
  names(res) <- conds
  
  org_dat <- data
  
  for (i in 1:length(conds)) {
    print(i)
    # subset to only run on significant res from previous iteration 
    if (i > 1) {
      sig.data.a <- abs.res %>% filter(we.eBH < alpha) %>% dplyr::select(conditionA) %>% distinct() %>% filter(conditionA != 'Other')
      sig.data.b <- abs.res %>% filter(we.eBH < alpha) %>% dplyr::select(conditionB) %>% distinct() %>% filter(conditionB != 'Other')
      sig.cond <- c(sig.data.a$conditionA, sig.data.b$conditionB)
      data <- data %>%
        filter(condition_column %in% sig.cond)
    }
    
    if (nrow(data) == 0) {
      return(res)
    }
    
    cond <- conds[i]
    all.conds <- paste(conds[1:i], collapse= '-')
    assign(paste("res", cond, sep='$'), list)
      
     # create condition column
     data <- data %>% unite(condition_column, conds[1:i], sep='_', remove=F)
     
     # sum parts vs baccol levels
     abs.res <- full_aldex2(
       data = data,
       parts = c(totcol, baccol),
       condition_column = condition_column
     )
     
     # relative level of each part compared
     rel.res <- full_aldex2(
       data = data,
       parts = parts,
       condition_column = condition_column
     )

     # add results to res
     res[[cond]] <- list(
       conditions=all.conds,
       abs.res = abs.res,
       rel.res = rel.res
       )
  }

  return(res)
}

pretty_table <- function(res, alpha = 0.05, star='⭐') {
  t <- res %>%
    filter(we.eBH < alpha) %>%
    mutate(star=star) %>%
    pivot_wider(
      id_cols = conditions.A.vs.B,
      names_from = parts,
      values_from = star
    ) %>% 
    replace(is.na(.), '') %>%
    gt()
  
  return(t)
}

make_effects <- function(xmax, zones=3){
  effect_zones <- as_tibble(list('diff.win'=seq(1, xmax)))
  
  for (z in 0:zones) {
    effect_zones[[paste0('e', 2^z)]] = effect_zones$diff.win * 2^z
    effect_zones[[paste0('em', 2^z)]] = effect_zones$diff.win * (-(2^z))
  }
  
  effect_zones <- effect_zones %>% 
    melt('diff.win', value.name = 'diff.btw') %>%
    mutate(
      'effect.zone' = case_when(
        str_ends(variable, '1') ~ '[0, 1]',
        str_ends(variable, '2') ~ '[1, 2]',
        str_ends(variable, '4') ~ '[2, 4]',
        str_ends(variable, '8') ~ '[4, 8]',
      )
    )
  
  return(effect_zones)
}

gene_formatter <- function(gene, sep='-'){
  return(paste('*', str_replace(gene, sep, '*-'), sep=''))
}

aldex2_plot <- function(res, parts, alpha=0.05, title = '', effect_zones = F, ptext=T, format_genes=T, xlabel='Dispersion', ylabel='Difference', legend_title='Gene', pallette=NA) {
  d <- res %>% 
    mutate(
      group = str_remove(
        str_remove(
          conditions.A.vs.B, 
          'Other vs '),
        ' vs Other'
        
      )
    )
  
  d <- d[which(d$parts %in% parts),]
  
  if (format_genes) {
    d$parts <- gene_formatter(d$parts)
  }
  
  d.sig <- d %>% filter(we.eBH < alpha)
  d.nsig <- d %>% filter(we.eBH >= alpha)

  # get min, max for axis
  xmax <- ceiling(max(d$diff.win))
  xmin <- floor(max(d$diff.win))
  ymax <- ceiling(max(abs(d$diff.btw)))
  
  # initialize plot
  g <- ggplot()
  
  if (effect_zones) {
    g <- g + 
      geom_abline(intercept=0, slope=1, linetype='dashed', alpha=.5) + 
      geom_abline(intercept=0, slope=-1, linetype='dashed', alpha=.5)
  }
  
  g <- g + 
    geom_point(
      data = d.nsig, 
      mapping = aes(x=diff.win, y=diff.btw), color='grey', alpha=.5
    ) +
    new_scale_color() + 
    geom_point(
      data = d.sig,
      mapping = aes(x=diff.win, y=diff.btw, color=group, shape=parts)
    )
  
  if(!is.na(pallette)){
    g <- g + 
      scale_color_manual(values=palette)
  }
  
  if (ptext) {
    g <- g + 
      geom_text(
        data = d.sig, 
        mapping = aes(x=diff.win, y=diff.btw, label=group, color=group), 
        size=2,
        vjust=2
      ) +
      guides(color=F)
  }
  
  # make it pretty
  g <- g +
    xlab(xlabel) + ylab(ylabel) +
    ggtitle(title) + 
    xlim(0, xmax) + ylim(-ymax, ymax) +
    theme_bw() +
    theme(
      plot.title=ggtext::element_markdown(),
      legend.text = ggtext::element_markdown()
    ) +
    guides(
      shape=guide_legend(title=legend_title)
    )
  
  return(g)
}

aldex_ma_plot <- function(res, alpha=0.05, rare=NULL, title="", simplify.cond = TRUE, xlab=NULL, ylab=NULL) {
  # Adaption of MA aldex.plot
  
  if (is.null(xlab)) {
    xlab <- expression("Median" ~ ~Log[2] ~ ~"relative abundance")
  }
  if (is.null(ylab)) {
    ylab <- expression("Median" ~ ~Log[2] ~ ~"Difference")
  }
  
  # simplify how condition legend is written
  if (simplify.cond) {
    res <- res %>% 
      mutate(conditions = str_replace(conditions.A.vs.B, ' vs Other', '')) %>% 
      mutate(conditions = str_replace(conditions, 'Other vs ', ''))
    
  } else {
    res$conditions <- res$conditions.A.vs.B
  }
  
  shapes = c("R" = 20, "S" = 8, "NS" = 1)

  # split data
  sig.data <- res %>% filter(we.eBH <= alpha)
  rare.data <- res %>% filter(rab.all < rare)
  nsig.data <- res %>% filter(we.eBH > alpha, rab.all > rare)
  
  # plot
  g <- ggplot() +
    geom_point(
      data = sig.data, aes(x = rab.all, y = diff.btw, color=conditions, shape="S"),
      size = 1.5
    ) + 
    geom_point(
      data = rare.data, aes(x = rab.all, y = diff.btw, color=conditions, shape="R"),
      size = 1
    ) +
    geom_point(
      data = nsig.data, aes(x = rab.all, y = diff.btw, color=conditions, shape="NS"),
      size = 1,
      alpha = 0.5
    ) + 
    scale_shape_manual(values=shapes) +
    facet_grid(. ~ parts, scales='free') +
    xlab(xlab) + 
    ylab(ylab) +
    ggtitle(title)
  
  return(g)
}

aldex_mw_plot <- function(res, alpha=0.05, rare=NULL, title="", xlab=NULL, ylab=NULL, lncol=1) {
  # Adaption of MW aldex.plot
  
  if (is.null(xlab)) {
    xlab <- expression("Median" ~ ~Log[2] ~ ~"Dispersion")
  }
  if (is.null(ylab)) {
    ylab <- expression("Median" ~ ~Log[2] ~ ~"Difference")
  }
  
  shapes = c("R" = 20, "S" = 8, "NS" = 1)

  sig.data <- res %>% filter(we.eBH <= alpha)
  nsig.data <- res %>% filter(we.eBH > alpha)
  
  g <- ggplot()
  
  g <- g + geom_point(
      data = sig.data, aes(x = diff.win, y = diff.btw, color=conditions.A.vs.B, shape= "S"),
      size = 1.5
    )
  if (!is.null(rare)) {
    rare.data <- res %>% filter(rab.all < rare)
    g <- g + geom_point(
      data = rare.data, aes(x = diff.win, y = diff.btw, color=conditions.A.vs.B, shape="R"),
      size = 1
    ) 
    nsig.data <- nsig.data %>% filter(rab.all > rare)
  }
  g <- g + geom_point(
      data = nsig.data, aes(x = diff.win, y = diff.btw, shape="NS"),
      size = 1,
      alpha = 0.5,
      color = 'grey'
      )
  
  # prettify
  first_facet <- sort(unique(res$parts))[1]
  x0 <- min(filter(res, parts == first_facet)$diff.win)
  y1 = min(res$diff.btw)
  y2 = max(res$diff.btw)
  
  d1 <- data.frame(parts=c(first_facet), x=c(x0), y=c(y1), c=c("condition A"))
  d2 <- data.frame(parts=c(first_facet), x=c(x0), y=c(y2), c=c("condition B"))
  g <- g + coord_cartesian(clip = 'off') + 
    geom_text(data = d1, aes(x=x, y=y, label=c), color='grey', size=3, angle=90, vjust=-.8, hjust=.3) +
    geom_text(data = d2, aes(x=x, y=y, label=c), color='grey', size=3, angle=90, vjust=-.8, hjust=.3)
  
  
  # labels
  g <- g + 
    scale_shape_manual(values=shapes) +
    facet_grid(. ~ parts, scales='free') +
    xlab(xlab) + 
    ylab(ylab) +
    ggtitle(title) +
    guides(color=guide_legend(ncol=lncol, title="conditions A vs B"), shape=guide_legend(title='Significance')) +
    theme_bw() +
    theme(axis.text.x=element_text(size=rel(0.5)))
  
  return(g)
}

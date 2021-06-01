add_arrows <- function(g, loadings, add_text=T, new_scale=T) {
  
  # add new coloring scheme for arrows
  if (new_scale) {
    g <- g + new_scale_color()
  }
  
  # draw arrows with geom_segment
  g <- g + 
    geom_segment(
      data = loadings,
      mapping = aes(x=0, y=0, xend=PC1, yend=PC2, color=name),
      arrow = arrow(length=unit(0.3, 'cm'))
    ) +
    scale_color_brewer(palette='Set1')
  
  # annotate arrows by putting a label at the arrow head
  if (add_text) {
    g <- g + geom_text(
      data = loadings,
      mapping = aes(x=PC1, y=PC2, label=name, color=name),
      size = 3, 
      nudge_x=-0.5,
      nudge_y=-1,
      vjust="inward",
      hjust="inward"
    ) +
      guides(
        color=F
      ) +
      coord_cartesian(clip='off')
  }
  
  return(g)
}

add_scores <- function(g, scores, x, y, color=NA, shape=16, add_text=F, alpha=.5,  new_scale=T, color_label=NA, n_colors=NA, color_wheel=NA) {
  
  # create new color scale
  if (new_scale) {
    g <- g + new_scale_color()
  }
  
  # if no color label is given, color all points grey
  # same goes for NA values, they will be colored grey
  if (is.na(color)){
    g <- g + geom_point(
      data = scores,
      mapping = aes_string(x=x, y=y), 
      alpha=alpha, shape=16, color='grey'
    )
  } else {
    
    if (is.na(color_wheel)) {
      if (is.na(n_colors)) {
        n_colors <- length(unique(scores[[color]])) 
      }
      color_wheel <- unname(c(pals::alphabet(), pals::alphabet2()))[1:n_colors]
    }
  
    # if shape label is given, give points different shapes according to label
    if (is.character(shape)) {
      n_shapes <- length(unique(scores[[shape]]))
      g <- g + geom_point(
        data = scores,
        mapping = aes_string(x=x, y=y, color=color, shape=shape), 
        alpha=alpha
      )+ 
        scale_shape_manual(values=1:n_shapes, na.value=16)
    } else {
      shapes <- 16:16
      g <- g + geom_point(
        data = scores,
        mapping = aes_string(x=x, y=y, color=color), 
        alpha=alpha, shape=16
      )
    }
    
    if (add_text) {
      # get "extreme" points to add text inside of plot.
      # hopefully they will be outliers, so it will be easier for the viewer
      # to see which color belongs to which label

      sample_score <- tryCatch({
        df <- scores %>% 
          group_by(get(color)) %>%
          mutate(PC=abs(get(x))+abs(get(y))) %>%
          filter(PC == max(PC))
        },
        error = function(err) {
          df <- scores %>% 
            group_by(color) %>%
            mutate(PC=abs(get(x))+abs(get(y))) %>%
            filter(PC == max(PC))
        }
      )
      
      # annotate
      g <- g + 
        scale_color_manual(guide='none', values=color_wheel, na.value='grey') + 
        geom_text_repel(
          data = sample_score,
          mapping = aes_string(x=x, y=y, color=color, label=color),
          size = 2,
          box.padding = unit(.75, "lines"),
          max.overlaps = Inf
        )
    } else {
      if (is.na(color_label)) {
        color_label <- str_to_title(str_replace(color, '_', ' '))
      }
      g <- g + scale_color_manual(color_label, values=color_wheel, na.value='grey')
    }
  }
  
  if(is.character(shape)){
    g <- g + guides(
      shape=guide_legend(
        title = str_to_title(str_replace(shape, '_', ' '))
      )
    )
  }
  return(g)
}

draw_ellipses <- function(g, scores, groupby, x, y, scale=5.991, n_min=3){

  g <- g + 
      new_scale_fill() + new_scale_color()
  
  group_counts <- scores %>% group_by(get(groupby)) %>% tally(sort=T)
  cgroups <- group_counts[['get(groupby)']] #unique(scores[[groupby]])
  n_groups <- nrow(cgroups)
  colors <- colorRampPalette(brewer.pal(name='Set3', n=8))(n_groups)

  for (cgroup in cgroups) {
    filtered_score <- scores[scores[[groupby]] == cgroup, c(x, y)]
    if(nrow(filtered_score) > n_min) {
      ei <- eigen(cov(filtered_score))
      lambda <- sqrt(ei$values)
      angle <- ei$vectors
      
      ellipse_row <- tibble(
        groupby = cgroup,
        a = lambda[1] * sqrt(scale),
        b = lambda[2] * sqrt(scale),
        x0 = mean(filtered_score[[x]]),
        y0 = mean(filtered_score[[y]]),
        angle = atan(angle[2,1] / angle[1, 1])
      )
      
      g <- g + ggforce::geom_ellipse(
        data = ellipse_row,
        mapping = aes(x0=x0, y0=y0, a=a, b=b, angle=angle, color=groupby, fill=groupby),
        alpha = .35
      )
    }
  }
  
  n_colors <- length(unique(scores[scores$host %in% as.list(sig.hosts)$host, 'REGION']))
  my_colors <- colorRampPalette(brewer.pal(name="Set3", n = 12))(n_colors)
  
  legend_title <- str_to_title(str_replace(groupby, '_', ' '))
  g <- g + scale_fill_manual(legend_title, values=my_colors) + 
      scale_color_manual(legend_title, values=my_colors)
  
  
  return(g)
}


make_biplot <- function(scores, loadings, filterby, color, x, y, sig_values=c(), ellipses=NA, shape=NULL, scores_text=F, arrows_text=T, title=waiver(), xlabel=waiver(), ylabel=waiver(), ellipses_nmin=3, scores_palette=NA) {
  g <- ggplot()
  
  if (!is.na(ellipses)) {
    g <- draw_ellipses(g = g, scores = scores, x=x, y=y, groupby = ellipses)
  }
  
  # add scatter points
  if (length(sig_values) > 0) {
    # non significant
    g <- add_scores(
      g = g, 
      scores = scores[!(scores[[filterby]] %in% sig_values), ], 
      x=x, y=y, color=NA, 
      new_scale=T
    )
    
    # significant
    g <- add_scores(
      g = g, 
      scores = scores[(scores[[filterby]] %in% sig_values), ], 
      x=x, y=y, color=color,  alpha = .75, shape=shape, 
      new_scale=T, add_text = scores_text,
      color_wheel = scores_palette
    )
    
  } else {
    g <- add_scores(
      g = g, 
      scores = scores, 
      x=x, y=y, color=color, shape=shape,
      new_scale=T, add_text=scores_text,
      color_wheel = scores_palette
    )
  }
  
  # add arrows
  g <- add_arrows(g = g, loadings = loadings, add_text=arrows_text, new_scale=T)
  
  g <- g + 
    theme_bw() + 
    labs(
      title = title,
      x = xlabel,
      y = ylabel
    )
  
  return(g)
}

make_adv_biplot <- function(scores, loadings, sigs, x, y, color, ellipses, shape, scores_text=F, arrows_text=T, title=waiver(), xlabel=waiver(), ylabel=waiver(), ellipses_nmin=3, alpha=.5) {
  
  g <- ggplot()
  
  # draw ellipses if specified
  if (!is.na(ellipses)) {
    filteron=ellipses
    if (ellipses == 'REGION') {
      filteron='country'
    } 
    
    g <- draw_ellipses(g = g, scores = scores[scores[[filteron]] %in% sigs[[filteron]], ], x=x, y=y, groupby = ellipses)
  }
  
  # add scatter points
  scores$color <- NA #'Not significant'
  scores$shape <- NA #'Not significant'
  
  scores[scores[[color]] %in% sigs[[color]], 'color'] <- scores[scores[[color]] %in% sigs[[color]], color] 
  scores[scores[[shape]] %in% sigs[[shape]], 'shape'] <- scores[scores[[shape]] %in% sigs[[shape]], shape]
  
  # add NA values (not significant)
  g <- add_scores(
    g = g,
    scores = scores[is.na(scores[['color']]),],
    x = x, y = y,
    shape='shape',
    new_scale=T
  ) + scale_color_discrete(guide = 'none')

  g <- add_scores(
    g = g,
    scores = scores[is.na(scores[['shape']]),],
    x = x, y = y,
    color='color',
    add_text = scores_text
  ) + scale_shape(guide='none')
  
  # add significant groups
  g <- add_scores(
    g = g,
    scores = filter(scores, !is.na(color), !is.na(shape)),
    x=x, y=y, color='color', shape='shape',
    new_scale=F, add_text=scores_text, 
    color_label = str_to_title(str_replace(color, '_', ' ')), 
    n_colors = length(sigs[[color]])
  )
  
  # add arrows
  g <- add_arrows(g = g, loadings = loadings, add_text=arrows_text, new_scale=T)
  
  g <- g + 
    theme_bw()
    labs(
      title = title,
      x = xlabel,
      y = ylabel
    )
  
  return(g)
}

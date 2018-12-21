
plot_interaction <- function(x, 
                             xlab = "", ylab = "default", main = "default",
                             them,
                             point.size = 2.5, line.size = 1,
                             palette = "Set1", ...){
  
  if(main == "default")
    main <- sprintf("Point estimates with %s%%-CI", (1-x$level)*100)

  if(ylab == "default")
    ylab <- switch(x$fitfunc,
                   "coxph" = "Treatment effect differences",
                   "Treatment effect differences")

  data <- x$trtEffDiff
  
  if(x$type == "bagged") {
    selected = which(x$boot_results$unadjusted$selected == 1)
    selGroup = paste(x$boot_results$unadjusted$subgroup[selected])
    if (selGroup == "none") return("No selected subgroup")
    data = data[which(data$Group == selGroup), ]
  }
  data$Group <- factor(data$Group, levels = unique(data$Group))
  # over <- x$overall[c(2,1,3)]
  if(missing(them)){
    them <- .subtee.them
  }
    p <- ggplot(data,
                aes_string(x=as.factor(1), y="trtEffDiff",
                           ymin="LB", ymax="UB")) + 
      facet_grid(Group ~ .) +
      geom_point(size=point.size) +
      geom_linerange(size=line.size, show.legend=FALSE) + 
      coord_flip() +
      xlab(xlab) + ylab(ylab) + ggtitle(main) +
      scale_colour_brewer(palette = palette) + 
      them + theme(axis.text.y=element_blank(),
                   strip.text.y=element_text(angle=0)) +
      guides(colour=FALSE)
  xmin <- 0.5;xmax <- 1.5
  p <- p + geom_hline(yintercept=0,linetype=2) 
  p
}

plot_effect <- function(x, show.compl = FALSE, 
                        xlab = "", ylab = "", main = "default",
                        them,
                        point.size = 2.5, line.size = 1,
                        palette = "Set1", ...){
  
  stopifnot(is.logical(show.compl), 
            is.numeric(point.size), is.numeric(line.size))
  
  if(main == "default")
    main <- sprintf("Point estimates with %s%%-CI", (1-x$level)*100)
  
  if(ylab == "default")
    ylab <- switch(x$fitfunc,
                   "coxph" = "Treatment effect estimates (log-hazard ratio)",
                   "Treatment effect estimates")
  
  data <- x$trtEff
  
  if(x$type == "bagged") {
    selected = which(x$boot_results$unadjusted$selected == 1)
    selGroup = paste(x$boot_results$unadjusted$subgroup[selected])
    if (selGroup == "none") return("No selected subgroup")
    data = data[which(data$Group == selGroup), ]
  }
  data$Group <- factor(data$Group, levels = unique(data$Group))
  over <- x$overall[c(2,1,3)]
  if(missing(them)){
    them <- .subtee.them
  }
  if(show.compl){
    p <- ggplot(data,
                aes_string(x="Subset", y="trtEff",
                           ymin="LB", ymax="UB", colour="Subset")) + 
      facet_grid(Group ~ .) +
      geom_point(size=point.size) +
      geom_linerange(size=line.size, show.legend=FALSE) + 
      coord_flip() +
      xlab(xlab) + ylab(ylab) + ggtitle(main) +
      scale_colour_brewer(palette = palette, name="Group")+
      them + theme(axis.text.y=element_blank(),
                   strip.text.y=element_text(angle=0))+
      guides(colour = guide_legend(reverse = TRUE))
  } else {
    p <- ggplot(data[data$Subset=="Subgroup",],
                aes_string(x="Subset", y="trtEff",
                           ymin="LB", ymax="UB", colour="Subset")) + 
      facet_grid(Group ~ .) +
      geom_point(size=point.size) +
      geom_linerange(size=line.size, show.legend=FALSE) + 
      coord_flip() +
      xlab(xlab) + ylab(ylab) + ggtitle(main) +
      scale_colour_brewer(palette = palette) + 
      them + theme(axis.text.y=element_blank(),
                   strip.text.y=element_text(angle=0)) +
      guides(colour=FALSE)
  }
  if(show.compl){
    xmin <- 0.5;xmax <- 2.5
  } else {
    xmin <- 0.5;xmax <- 1.5
  }
  p <- p + geom_hline(yintercept=over[2],linetype=2) + 
    geom_rect(aes(xmin=xmin,
                  xmax=xmax,
                  ymin=over[1],ymax=over[3]),
              fill=rgb(0, 0, 0, 0.05), color=rgb(0, 0, 0, 0.05))
  p
}

plot_compare_interaction <- function(x, 
                                     xlab = "", ylab = "", main = "default",
                                     them,
                                     point.size = 2.5, line.size = 1,
                                     palette = "Dark2", ...){
  n_comp = length(x)
  if (!(n_comp %in% 2:3)) {
    stop("The argument x should be a list of 2 or 3 components.")
  }
  types = sapply(x, function(res) {
    res$type
  })
  levels = sapply(x, function(res) {
    res$level
  })
  if (length(unique(levels)) != 1){
    stop("Different levels of confidence were used.")
  }
  if(main == "default")
    main <- sprintf("Point estimates with %s%%-CI", (1-x[[1]]$level)*100)
  
  
  if(ylab == "default")
    ylab <- switch(x[[1]]$fitfunc,
                   "coxph" = "Treatment effect differences",
                   "Treatment effect differences")
  
  data = do.call(rbind, lapply(x, function(results) {
    data.frame(results$trtEffDiff, type = results$type)
  }))
  
  if (any(types == "bagged")){
    selected = 0
    for (i in 1:n_comp){
      if (x[[i]]$type != "bagged") next()
      selected = which(x[[i]]$boot_results$unadjusted$selected == 1)
      selGroup = paste(x[[i]]$boot_results$unadjusted$subgroup[selected])
      if (selGroup == "none") return("No selected subgroup for the bagged method")
      data = data[which(data$Group == selGroup), ]
    }
  }
  
  data$Group <- factor(data$Group, levels = unique(data$Group))
  
  if(missing(them)){
    them <- .subtee.them
  }
  p <- ggplot(data,
              aes_string(x="type", y="trtEffDiff",
                         ymin="LB", ymax="UB", colour = "type")) + 
    facet_grid(Group ~ .) +
    geom_point(size=point.size) +
    geom_linerange(size=line.size, show.legend=FALSE) + 
    coord_flip() +
    xlab(xlab) + ylab(ylab) + ggtitle(main) +
    scale_colour_brewer(palette = palette) + 
    them + theme(axis.text.y=element_blank(),
                 strip.text.y=element_text(angle=0)) +
    guides(colour = guide_legend(reverse = TRUE))
  if(n_comp == 3){
    xmin <- 0.5;xmax <- 3.5
  } else {
    xmin <- 0.5;xmax <- 2.5
  }
  p <- p + geom_hline(yintercept=0,linetype=2) 
  p
}


plot_compare_effect <- function(x, show.compl = FALSE, 
                                     xlab = "", ylab = "", main = "default",
                                     them,
                                     point.size = 2.5, line.size = 1,
                                     palette = "Dark2", ...){
  
  stopifnot(is.logical(show.compl), 
            is.numeric(point.size), is.numeric(line.size))
  
  n_comp = length(x)
  if (!(n_comp %in% 2:3)) {
    stop("The argument x should be a list of 2 or 3 components.")
  }
  types = sapply(x, function(res) {
    res$type
  })
  levels = sapply(x, function(res) {
    res$level
  })
  if (length(unique(levels)) != 1){
    stop("Different levels of confidence were used.")
  }
  if(main == "default")
    main <- sprintf("Point estimates with %s%%-CI", (1-x[[1]]$level)*100)
  
  if(ylab == "default")
    ylab <- switch(x[[1]]$fitfunc,
                   "coxph" = "Treatment effect estimates (log-hazard ratio)",
                   "Treatment effect estimates")
  
  data = do.call(rbind, lapply(x, function(results) {
    data.frame(results$trtEff, type = results$type)
  }))
  
  if (any(types == "bagged")){
    selected = 0
    for (i in 1:n_comp){
      if (x[[i]]$type != "bagged") next()
      selected = which(x[[i]]$boot_results$unadjusted$selected == 1)
      selGroup = paste(x[[i]]$boot_results$unadjusted$subgroup[selected])
      if (selGroup == "none") return("No selected subgroup for the bagged method")
      data = data[which(data$Group == selGroup), ]
    }
  }
  
  data$Subgroup <- paste0(data$Group, " (", substr(data$Subset, 1, 1), ")")  
  data$Subgroup <- factor(data$Subgroup, levels = unique(data$Subgroup))
  data$Group <- factor(data$Group, levels = unique(data$Group))
  over <- x[[1]]$overall[c(2,1,3)]
  if(missing(them)){
    them <- .subtee.them
  }
  
  if(show.compl){
    p <- ggplot(data,
                aes_string(x="type", y="trtEff",
                           ymin="LB", ymax="UB", colour="type")) + 
      facet_grid(Subgroup ~ .) +
      geom_point(size=point.size) +
      geom_linerange(size=line.size, show.legend=FALSE) + 
      coord_flip() +
      xlab(xlab) + ylab(ylab) + ggtitle(main) +
      scale_colour_brewer(palette = palette, name="Group")+
      them + theme(axis.text.y=element_blank(),
                   strip.text.y=element_text(angle=0))+
      guides(colour = guide_legend(reverse = TRUE))
  } else {
    p <- ggplot(data[data$Subset=="Subgroup",],
                aes_string(x="type", y="trtEff",
                           ymin="LB", ymax="UB", colour="type")) + 
      facet_grid(Group ~ .) +
      geom_point(size=point.size) +
      geom_linerange(size=line.size, show.legend=FALSE) + 
      coord_flip() +
      xlab(xlab) + ylab(ylab) + ggtitle(main) +
      scale_colour_brewer(palette = palette) + 
      them + theme(axis.text.y=element_blank(),
                   strip.text.y=element_text(angle=0)) +
      guides(colour = guide_legend(reverse = TRUE))
  }
  if(n_comp == 3){
    xmin <- 0.5;xmax <- 3.5
  } else {
    xmin <- 0.5;xmax <- 2.5
  }
  p <- p + geom_hline(yintercept=over[2],linetype=2) + 
    geom_rect(aes(xmin=xmin,
                  xmax=xmax,
                  ymin=over[1],ymax=over[3]),
              fill=rgb(0, 0, 0, 0.05), color=rgb(0, 0, 0, 0.05))
  p
}


.subtee.them <- theme(strip.text = element_text(size = rel(1.0),
                                        face = "bold"),
              legend.title = element_text(size=rel(1.15),
                                          face = "bold"),
              legend.text = element_text(size=rel(1.00),
                                         face = "bold"),
              axis.ticks = element_blank(),
              axis.title.x = element_text(size=rel(1.15),
                                          face = "bold"),
              axis.title.y = element_text(size=rel(1.15),
                                          face = "bold"),
              axis.text.x = element_text(size=rel(1.00),
                                         face = "bold"),
              axis.text.y = element_text(size=rel(1.00),
                                         face = "bold"),
              plot.title = element_text(size=rel(1.55),
                                        lineheight = 1,
                                        face = "bold"))

ideo_single <- function(genome="hg38",chr="chr5",yrange=c(0,10))
{

  k <-  biovizBase::getIdeogram(genome = genome,subchr = chr,cytobands = TRUE)
  df <- k %>% as.data.frame()

  cytobandColor <- getOption("biovizBase")$cytobandColor
  df.rect <- subset(df, gieStain != "acen")
  df.tri <- subset(df, gieStain == "acen")
  df.tri.p <- df.tri[substr(df.tri$name, 1, 1) == "p",]
  df.tri.q <- df.tri[substr(df.tri$name, 1, 1) == "q",]

  ## main
  .ideo.range <- yrange
  p.ideo <- list(do.call(ggplot2::geom_rect, c(list(data = df.rect),
                                               list(do.call(aes,list(xmin = as.name("start"),
                                                                     ymin =.ideo.range[1],
                                                                     xmax = as.name("end"),
                                                                     ymax = .ideo.range[2],
                                                                     fill = as.name("gieStain")))),
                                               list(color = NA, alpha = 0.7))))

  ## draw line
  df.p <- df.rect[substr(df.rect$name, 1, 1) == "p",]
  df.q <- df.rect[substr(df.rect$name, 1, 1) == "q",]

  if(nrow(df.p)){

    df.p.d <- do.call(rbind, by(df.p, df.p$seqnames, function(dd){
      data.frame(x = min(dd$start),
                 y = .ideo.range[1],
                 y2 = .ideo.range[2],
                 xend = max(dd$end),
                 yend = .ideo.range[1],
                 yend2 = .ideo.range[2],
                 seqnames = unique(dd$seqnames))
    }))


    p.ideo <- c(p.ideo, list(do.call(ggplot2::geom_segment, c(list(data = df.p.d),
                                                              list(aes(x = x, y = y, xend = xend, yend = yend)),
                                                              list(color = "black",
                                                                   alpha = 1, size = 0.3)))))
    p.ideo <- c(p.ideo, list(do.call(ggplot2::geom_segment, c(list(data = df.p.d),
                                                              list(aes(x = x, y = y2, xend = xend, yend = yend2)),
                                                              list(color = "black",
                                                                   alpha = 1, size = 0.3)))))
    p.ideo <- c(p.ideo, list(do.call(ggplot2::geom_segment, c(list(data = df.p.d),
                                                              list(aes(x = x, y = y, xend = x, yend = y2)),
                                                              list(color = "black",
                                                                   alpha = 1, size = 0.3)))))
  }

  if(nrow(df.q)){
    df.q.d <- do.call(rbind, by(df.q, df.q$seqnames, function(dd){
      data.frame(x = min(dd$start),
                 y = .ideo.range[1],
                 y2 = .ideo.range[2],
                 xend = max(dd$end),
                 yend = .ideo.range[1],
                 yend2 = .ideo.range[2],
                 seqnames = unique(dd$seqnames))
    }))


    p.ideo <- c(p.ideo, list(do.call(ggplot2::geom_segment, c(list(data = df.q.d),
                                                              list(aes(x = x, y = y, xend = xend, yend = yend)),
                                                              list(color = "black",
                                                                   alpha = 1, size = 0.3)))))
    p.ideo <- c(p.ideo, list(do.call(ggplot2::geom_segment, c(list(data = df.q.d),
                                                              list(aes(x = x, y = y2, xend = xend, yend = yend2)),
                                                              list(color = "black",
                                                                   alpha = 1, size = 0.3)))))

    p.ideo <- c(p.ideo, list(do.call(ggplot2::geom_segment, c(list(data = df.q.d),
                                                              list(aes(x = xend, y = y, xend = xend, yend = y2)),
                                                              list(color = "black",
                                                                   alpha = 1, size = 0.3)))))





  }

  if(nrow(df.q)){
    df.q.d <- do.call(rbind, by(df.q, df.q$seqnames, function(dd){
      data.frame(x = min(dd$start),
                 y = .ideo.range[1],
                 y2 = .ideo.range[2],
                 xend = max(dd$end),
                 yend = .ideo.range[1],
                 yend2 = .ideo.range[2],
                 seqnames = unique(dd$seqnames))
    }))


    p.ideo <- c(p.ideo, list(do.call(ggplot2::geom_segment, c(list(data = df.q.d),
                                                              list(aes(x = x, y = y, xend = xend, yend = yend)),
                                                              list(color = "black",
                                                                   alpha = 1, size = 0.3)))))
    p.ideo <- c(p.ideo, list(do.call(ggplot2::geom_segment, c(list(data = df.q.d),
                                                              list(aes(x = x, y = y2, xend = xend, yend = yend2)),
                                                              list(color = "black",
                                                                   alpha = 1, size = 0.3)))))

    p.ideo <- c(p.ideo, list(do.call(ggplot2::geom_segment, c(list(data = df.q.d),
                                                              list(aes(x = xend, y = y, xend = xend, yend = y2)),
                                                              list(color = "black",
                                                                   alpha = 1, size = 0.3)))))





  }

  df.tri.p2 <- with(df.tri.p,
                    data.frame(x=start,
                               y=rep(0, nrow(df.tri.p)),
                               xend=start,
                               yend=rep(10, nrow(df.tri.p)),
                               height=abs(start - end),
                               seqnames=seqnames, strand=strand,
                               name=name, gieStain=gieStain))

  df.tri.q2 <- with(df.tri.q,
                    data.frame(x=end,
                               y=rep(0, nrow(df.tri.q)),
                               xend=end,
                               yend=rep(10, nrow(df.tri.q)),
                               height=-abs(start - end),
                               seqnames=seqnames, strand=strand,
                               name=name, gieStain=gieStain))

  if(nrow(df.tri.p2) > 0L)
    p.ideo <- c(p.ideo,
                list(geom_arch_flip2(df.tri.p2,
                                     aes(x = x,
                                         y = y ,
                                         xend = xend,
                                         yend = yend,
                                         height = height),
                                     color = "black", size = 0.5),
                     geom_arch_flip(df.tri.p2,
                                    aes(x = x,
                                        y = y ,
                                        xend = xend,
                                        yend = yend,
                                        height = height,
                                        fill = gieStain))))

  if(nrow(df.tri.p2) > 0L)
    p.ideo <- c(p.ideo,
                list(geom_arch_flip2(df.tri.q2,
                                     aes(x = x,
                                         y = y ,
                                         xend = xend,
                                         yend = yend,
                                         height = height),
                                     color = "black",
                                     size = 0.5),
                     geom_arch_flip(df.tri.q2,
                                    aes(x = x,
                                        y = y ,
                                        xend = xend,
                                        yend = yend,
                                        height = height,
                                        fill = gieStain))))

  p.ideo1 <- c(p.ideo,
               list(theme(axis.text.y = element_blank(),
                          axis.title.y=element_blank(),
                          axis.ticks = element_blank(),
                          axis.title.x=element_blank(),
                          axis.text.x=element_blank(),
                          axis.ticks.x=element_blank(),
                          legend.position="none",
                          panel.background = element_rect(fill = "white"),
                          panel.grid.minor = element_line(colour = NA),
                          panel.grid.major = element_line(colour = NA)),
                    scale_fill_manual(values = cytobandColor)))
}


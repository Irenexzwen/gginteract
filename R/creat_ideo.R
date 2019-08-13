#' Title Create customized ideogram ggplot object
#'
#' @param genome String. Genome version eg."hg38","mm10","mm20". default("hg38")
#' @param chr String. Subchr of \code{geome} eg."chr4"
#' @param ideo.width Numeric. Width of the ideogram ggplot object.
#' @param ydrift Numeric. x axis dirft from 0.
#' @param xdrift Numeric. y axis dirft from 0.
#' @param whratio Numeric, the width and height ratio of the ideogram.
#'
#' @return ideo class
#' @export
#'
#' @examples
#' ideo <- create_ideo(chr="chr4",ideo.width=400)
#'
#' @include gginter_utils.R
#' @include AllClass.R
create_ideo <- function(genome="hg38",chr,ideo.width,ydrift=0,xdrift=0,whratio=16)
{

  k <-  biovizBase::getIdeogram(genome = genome,subchr = chr,cytobands = TRUE)
  df <- k %>% as.data.frame()

  chr_start <- 1
  chr_end <- max(df$end)

  ideo.height <- ideo.width/whratio

  if(ydrift>0){
    .ideo.range <- c(ydrift,ydrift+ideo.height)
    }
  else{
    .ideo.range <- c(ydrift-ideo.height,ydrift)
  }

  scale <- max(df$end)/ideo.width
  df$start <- df$start/scale + xdrift
  df$end <-   df$end/scale + xdrift
  df$width <- df$width/scale

  left <- min(df$start)
  right <- max(df$end)

  cytobandColor <- getOption("biovizBase")$cytobandColor
  df.rect <- subset(df, gieStain != "acen")
  df.tri <- subset(df, gieStain == "acen")
  df.tri.p <- df.tri[substr(df.tri$name, 1, 1) == "p",]
  df.tri.q <- df.tri[substr(df.tri$name, 1, 1) == "q",]

  ## main
  p.ideo <- list(do.call(ggplot2::geom_rect, c(list(data = df.rect),
                                               list(do.call(ggplot2::aes,list(xmin = as.name("start"),
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
                               y=rep(.ideo.range[1], nrow(df.tri.p)),
                               xend=start,
                               yend=rep(.ideo.range[2], nrow(df.tri.p)),
                               height=abs(start - end),
                               seqnames=seqnames, strand=strand,
                               name=name, gieStain=gieStain))

  df.tri.q2 <- with(df.tri.q,
                    data.frame(x=end,
                               y=rep(.ideo.range[1], nrow(df.tri.q)),
                               xend=end,
                               yend=rep(.ideo.range[2], nrow(df.tri.q)),
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

  return(new("ideo",
             geom_ideobody = p.ideo1,
             chr = chr,
             height = ideo.height,
             xdrift = xdrift,
             ydrift = ydrift,
             chr_start = chr_start,
             chr_end = chr_end,
             .plot_scale = scale,
             .ideo_left = left,
             .ideo_right = right
             ))

}


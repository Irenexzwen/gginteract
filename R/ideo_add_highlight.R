# ideo highlight add a highlight region on an ideaogram geom_list,height,xdrift,ydrift

#' Title
#'
#' @param ideo ideo object, object to add a highlight region.
#' @param region list, the region to highlight,
#' @param tick_len numeric, the vertical height of the highligh region.
#'
#' @return ideo object.
#' @export
#'
#' @examples
#' ideo <- create_ideo(chr="chr4",ideo.width=400)
#' ideo <- ideo_add_highlight(ideo, region=c(50000,60000))
#'
#'
#' @include AllClass.R
ideo_add_highlight <- function(ideo, region, tick_len = 15) {
    # region is a list
    if (!(class(ideo) == "ideo"))
        stop("First parameter is not ideo object")

    if (!(region[1] < region[2]))
        stop("region is expected to be left_loc < right_loc")

    if(ideo@ydrift >= 0){
        y_min <- ideo@ydrift - tick_len/2 + ideo@height/2
        y_max <- ideo@ydrift + tick_len/2 + ideo@height/2
    }
    else{
        y_min <- ideo@ydrift - tick_len/2 - ideo@height/2
        y_max <- ideo@ydrift + tick_len/2 - ideo@height/2
    }

    df = data.frame(xmin = region[1]/ideo@.plot_scale + ideo@xdrift, xmax = region[2]/ideo@.plot_scale +
        ideo@xdrift, ymin = y_min, ymax = y_max)

    g <- geom_rect(data = df, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), alpha = 0.5,
        fill = "red", color = "red")


    g1 <- c(g, list(theme(axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks = element_blank(),
        axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        legend.position = "none", panel.background = element_rect(fill = "white"), panel.grid.minor = element_line(colour = NA),
        panel.grid.major = element_line(colour = NA))))


    ideo@geom_tick <- g1
    ideo@.tick_bot <- y_min
    ideo@.tick_top <- y_max
    ideo@.leftregion <- df$xmin
    ideo@.rightregion <- df$xmax

    return(ideo)
}

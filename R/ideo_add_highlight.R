# ideo highlight add a highlight region on an ideaogram geom_list,height,xdrift,ydrift

#' Title
#'
#' @param ideo
#' @param region
#' @param tick_len
#'
#' @return
#' @export
#'
#' @examples
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

    # tick_bot <- ideo@ydrift - tick_len/2 + ideo@height/2
    # tick_top <- ideo@ydrift + tick_len/2 + ideo@height/2

    ideo@geom_tick <- g1
    ideo@.tick_bot <- y_min
    ideo@.tick_top <- y_max
    ideo@.leftregion <- df$xmin
    ideo@.rightregion <- df$xmax

    return(ideo)
}

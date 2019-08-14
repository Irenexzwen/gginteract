#' Title Ceate single ideogram ggplot object.
#'
#' @slot geom_ideobody list. ggplot layer list for ideogram presentation.
#' @slot geom_tick list. ggplot layer list for hightlight region of ideogram.
#' @slot chr character. chromosme name of ideogram, starts with "chr".
#' @slot height numeric. single ideogram height (horizontal layout).
#' @slot xdrift numeric. x axis drift from 0.
#' @slot ydrift numeric. y axis drift from 0.
#' @slot chr_start numeric. chromosome start based on annotaion.
#' @slot chr_end   numeric. chromosome end based on annotaion.
#' @slot .plot_scale numeric. scale compression level from real annotation to ggplot presented object.
#' @slot .ideo_left numeric. leftmost coordinate of ideogram object (horizontal layout).
#' @slot .ideo_right numeric. rightmost coordinate of ideogram object (horizontal layout).
#' @slot .tick_top numeric. topmost coordinate of ideogram object (horizontal layout, if highlight region/tick added).
#' @slot .tick_bot numeric. bottommost coordinate of ideogram object (horizontal layout, if highlight region/tick added).
#' @slot .leftregion numeric. start of highlight region (orginal coordinate).
#' @slot .rightregion numeric. end of highlight region (orginal coordinate).
#'
#' @return ideo object with ggplot layer and meta information about the ideogram,
#' the highlight region and the plot object.
#' @export
#'
#' @examples NULL
setClass(Class="ideo",
         representation(
           geom_ideobody = "list",
           geom_tick = "list",
           chr = "character",
           height = "numeric",
           xdrift = "numeric",
           ydrift = "numeric",
           chr_start = "numeric",
           chr_end = "numeric",
           .plot_scale = "numeric",
           .ideo_left = "numeric",
           .ideo_right = "numeric",
           .tick_top = "numeric",
           .tick_bot = "numeric",
           .leftregion = "numeric",
           .rightregion = "numeric"
         )
)


#' Title Gene trascript annotation class
#'
#' @slot name character. gene name.
#' @slot chr_num character. the number of chromosme. eg("1", "5", "X").
#' @slot chr character. chromosme name.
#' @slot chromstart numeric. numeric. chromosome start based on annotaion.
#' @slot chromend numeric. numeric. chromosome end based on annotaion.
#' @slot genelen numeric. length of transcript.
#' @slot center numeric. center coordinate of the transcript.
#' @slot anno data.frame. A datafram contain all exon coordinates.
#'
#' @return gene_anno class.
#' @export
#'
#' @examples NULL
setClass(Class = "gene_anno",
         representation(
           name = "character",
           chr_num = "character",
           chr = "character",
           chromstart = "numeric",
           chromend = "numeric",
           genelen = "numeric",
           center = "numeric",
           anno = "data.frame"
           )
)


#' Title Parallel fashion interaction plot
#' A plot skeleton. Two interacted genes lie in up and bottom in a parallel fashion,
#' where interaction identified from pairend sequence will be presented as polygon connected regions.
#'
#' @slot geom_para list. geom layer list.
#' @slot genetop gene_anno. gene annotation object for the top gene.
#' @slot genebot gene_anno. gene annotation object for the bottom gene.
#' @slot .TopLeft_x  numeric. top left coordinates of the ggplot object. Here top left also indicates the left coord of the gene top.
#' @slot .TopRight_x numeric. top right coordinates of the ggplot object. Here top right also indicates the right coord of the gene top.
#' @slot .BotLeft_x  numeric. bottom left coordinates of the ggplot object. Here botleft also indicates the left coord of the bottom gene.
#' @slot .BotRight_x numeric. bottom right coordinates of the ggplot object. Here botright also indicates the right coord of the bottom gene.
#' @slot .Top_y numeric. top coordinates of the ggplot object.
#' @slot .Bot_y numeric. bottom coordinates of the ggplot object.
#' @slot .Topcolor character. top gene color.
#' @slot .Botcolor character. bottom gene color.
#' @slot .VEXON numeric. exon veriticle height.
#' @slot .VGAP numeric. verticle distance between genetop and gene bottom.
#'
#' @return para obj.
#' @export
#'
#' @examples Null
setClass(Class="para",
         representation(
           geom_para = "list",
           genetop = "gene_anno",
           genebot = "gene_anno",
           .TopLeft_x = "numeric",
           .TopRight_x = "numeric",
           .BotLeft_x = "numeric",
           .BotRight_x = "numeric",
           .Top_y = "numeric",
           .Bot_y = "numeric",
           .Topcolor = "character",
           .Botcolor = "character",
           .VEXON = "numeric",
           .VGAP = "numeric"
         )
)


#' Title Pairend fashion interaction plot
#'
#' @slot geom_pair list. geom layer list.
#' @slot geneleft gene_anno.  gene_anno. gene annotation object for the left gene.
#' @slot generight gene_anno. gene_anno. gene annotation object for the right gene.
#' @slot .TopLeft_x numeric.  top left coordinates of the ggplot object. Here top left also indicates the left coord of the gene left
#' @slot .TopRight_x numeric. top right coordinates of the ggplot object. Here top right also indicates the right coord of the gene right.
#' @slot .BotLeft_x numeric.  bottom left coordinates of the ggplot object. Here botleft also indicates the leftmost coord of the reads mapping to leftgene.
#' @slot .BotRight_x numeric. bottom right coordinates of the ggplot object. Here botright also indicates the rightmost coord of the reads mapping to rightgene.
#' @slot .Top_y numeric. top coordinates of the ggplot object.
#' @slot .Bot_y numeric. bottom coordinates of the ggplot object.
#' @slot .gr1_right numeric. rightmost coordinate of geneleft.
#' @slot .gr2_left numeric. leftmost coordinate of generight.
#' @slot .VEXON numeric. vertical height of exon.
#'
#' @return pairend object
#' @export
#'
#' @examples Null
#'
setClass(Class="pairend",
         representation(
           geom_pair = "list",
           geneleft = "gene_anno",
           generight = "gene_anno",
           .TopLeft_x = "numeric",
           .TopRight_x = "numeric",
           .BotLeft_x = "numeric",
           .BotRight_x = "numeric",
           .Top_y = "numeric",
           .Bot_y = "numeric",
           .gr1_right = "numeric",
           .gr2_left = "numeric",
           .VEXON = "numeric"
         )
)

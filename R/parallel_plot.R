#' Title Parallel interaction plot
#' One-station function for parallel interaction plot
#'
#' @param GENE1_anno Dataframe, a bed file format dataframe for gene1. each row represents an exon with chr, start, end.
#' @param GENE2_anno Dataframe, a bed file format dataframe for gene2. each row represents an exon with chr, start, end.
#' @param R1 Dataframe, a bed file format dataframe for Read1 (first end) bed file.
#' @param R2 Dataframe, a bed file format dataframe for Read2 (second end) bed file.
#' @param GENE1_COLOR character, In R, colors can be specified either by name e.g col = 'red' or as a hexadecimal RGB triplet
#' @param GENE2_COLOR character.In R, colors can be specified either by name e.g col = 'red' or as a hexadecimal RGB triplet
#'
#' @return ggplot object
#' @export
#'
#' @examples
#' data(GENE1_anno,GENE2_anno,R1,R2)
#' para <- parallel_plot(GENE1_anno,GENE2_anno,R1,R2)
#'
parallel_plot <- function(GENE1_anno,GENE2_anno,R1,R2,GENE1_COLOR="#deb210",GENE2_COLOR="#668ed1",genome="hg38",genename1="",genename2=""){

  k <- parallel_inter(GENE1_anno,GENE2_anno,R1,R2,GENE1_COLOR="#deb210",GENE2_COLOR="#668ed1",genename1 = genename1,genename2 = genename2)

  ideo_width <- k@.BotRight_x/2.5
  xd <- k@.BotRight_x/2 - ideo_width/2

  # create ideogram based on parallel interaction skeleton ----
  upideo <- create_ideo(genome=genome,
                       k@genetop@chr,
                       ideo.width=ideo_width,
                       ydrift=k@.Top_y + k@.VEXON,
                       xdrift=xd)

  botideo <- create_ideo(genome=genome,
                        k@genebot@chr,
                        ideo.width=ideo_width,
                        ydrift=k@.Bot_y - k@.VEXON,
                        xdrift=xd)

  # add high light region on -----------------------------------
  upideo <- ideo_add_highlight(ideo = upideo,
                               region = c(k@genetop@chromstart, k@genetop@chromend),
                               tick_len = ideo_width/16*1.6
                               )

  botideo <- ideo_add_highlight(ideo = botideo,
                               region = c(k@genebot@chromstart, k@genebot@chromend),
                               tick_len = ideo_width/16*1.6
  )

  # add-dash line -----------------------------------------------
  dash <- data.frame(x = c(upideo@.leftregion,
                           k@.TopLeft_x,
                           upideo@.rightregion,
                           k@.TopRight_x,
                           botideo@.leftregion,
                           k@.BotLeft_x,
                           botideo@.rightregion,
                           k@.BotRight_x),
                     y = c(upideo@.tick_bot,
                           k@.Top_y,
                           upideo@.tick_bot,
                           k@.Top_y,
                           botideo@.tick_top,
                           k@.Bot_y,
                           botideo@.tick_top,
                           k@.Bot_y))

  dash['group'] <- rep(seq(1,4),each=2)

  # plot transcript text ---------------------------------------------------
  nudge <- k@.VEXON
  TextDF <- data.frame(x=c(k@.TopLeft_x,k@.BotLeft_x),
                       y=c(k@.Top_y-k@.VEXON,k@.Bot_y+k@.VEXON),
                       label=c(k@genetop@name,k@genebot@name))


  # plot ideogram text -----------------------------------------------------
  format_ <- function(x){
    return(formatC(x,format = "f", big.mark = ",",drop0trailing = T))
  }
  text_top <- paste0(k@genetop@chr,":",format_(k@genetop@chromstart),"-",format_(k@genetop@chromend))
  text_bot <- paste0(k@genebot@chr,":",format_(k@genebot@chromstart),"-",format_(k@genebot@chromend))
  center <- (k@.TopLeft_x+k@.TopRight_x)/2
  ideotext <- data.frame(x=c(center,center),
                         y=c(upideo@.tick_top + nudge, botideo@.tick_bot - nudge),
                         labels = c(text_top,text_bot))

  # plot all ----------------------------------------------------
  ggplot() + k@geom_para + upideo@geom_ideobody +
    botideo@geom_ideobody +upideo@geom_tick + botideo@geom_tick +
    geom_line(data = dash, aes(x=x, y=y, group=group), linetype="dashed",color = "red")+
    geom_text(data=TextDF,aes(x=x-nudge*3,y=y),label=TextDF$label)+
    geom_text(data=ideotext, aes(x=x, y=y), label=ideotext$labels)


}






#' Title plot parallel interaction plot skeleton
#'
#' @param GENE1_anno Dataframe, gene annotaion bed file. eg could be seen \code{data(GENE1_anno)}
#' @param GENE2_anno Dataframe, gene annotaion bed file. eg could be seen \code{data(GENE2_anno)}
#' @param R1 Dataframe, Read1 annotation file. eg could be loaded with \code{data(R1)}
#' @param R2 Dataframe, Read2 annotation file. eg could be loaded with \code{data(R2)}
#' @param GENE1_COLOR String, with default "#deb210"
#' @param GENE2_COLOR string, with default "#668ed1"
#' @param xdrift numeric. x axis drift from 0.
#' @param ydrift numeric. y axis drift from 0.
#'
#' @return para object
#' @export
#'
#' @examples
#' data(GENE1_anno,GENE2_anno,R1,R2)
#' para <- parallel_inter(GENE1_anno,GENE2_anno,R1,R2)
#' ggplot() + para@geom_para
parallel_inter <- function(GENE1_anno,GENE2_anno,R1,R2,GENE1_COLOR="#deb210",GENE2_COLOR="#668ed1",xdrift=0,ydrift=0,genename1="",genename2=""){


  # reorganize R1 and R2 pair
  R1_bed <- R1 %>% dplyr::arrange(V4)
  R2_bed <- R2 %>% dplyr::arrange(V4)

  ppi <- rbind(R1_bed,R2_bed) %>% dplyr::arrange(V4) %>% as.data.frame()
  colnames(ppi) <- c("chr","xstart",'xend','readsname','height','strand') # need refine

  #filter unpaired reads
  filter_unpair <- function(ppi){
    count <- ppi %>% dplyr::group_by(readsname) %>% dplyr::summarize(n=length(readsname))
    test <- ppi[ppi$readsname %in% count$readsname[count$n==2], ]
    return(test)
  }

  ppi <- filter_unpair(ppi)

  # reform chr name accepts both NC* or chr*
  ppi$chr <- ppi$chr %>% gsub("^NC_0{4}0?","",.) %>% gsub("\\..*$","",.) %>% gsub("[Cc]hr","",.)
  ppi['pair'] <- factor(rep(seq(1,nrow(ppi)/2),each=2))

  # genrate gene_anno class
  GENE1_anno <-  gene_anno(GENE1_anno,genename1)
  GENE2_anno <-  gene_anno(GENE2_anno,genename2)

  VGAP = max(GENE1_anno@genelen,GENE2_anno@genelen)/3
  VEXON <- VGAP/20

  # distance between two transcript centers
  # adjust horizontal coordinates
  MOVE = GENE1_anno@center-GENE2_anno@center

  for (i in seq(1:nrow(ppi))){
    if(ppi[i,1]==GENE1_anno@chr_num){
      ppi[i,'xstart'] <- ppi[i,'xstart'] - GENE1_anno@chromstart
      ppi[i,'xend'] <- ppi[i,'xend'] - GENE1_anno@chromstart
    }
    else {
      ppi[i,'xstart'] <- ppi[i,'xstart'] - GENE2_anno@chromstart + MOVE  # move to center alignment
      ppi[i,'xend'] <- ppi[i,'xend'] - GENE2_anno@chromstart + MOVE
    }
  }

  #### draw the up and bottom transcript anno
  # GENE1 at the bottom
  # GENE2 at the top

  GENE1_anno@anno['yvalue'] <- 0
  GENE1_anno@anno['height'] <- VEXON
  GENE2_anno@anno['yvalue'] <- 0
  GENE2_anno@anno['height'] <- VEXON

  gr1 <- GENE1_anno@anno
  gr2 <- GENE2_anno@anno

  # horizon adjust
  gr1$start <- gr1$start-GENE1_anno@chromstart
  gr1$end <- gr1$end-GENE1_anno@chromstart
  gr2$start <- gr2$start-GENE2_anno@chromstart+MOVE
  gr2$end <- gr2$end-GENE2_anno@chromstart+MOVE

  # verticle adjust
  gr2$yvalue <- VGAP


  ################# draw ppi in polygon ####################
  # generate polygon dataframe
  # x axis
  poly_df_x <- ppi[,c(2,3)] %>% t() %>%
    as.matrix() %>% matrix(.,nrow=nrow(ppi)/2,byrow = T) %>%
    t() %>% as.data.frame() %>% reshape2::melt()

  colnames(poly_df_x) <- c("chr","x")
  poly_df_x$chr <- rep(ppi$chr,each=2)
  poly_df_x['pair'] <- rep(seq(1,nrow(ppi)/2),each=4)

  for(i in seq(1,nrow(poly_df_x))){ # add y value according to chr
    if(poly_df_x[i,1] == GENE1_anno@chr_num){
      poly_df_x[i,'y'] <- 0+VEXON/2+VEXON/10
    }
    else{
      poly_df_x[i,'y'] <- (VGAP - VEXON/2 - VEXON/10)
    }
  }

  for(i in seq(1,nrow(poly_df_x),by=4)){
    row3 <- poly_df_x$x[i+2]
    row4 <- poly_df_x$x[i+3]
    poly_df_x$x[i+2] <- row4
    poly_df_x$x[i+3] <- row3
  }

  # adjust if drift is assign
  poly_df_x$x <- poly_df_x$x + xdrift
  poly_df_x$y <- poly_df_x$y + ydrift
  gr1$start <- gr1$start + xdrift
  gr1$end <- gr1$end + xdrift
  gr2$start <- gr2$start + xdrift
  gr2$end <- gr2$end + xdrift
  gr1$yvalue <- gr1$yvalue + ydrift
  gr2$yvalue <- gr2$yvalue + ydrift


  # return ggplot layer list
    p <- list(geom_polygon(data=poly_df_x,aes(x=x,y=y,group=pair),fill="darkgrey",alpha=0.1),
    geom_rect(data=gr1,aes(xmin=start,ymin=yvalue-height/2,xmax=end,ymax=yvalue+height/2),color=GENE1_COLOR,fill=GENE1_COLOR), #annotation
    geom_line(data=gr1,aes(start,y=yvalue),size=1,color=GENE1_COLOR),
    geom_rect(data=gr2,aes(xmin=start,ymin=yvalue-height/2,xmax=end,ymax=yvalue+height/2),color=GENE2_COLOR,fill=GENE2_COLOR), #annotation
    geom_line(data=gr2,aes(start,y=yvalue),size=1,color=GENE2_COLOR),
    theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),legend.position="none",
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank())
  )

  return(new("para",
             geom_para = p,
             genetop = GENE2_anno,
             genebot = GENE1_anno,
             .TopLeft_x = min(gr2$start),
             .TopRight_x = max(gr2$end),
             .BotLeft_x = min(gr1$start),
             .BotRight_x = max(gr1$end),
             .Top_y = gr2$yvalue[1] + VEXON,
             .Bot_y = gr1$yvalue[1] - VEXON,
             .Topcolor = GENE2_COLOR,
             .Botcolor = GENE1_COLOR,
             .VEXON = VEXON,
             .VGAP = VGAP
  ))
}





#' Title Pairend interaction plot
#'
#' @param GENE1_anno Dataframe, a bed file format dataframe for gene1. each row represents an exon with chr, start, end.
#' @param GENE2_anno Dataframe, a bed file format dataframe for gene2. each row represents an exon with chr, start, end.
#' @param R1 Dataframe, a bed file format dataframe for Read1 (first end) bed file.
#' @param R2 Dataframe, a bed file format dataframe for Read2 (second end) bed file.
#' @param GENE1_COLOR character, color of left gene. In R, colors can be specified either by name e.g col = "red" or as a hexadecimal RGB triplet
#' @param GENE2_COLOR character, color of right gene. In R, colors can be specified either by name e.g col = "red" or as a hexadecimal RGB triplet
#' @param xdrift numeric. x axis drift from 0.
#' @param ydrift numeric. y axis drift from 0.
#' @param VEXON numeric. verticle height of exon.
#' @param genome String. Genome version eg."hg38","mm10","mm20". default "hg38"
#'
#' @return ggplot object of pairend interaction plot
#' @export pairend_plot
#'
#' @examples
#' data(GENE1_anno,GENE2_anno,R1,R2)
#' pairend <- pairend_plot(GENE1_anno,GENE2_anno,R1,R2)
#'
pairend_plot <- function(GENE1_anno,GENE2_anno,R1,R2,GENE1_COLOR="#deb210",GENE2_COLOR="#668ed1",genename1 = "", genename2 = "", xdrift=0,ydrift=0,VEXON=10,genome="hg38"){


  k <- pairend_inter(GENE1_anno,GENE2_anno,R1,R2,GENE1_COLOR="#deb210",GENE2_COLOR="#668ed1",xdrift=0,ydrift=0,VEXON=10, genename1 = genename1, genename2 = genename2)

  ideo_width <- (k@.TopRight_x - k@.TopLeft_x)/3
  yd <- k@.Top_y + k@.VEXON # ydrift
  width_ratio <- ideo_width/k@.VEXON

  # ideogram
  leftideo <- create_ideo(genome=genome,
                         k@geneleft@chr,
                         ideo.width=ideo_width,
                         ydrift=yd,
                         xdrift=k@.TopLeft_x,
                         whratio = width_ratio)

    rightideo <- create_ideo(genome=genome,
                         k@generight@chr,
                         ideo.width=ideo_width,
                         ydrift=yd,
                         xdrift=k@.TopRight_x - ideo_width,
                         whratio = width_ratio)

  # add tick region ------------------------------------------------------
  leftideo <- ideo_add_highlight(ideo = leftideo,
                               region = c(k@geneleft@chromstart, k@geneleft@chromend),
                               tick_len = k@.VEXON*1.6
  )

  rightideo <- ideo_add_highlight(ideo = rightideo,
                                region = c(k@generight@chromstart, k@generight@chromend),
                                tick_len = k@.VEXON*1.6
  )

  # add-dash line -----------------------------------------------
  dash <- data.frame(x = c(leftideo@.leftregion,
                           k@.TopLeft_x,
                           leftideo@.rightregion,
                           k@.gr1_right,
                           rightideo@.leftregion,
                           k@.gr2_left,
                           rightideo@.rightregion,
                           k@.TopRight_x),
                     y = c(leftideo@.tick_bot,
                           k@.Top_y,
                           leftideo@.tick_bot,
                           k@.Top_y,
                           rightideo@.tick_bot,
                           k@.Top_y,
                           rightideo@.tick_bot,
                           k@.Top_y))

  dash['group'] <- rep(seq(1,4),each=2)


  # plot ideogram text -----------------------------------------------------
  format_ <- function(x){
    return(formatC(x,format = "f", big.mark = ",",drop0trailing = T))
  }
  text_right <- paste0(k@geneleft@chr,":",format_(k@geneleft@chromstart),"-",format_(k@geneleft@chromend))
  text_left <- paste0(k@generight@chr,":",format_(k@generight@chromstart),"-",format_(k@generight@chromend))
  left_center <- (leftideo@.ideo_left+leftideo@.ideo_right)/2
  right_center <- (rightideo@.ideo_left+rightideo@.ideo_right)/2
  ideotext <- data.frame(x=c(left_center,right_center),
                         y=c(rightideo@.tick_top,leftideo@.tick_top),
                         labels = c(text_left,text_right))

  nudge <- k@.VEXON
  # plot all ----------------------------------------------------

  ggplot() +k@geom_pair+ leftideo@geom_ideobody + rightideo@geom_ideobody+
    rightideo@geom_tick + leftideo@geom_tick+
    geom_line(data = dash, aes(x=x, y=y, group=group), linetype="dashed",color = "red")+
    geom_text(data=ideotext, aes(x=x, y=y+nudge), label=ideotext$labels)
  }





# pair-end interaction skeleton--------------------------

#' Title pair-end inter skeleton
#'
#' @param GENE1_anno Dataframe, a bed file format dataframe for gene1. each row represents an exon with chr, start, end.
#' @param GENE2_anno Dataframe, a bed file format dataframe for gene2. each row represents an exon with chr, start, end.
#' @param R1 Dataframe, a bed file format dataframe for Read1 (first end) bed file.
#' @param R2 Dataframe, a bed file format dataframe for Read2 (second end) bed file.
#' @param GENE1_COLOR String, with default "#deb210".
#' @param GENE2_COLOR string, with default "#668ed1".
#' @param xdrift numeric. x axis drift from 0.
#' @param ydrift numeric. y axis drift from 0.
#' @param VEXON numeric. verticle height of exon.
#'
#' @return pair-end inter skeleton, ggplot object
#' @export
#'
#' @examples
#' data(GENE1_anno,GENE2_anno,R1,R2)
#' pairend_skeleton <- pairend_inter(GENE1_anno,GENE2_anno,R1,R2)
#' ggplot()+pairend_skeleton@geom_pair

pairend_inter <- function(GENE1_anno,GENE2_anno,R1,R2,GENE1_COLOR="#deb210",GENE2_COLOR="#668ed1",xdrift=0,ydrift=0,VEXON=10,genename1="",genename2=""){


  # reorganize R1 and R2 pair
  R1_bed <- R1 %>% dplyr::arrange(V4)
  R2_bed <- R2 %>% dplyr::arrange(V4)

  ppi <- rbind(R1_bed,R2_bed) %>% dplyr::arrange(V4) %>% as.data.frame()
  colnames(ppi) <- c("chr","xstart",'xend','readsname','height','strand') # need refine

  #filter unpaired reads
  filter_unpair <- function(ppi){
    count <- ppi %>% dplyr::group_by(readsname) %>% dplyr::summarise(n=n())
    test <- ppi[ppi$readsname %in% count$readsname[count$n==2], ]
    return(test)
  }

  ppi <- filter_unpair(ppi)

  # reform chr name accepts both NC* or chr*
  ppi$chr <- ppi$chr %>% gsub("^NC_0{4}0?","",.) %>% gsub("\\..*$","",.) %>% gsub("[Cc]hr","",.)
  ppi['pair'] <- factor(rep(seq(1,nrow(ppi)/2),each=2))
  ppi['yvalue'] <- rep(seq(1,nrow(ppi)/2)*VEXON,each=2)

  # genrate gene_anno class
  GENE1_anno <-  gene_anno(GENE1_anno,genename1)
  GENE2_anno <-  gene_anno(GENE2_anno,genename2)

  # set meta data
  HGAP = sum(GENE1_anno@genelen,GENE2_anno@genelen)*1.2
  HEIGHT=VEXON

  # READ1 lands in left and READ2 in right
  # adjust horizontal coordinates
  MOVE = HGAP

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

  # split ppi into read1 DF and read2 DF to avoid color fill conflict with ideogram
  ppi_left <- ppi[ppi$chr == GENE1_anno@chr_num,]
  ppi_right <- ppi[ppi$chr == GENE2_anno@chr_num,]

  # draw gene anno track
  # GENE1 on the left and GENE2 on the right
  GENE1_anno@anno['yvalue'] <- max(ppi$yvalue)+HEIGHT*3
  GENE1_anno@anno['height'] <- HEIGHT
  GENE2_anno@anno['yvalue'] <- max(ppi$yvalue)+HEIGHT*3
  GENE2_anno@anno['height'] <- HEIGHT

  gr1 <- GENE1_anno@anno
  gr2 <- GENE2_anno@anno

  # horizon adjust
  gr1$start <- gr1$start - GENE1_anno@chromstart
  gr1$end <- gr1$end - GENE1_anno@chromstart
  gr2$start <- gr2$start - GENE2_anno@chromstart + MOVE
  gr2$end <- gr2$end - GENE2_anno@chromstart + MOVE

  # gene name text --------------------------------------------------------
  nudge <- (min(gr2$start)-max(gr1$end))/6 # need refine
  TextDF <- data.frame(x=c(max(gr1$end)+nudge,min(gr2$start)-nudge),
                       y=c(gr1$yvalue[1],gr2$yvalue[1]),
                       label=c(GENE1_anno@name,GENE2_anno@name))


  p <- list(
    geom_line(data=ppi,aes(xstart,y=yvalue+HEIGHT/2,group=pair),size=0.1,color="darkgrey"),
    geom_rect(data=ppi_left,aes(xmin=xstart,ymin=yvalue,xmax=xend,ymax=yvalue+HEIGHT),color=GENE1_COLOR,fill=GENE1_COLOR),
    geom_rect(data=ppi_right,aes(xmin=xstart,ymin=yvalue,xmax=xend,ymax=yvalue+HEIGHT),color=GENE2_COLOR,fill=GENE2_COLOR),
    geom_line(data=gr1,aes(start,y=yvalue),size=1,color="darkgrey"),
    geom_rect(data=gr1,aes(xmin=start,ymin=yvalue-HEIGHT/2,xmax=end,ymax=yvalue+HEIGHT/2),color="darkgrey",fill="darkgrey"),
    geom_line(data=gr2,aes(start,y=yvalue),size=1,color="darkgrey"),
    geom_rect(data=gr2,aes(xmin=start,ymin=yvalue-HEIGHT/2,xmax=end,ymax=yvalue+HEIGHT/2),color="darkgrey",fill="darkgrey"),
    geom_text(data=TextDF,aes(x=x,y=y),label=TextDF$label,color="darkgrey"),
    theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),legend.position="none",
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank())

  )

  # return --------------------------------------------

  return(new("pairend",
             geom_pair = p,
             geneleft =  GENE1_anno,
             generight = GENE2_anno,
             .TopLeft_x = min(gr1$start),
             .TopRight_x = max(gr2$end),
             .BotLeft_x = min(ppi$xstart[ppi$chr == GENE1_anno@chr_num]),
             .BotRight_x = max(ppi$xend[ppi$chr == GENE2_anno@chr_num]),
             .Top_y = gr1$yvalue[1]+HEIGHT/2,
             .Bot_y = ppi$yvalue[1]-HEIGHT/2,
             .gr1_right = max(gr1$end),
             .gr2_left = min(gr2$start),
             .VEXON = VEXON))

}

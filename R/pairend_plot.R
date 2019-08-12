pairend_plot <- function(GENE1_anno,GENE2_anno,R1,R2,GENE1_COLOR="#deb210",GENE2_COLOR="#668ed1",xdrift=0,ydrift=0,VEXON=10,genome="hg38"){

  k <- pairend_inter(GENE1_anno,GENE2_anno,R1,R2,GENE1_COLOR="#deb210",GENE2_COLOR="#668ed1",xdrift=0,ydrift=0,VEXON=10)

  ideo_width <- (k@.TopRight_x - k@.TopLeft_x)/3
  yd <- k@.Top_y + k@.VEXON*3 # ydrift
  width_ratio <- ideo_width/k@.VEXON

  genome_ <- genome
  # ideogram
  leftideo <- creat_ideo(genome="hg38",
                         k@geneleft@chr,
                         ideo.width=ideo_width,
                         ydrift=yd,
                         xdrift=k@.TopLeft_x,
                         whratio = width_ratio)

  rightideo <- creat_ideo(genome="hg38",
                         k@generight@chr,
                         ideo.width=ideo_width,
                         ydrift=yd,
                         xdrift=k@.TopRight_x - ideo_width,
                         whratio = width_ratio)


  ggplot() + leftideo@geom_ideobody + rightideo@geom_ideobody+ k@geom_pair
  }





# pair-end inter skeleton--------------------------

pairend_inter <- function(GENE1_anno,GENE2_anno,R1,R2,GENE1_COLOR="#deb210",GENE2_COLOR="#668ed1",xdrift=0,ydrift=0,VEXON=10){


  # reorganize R1 and R2 pair
  R1_bed <- read.table(R1,header = F,stringsAsFactors = F) %>% dplyr::arrange(V4)
  R2_bed <- read.table(R2,header = F,stringsAsFactors = F) %>% dplyr::arrange(V4)

  ppi <- rbind(R1_bed,R2_bed) %>% dplyr::arrange(V4) %>% as.data.frame()
  colnames(ppi) <- c("chr","xstart",'xend','readsname','height','strand') # need refine

  #filter unpaired reads
  filter_unpair <- function(ppi){
    count <- ppi %>% group_by(readsname) %>% dplyr::summarise(n=n())
    test <- ppi[ppi$readsname %in% count$readsname[count$n==2], ]
    return(test)
  }

  ppi <- filter_unpair(ppi)

  # reform chr name accepts both NC* or chr*
  ppi$chr <- ppi$chr %>% gsub("^NC_0{4}0?","",.) %>% gsub("\\..*$","",.) %>% gsub("[Cc]hr","",.)
  ppi['pair'] <- factor(rep(seq(1,nrow(ppi)/2),each=2))
  ppi['yvalue'] <- rep(seq(1,nrow(ppi)/2)*VEXON,each=2)

  # genrate gene_anno class
  GENE1_anno <-  read_anno(GENE1_anno)
  GENE2_anno <-  read_anno(GENE2_anno)

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

  # reads fill color ------------------------------------------------------
  col <- ifelse(ppi$chr == GENE1_anno@chr, GENE1_COLOR, GENE2_COLOR)

  p <- list(
    geom_line(data=ppi,aes(xstart,y=yvalue+HEIGHT/2,group=pair),size=0.1,color="darkgrey"),
    geom_rect(data=ppi,aes(xmin=xstart,ymin=yvalue,xmax=xend,ymax=yvalue+HEIGHT,fill=col)),
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
             geneleft =  GENE2_anno,
             generight = GENE1_anno,
             .TopLeft_x = min(gr1$start),
             .TopRight_x = max(gr2$end),
             .BotLeft_x = min(ppi$xstart[ppi$chr == GENE1_anno@chr_num]),
             .BotRight_x = max(ppi$xend[ppi$chr == GENE2_anno@chr_num]),
             .Top_y = gr1$yvalue[1]+HEIGHT/2,
             .Bot_y = ppi$yvalue[1]-HEIGHT/2,
             .VEXON = VEXON))

}

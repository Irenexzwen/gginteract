# Sankey transcript interaction style


setwd(WORK_PATH)

GENE1_COLOR <- "#deb210"
GENE2_COLOR <- "#668ed1"
VGAP <- 200  # verticle gap between two transcript
VEXON <- VGAP/20  # transcript exon height

### deal with ppi pair info
R1_bed <- read.table(READ1,header = F,stringsAsFactors = F) %>% arrange(V4)
R2_bed <- read.table(READ2,header = F,stringsAsFactors = F) %>% arrange(V4)
ppi <- rbind(R1_bed,R2_bed) %>% arrange(V4) %>% as.data.frame()
colnames(ppi) <- c("chr","xstart",'xend','readsname','height','strand')

filter_unpair <- function(ppi){
  count <- ppi %>% group_by(readsname) %>% dplyr::summarise(n=n())
  test <- ppi[ppi$readsname %in% count$readsname[count$n==2], ]
  return(test)
}

ppi <- filter_unpair(ppi)
ppi$chr <- ppi$chr %>% gsub("^NC_0{4}0?","",.) %>% gsub("\\..*$","",.)
ppi['pair'] <- factor(rep(seq(1,nrow(ppi)/2),each=2))


read_anno <- function(file_name){
  gene_anno1 <- read.table(file_name,header = F,stringsAsFactors = F)
  anno <- data.frame(x='hi')
  anno['name'] <-  gsub("\\..*$","",file_name)
  anno['x'] <- gene_anno1$V1[1] %>% gsub("chr","",.)
  anno['chr'] <-  gene_anno1$V1[1]  # chr name with chrID
  anno['chromstart'] <-  min(gene_anno1$V2) # from IGV bed start MIN
  anno['chromend'] <- max(gene_anno1$V3)
  anno['genelen'] <- max(gene_anno1$V3) - min(gene_anno1$V2)
  anno['center'] <- (max(gene_anno1$V3) - min(gene_anno1$V2))/2
  return(anno)
}


n1 <-  read_anno(GENE1)
n2 <-  read_anno(GENE2)

MOVE = n1$center-n2$center # distance between two centers

for (i in seq(1:nrow(ppi))){
  if(ppi[i,1]==n1$x){
    ppi[i,'xstart'] <- ppi[i,'xstart']-n1$chromstart
    ppi[i,'xend'] <- ppi[i,'xend']-n1$chromstart
  }
  else {
    ppi[i,'xstart'] <- ppi[i,'xstart']-n2$chromstart+MOVE  # move to center alignment
    ppi[i,'xend'] <- ppi[i,'xend']-n2$chromstart+MOVE
  }
}



#### draw the up and bottom transcript anno
# gene1 in the bottom
# gene2 at the top

read_trans <- function(GENE1){
  t <- read.table(GENE1,header = F,stringsAsFactors = F)
  colnames(t) <- c("chr","start","end","exon","strand")
  t['yvalue'] <- 0
  t['height'] <- VEXON
  return(t)
}

gr1 <- read_trans(GENE1)
gr2 <- read_trans(GENE2)

# horizon adjust
gr1$start <- gr1$start-n1$chromstart
gr1$end <- gr1$end-n1$chromstart
gr2$start <- gr2$start-n2$chromstart+MOVE
gr2$end <- gr2$end-n2$chromstart+MOVE

# verticle adjust
gr2$yvalue <- VGAP

################# draw ppi in polygon ####################

# generate polygon df
# x axis
poly_df_x <- ppi[,c(2,3)] %>% t() %>%
  as.matrix() %>% matrix(.,nrow=nrow(ppi)/2,byrow = T) %>%
  t() %>% as.data.frame() %>% reshape2::melt()

colnames(poly_df_x) <- c("chr","x")
poly_df_x$chr <- rep(ppi$chr,each=2)
poly_df_x['pair'] <- rep(seq(1,nrow(ppi)/2),each=4)

for(i in seq(1,nrow(poly_df_x))){ # add y value according to chr
  if(poly_df_x[i,1] == n1$x){
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


ggplot()+geom_polygon(data=poly_df_x,aes(x=x,y=y,group=pair),fill="darkgrey",alpha=0.1)+
  geom_rect(data=gr1,aes(xmin=start,ymin=yvalue-height/2,xmax=end,ymax=yvalue+height/2),color=GENE1_COLOR,fill=GENE1_COLOR)+ #annotation
  geom_line(data=gr1,aes(start,y=yvalue),size=1,color=GENE1_COLOR)+
  geom_rect(data=gr2,aes(xmin=start,ymin=yvalue-height/2,xmax=end,ymax=yvalue+height/2),color=GENE2_COLOR,fill=GENE2_COLOR)+ #annotation
  geom_line(data=gr2,aes(start,y=yvalue),size=1,color=GENE2_COLOR)+
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())






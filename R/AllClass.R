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

# datastructure for one transcript
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


setClass(Class="pairend",
         representation(
           geom_pair = "list",
           geneleft = "gene_anno",
           generight = "gene_anno",
           .TopLeft_x = "numeric",  #gene1 left x
           .TopRight_x = "numeric", #gene2 right x
           .BotLeft_x = "numeric",  #read1 left x
           .BotRight_x = "numeric", #read2 right x
           .Top_y = "numeric",      #gene1/2 top
           .Bot_y = "numeric",      #read1/2 bot
           .gr1_left = "numeric",   #leftgene anno right
           .gr2_right = "numeric",  #rightgene anno left
           .VEXON = "numeric"
         )
)

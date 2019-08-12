read_anno <- function(file_name){
  gene_anno <- data.table::fread(file_name) %>% as.data.frame()
  colnames(gene_anno) <- c("chr","start","end")
  # gene_anno['yvalue'] <- 0       should add in other method of gene_anno
  # gene_anno['height'] <- VEXON

  return(new("gene_anno",
             name = basename(file_name) %>% gsub("\\..*$","",.),
             chr_num = gene_anno$chr[1] %>% gsub("chr","",.),
             chr = gene_anno$chr[1],
             chromstart = min(gene_anno$start),
             chromend = max(gene_anno$end),
             genelen = max(gene_anno$end) - min(gene_anno$start),
             center = (max(gene_anno$end) - min(gene_anno$start))/2,
             anno = gene_anno
             ))
}





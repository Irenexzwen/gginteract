library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(biovizBase)
library(ggsci)
library(export)
# WORK_PATH <- 'G:/zhongsheng/PPI/20180804plot/RPL35_NPM1_plotFiles/RPL35_NPM1_plotFiles'
READ1 <- "Read1.bed"
READ2 <- "Read2.bed"
GENE1 <- "NPM1.igv.bed"
GENE2 <- "RPL35.igv.bed"
# setwd(WORK_PATH)

READ1 <- "G:/zhongsheng/PPI/20180804plot/RPL35_NPM1_plotFiles/RPL35_NPM1_plotFiles/Read1.bed"
READ2 <- "G:/zhongsheng/PPI/20180804plot/RPL35_NPM1_plotFiles/RPL35_NPM1_plotFiles/Read2.bed"
GENE1 <- "G:/zhongsheng/PPI/20180804plot/RPL35_NPM1_plotFiles/RPL35_NPM1_plotFiles/NPM1.igv.bed"
GENE2 <- "G:/zhongsheng/PPI/20180804plot/RPL35_NPM1_plotFiles/RPL35_NPM1_plotFiles/RPL35.igv.bed"

# test pair end
R1 <- "G:/zhongsheng/PPI/073120196_Plot/RPS7_DDX23_plotFiles/RPS7_DDX23_plotFiles/Read1.bed"
R2 <- "G:/zhongsheng/PPI/073120196_Plot/RPS7_DDX23_plotFiles/RPS7_DDX23_plotFiles/Read2.bed"
GENE1_anno <- "G:/zhongsheng/PPI/073120196_Plot/RPS7_DDX23_plotFiles/RPS7_DDX23_plotFiles/DDX23.igv.bed"
GENE2_anno <- "G:/zhongsheng/PPI/073120196_Plot/RPS7_DDX23_plotFiles/RPS7_DDX23_plotFiles/RPS7.igv.bed"



### deal with ppi pair info
k <- parallel_inter(GENE1_anno = GENE1, GENE2_anno = GENE2, R1 = READ1, R2 = READ2)
q <- parallel_plot(GENE1_anno = GENE1, GENE2_anno = GENE2,R1 = READ1, R2 = READ2,GENE1_COLOR="#deb210",GENE2_COLOR="#668ed1",genome="hg38")


# test pair-end
test <- pairend_inter(GENE1_anno = GENE1_anno,GENE2_anno = GENE2_anno,R1 = R1,R2 = R2)
pairend_plot(GENE1_anno = GENE1_anno,GENE2_anno = GENE2_anno,R1 = R1,R2 = R2)

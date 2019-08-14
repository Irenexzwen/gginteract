# gginteract
**gginteract**(Pairend interaction visualization with customizable ggplot object) is an in-house visualization tool for pairend reads interaction visualization.
These "pairend reads" can indirectly suggest the molecules interaction with next generation high throughput RNA-seq method (methods include MARGI, MARIO, PPI etc).
This package is specifically designed to present **two molecules** interaction in two fashions.
![](images/two_style.png)
These fashions trying to capture the R1 and R2 mapping details in a condensed way.
These plots are achieved by aggregating the building blocks of gginteract packages, they are:
- ideograms
- gene annotaion tracks
- parallel pairend skeleton
- pairend horizontal skeleton
- text and marks

## Motivation
Existing solutions ([gviz](https://bioconductor.riken.jp/packages/3.0/bioc/html/Gviz.html), [ggbio](http://bioconductor.org/packages/release/bioc/html/ggbio.html)) 
for merging multifacet genome annotation (eg, ideogram or transcript annotation) using a *track stack* strategy, where each object occupys a track and different tracks stacked vertically share the same horizontal coordinates.
However, these tools fail to offer a flexible way to manipulate details of the tracks. For example it could be difficult to:
* Draw ideograms at random positions while control their width and height.
* Add other customizable information on existing track sharing the same coordinate. (eg. a zoom-in view, a gene name list of a gene dense region highlighted on an ideogram).


## Design philosophy
gginteract build every object on top of ggplot2 buildings to ensure that every building blocks of the plot is a ggplot layer, which means
you could add any other layers on top of the sckeleton exported from gginteract package. 

To ensure users could easily aquire the accurate position of the skeleton, each skeleton object is an S4 object with location details store in their slots.

## Installation
Install ggplot R pakcage from github. 

```
$ R
> library(devtools)
> install_github("irenexzwen/gginteract")
```

## Basic Usage
### 1. Required input data
To generate interaction plot, at least four inputs are required:

first two are gene annotation file, at least first three columns are needed (_chr name, exon start, exon end_). :
- `GENE1_anno` Dataframe, a bed file format dataframe for gene1. each row represents an exon with _chr, start, end_.
- `GENE2_anno` Dataframe, a bed file format dataframe for gene2. each row represents an exon with _chr, start, end_.

Next, reads information in bed format. At least first four columns are needed (_chr name, read start, read end, read name_). : 
- `Read1` Dataframe, a bed file format dataframe for Read1 (first end) bed file.
- `Read2` Dataframe, a bed file format dataframe for Read2 (second end) bed file.

Example data could be checked with:
```R
library(gginteract)
data("R1","R2","GENE1_anno","GENE2_anno")
```
### 2. Parallel interaction plot
To generate a parallel interaction plot, you could simply use one function:
```R
data("R1","R2","GENE1_anno","GENE2_anno")
para <- gginteract::parallel_plot(GENE1_anno = GENE1_anno,
                                  GENE2_anno = GENE2_anno,
                                  R1 = R1,
                                  R2 = R2,
                                  genename1 = "DDX23",
                                  genename2 = "RPS7")
para
```
![](images/para.png)

### 3. Pairend interaction plot
To generate a parallel interaction plot, you could simply use one function:
```R
data("R1","R2","GENE1_anno","GENE2_anno")
pair <- gginteract::pairend_plot(GENE1_anno = GENE1_anno,
                                 GENE2_anno = GENE2_anno,
                                 R1 = R1,
                                 R2 = R2,
                                 genename1 = "DDX23",
                                 genename2 = "RPS7")
pair
```
![](images/pair.png)

:neutral_face: **Notice!** We suggest use pairend interact plot fashion when the reads pair less than < 40 to avoid pixel collision.   
If you have more than 40 read pairs to show, please switch to the parallel style.

## Building blocks details:
The above functions are highly integrated and easy to use. However, the following details of each building blocks might be helpful if you want to create your own customized plot. 
### 1. Ideogram 
In gginteract you could creat an ideogram ggplot object using the following code. You need to specify the genome you're using and sub chr you're looking. This function rewrite the ideogram object in [ggbio](http://bioconductor.org/packages/release/bioc/html/ggbio.html)). 
```R
？create_ideo
ideo <- create_ideo(genome = "hg38",
                    chr = "chr4",
                    ideo.width = 400)  
```





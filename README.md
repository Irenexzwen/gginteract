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


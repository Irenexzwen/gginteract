% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/creat_ideo.R
\name{create_ideo}
\alias{create_ideo}
\title{Title Create customized ideogram ggplot object}
\usage{
create_ideo(genome = "hg38", chr, ideo.width, ideo.height, ydrift = 0,
  xdrift = 0, whratio = 16)
}
\arguments{
\item{genome}{String. Genome version eg."hg38","mm10","mm20". default("hg38")}

\item{chr}{String. chromosome name eg."chr4".}

\item{ideo.width}{Numeric. Width of the ideogram ggplot object.}

\item{ideo.height}{Numeric. Height of the ideogram, if missing, then height would be calculated by ideo.width/whratio.}

\item{ydrift}{Numeric. x axis dirft from 0.}

\item{xdrift}{Numeric. y axis dirft from 0.}

\item{whratio}{Numeric, the width/height ratio of the ideogram.}
}
\value{
ideo class
}
\description{
Title Create customized ideogram ggplot object
}
\examples{
ideo <- create_ideo(chr="chr4",ideo.width=400)

}

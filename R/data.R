#' GENE1 annotation dataframe
#'
#' A dataset containing the basic format of a gene annotation dataframe.
#' each row represents an exon with start and end coordinates of it.
#' for the use of function \{code}gene_anno it requires \strong{at least first three columns}.
#'
#' @format A data frame with 17 rows and 5 variables:
#' \describe{
#'   \item{v1}{chromosome name (could either be NC^ or Chr^)}
#'   \item{V2}{exon start}
#'   \item{V3}{enxon end}
#'   \item{V4}{exon id}
#'   \item{V5}{strand}
#' }
"GENE1_anno"

#' GENE2 annotation dataframe
#'
#' A dataset containing the basic format of a gene annotation dataframe.
#' each row represents an exon with start and end coordinates of it.
#' for the use of function \{code}gene_anno it requires \strong{at least first three columns}.
#'
#' @format A data frame with 7 rows and 5 variables:
#' \describe{
#'   \item{v1}{chromosome name (could either be NC^ or Chr^)}
#'   \item{V2}{exon start}
#'   \item{V3}{enxon end}
#'   \item{V4}{exon id}
#'   \item{V5}{strand}
#' }
"GENE2_anno"

#' READ1 bed file dataframe
#'
#' A dataset containing the read1(end1) information in bed format:
#' chr, start, end, reads_name, score, strand. Only \strong{first 4 columns} will be used later.
#' column4 reads id should be corporated with R2 file reads name, the same name indicates a pair relationship.
#'
#' @format A data frame with 25 rows and 6 variables:
#' \describe{
#'   \item{v1}{chromosome name (could either be NC^ or Chr^)}
#'   \item{V2}{reads start}
#'   \item{V3}{reads end}
#'   \item{V4}{reads id}
#'   \item{V5}{score}
#'   \item{V6}{strand}
#' }
"R1"

#' READ2 bed file dataframe
#'
#' A dataset containing the read1(end1) information in bed format:
#' chr, start, end, reads_name, score, strand. Only \strong{first 4 columns} will be used later.
#' column4 reads id should be corporated with R2 file reads name, the same name indicates a pair relationship.
#'
#' @format A data frame with 25 rows and 6 variables:
#' \describe{
#'   \item{v1}{chromosome name (could either be NC^ or Chr^)}
#'   \item{V2}{reads start}
#'   \item{V3}{reads end}
#'   \item{V4}{reads id}
#'   \item{V5}{score}
#'   \item{V6}{strand}
#' }
"R2"

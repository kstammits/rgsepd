\name{Name_to_RefSeq}
\alias{Name_to_RefSeq}

\title{
Name to RefSeq}
\description{
Lookup a HGNC symbol and return an appropriate NM##.}
\usage{
Name_to_RefSeq(x)
}
\arguments{
  \item{x}{
The HGNC symbol(s) you wish to convert.}
}
\details{
Not found gene symbols will return NA or the empty string.  }
\value{
The NM_### id numbers corresponding to the input gene names (HGNC symbols.)
}

\note{
This routine relies on bioconductor annotation package org.Hs.eg.db to ensure the most up to date mappings. }

\examples{

 Name_to_RefSeq("LSMEM2")
 #should return NM_153215
}

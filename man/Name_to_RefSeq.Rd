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
This routine relies on bioconductor annotation package org.Hs.eg.db to ensure the most up to date mappings. 
A few years after writing this function, some gene names have had their primary transcript ID changed to a new one. It's not in the old bundled dataset, and some scripts are not running properly anymore. Will have to update the dataset or use a more modern version. }

\examples{

 Name_to_RefSeq("LSMEM2")
 #should return NM_153215
}

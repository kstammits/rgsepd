\name{ExtractProjection}
\alias{ExtractProjection}
\title{
ExtractProjection
}
\description{
This function takes a completed GSEPD object with sample data, and a set of gene identifiers and produces the projection of sample expression in the sub-space. 
}
\usage{
ExtractProjection(GSEPD, txids, DRAWING=FALSE, GN=c(1,2), PRINTING=FALSE, plotTitle="")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{GSEPD}{The GSEPD parameter object. Must be post-Process.}
  \item{txids}{The transcript IDs, generally REFSEQ identifiers corresponding to rows of the counts table for this a projection is desired. In normal usage these are based on a GO Term. }
  \item{DRAWING}{ Boolean flag to draw a plot of the projection.  }
  \item{GN}{ The gene numbers: which items of the 'txids' list are to be drawn. Only the first two are used. If Drawing=FALSE, this parameter is irrelevant.}
  \item{PRINTING}{ Boolean flag to print some debug information. }
  \item{plotTitle}{ A name for this set of genes, serves as the plot's main title. }
  
}

\details{

 Primary gene set projection tool. This function calculates the vector projection and axis in a N-dimensional space of gene expression for a set of samples. When DRAWING=TRUE you will get some diagrams of the expression normalized counts. 
}

\value{
Returns a list object with four values for each sample. 
  \item{alpha}{ Distance along the axis from group1 to group2, generally 0-1, as in percent. Samples within group 1 should average zero, and samples in group 2 should average one. }
  \item{beta}{ Distance from the samples to the axis. This is a measure of goodness of fit, when the value is zero it means the sample is a linear interpolation between the comparison groups. When the value is high, the sample is not along the n-dimensional axis. }
  \item{gamma1}{ Distance from the samples to the center of group1 }
  \item{gamma2}{ Distance from the samples to the center of group2 }
  \item{Validity.Score}{ A score, 0\% through 100\%, of the segregation validity for this gene set among the two sample test groups. }
  \item{Validity.P}{ The validity score's associated p-value, empirically calculated chance of a random sample assignment creating such a strong score. }

}

\examples{

  data("IlluminaBodymap")
  data("IlluminaBodymapMeta")
  set.seed(1000) #fixed randomness
  x <- Name_to_RefSeq(c("HIF1A","EGFR","MYH7","CD33","BRCA2"))
  isoform_ids = intersect(x, rownames(IlluminaBodymap))
  rows_of_interest <- unique( c( isoform_ids ,
                                 sample(rownames(IlluminaBodymap),
                                        size=2000,replace=FALSE)))
  G <- GSEPD_INIT(Output_Folder="OUT",
                finalCounts=round(IlluminaBodymap[rows_of_interest , ]),
                sampleMeta=IlluminaBodymapMeta,
                COLORS=c("green","black","red"))
                
  G <- GSEPD_ChangeConditions( G, c("A","B")) #set testing groups first!    
  G <- GSEPD_Process( G ) #have to have processed results to plot them
   
   # looking at genes 2 and 3 will show us a view in dimensions "EGFR" and "MYH7"
   # and an axis through five dimensional space.
  ExtractProjection(GSEPD=G, txids=isoform_ids, 
    DRAWING=TRUE, PRINTING=TRUE, GN=c(2,3))
  
}

\keyword{ plot }



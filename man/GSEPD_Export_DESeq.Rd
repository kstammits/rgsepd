\name{GSEPD_Export_DESeq}
\alias{GSEPD_Export_DESeq}
\title{
Export DESeqDataSet object}
\description{
Converts from the internal matrices to a DESeq standard object. 
}
\usage{
GSEPD_Export_DESeq(G)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{G}{
 The GSEPD list object to extract a DeseqDataSet from.}
}
\details{
 Using the given GSEPD object's finalCounts and sampleMeta, a simple DESeqDataSet object is created with the default design matrix.}
\value{
a DESeqDataSet
}
\references{
DESeq2 }

\examples{
  data("IlluminaBodymap")
  data("IlluminaBodymapMeta")
  set.seed(1000) #fixed randomness
  isoform_ids <- Name_to_RefSeq(c("HIF1A","EGFR","MYH7","CD33","BRCA2"))
  rows_of_interest <- unique( c( isoform_ids ,
                                 sample(rownames(IlluminaBodymap),
                                        size=500,replace=FALSE)))
  G <- GSEPD_INIT(Output_Folder="OUT",
                finalCounts=round(IlluminaBodymap[rows_of_interest , ]),
                sampleMeta=IlluminaBodymapMeta,
                COLORS=c("green","black","red"))    
  G <- GSEPD_ChangeConditions( G, c("A","B")) #set testing groups first!           
  dds <- GSEPD_Export_DESeq(G)
  print(dds)
}

\keyword{ DESeq2 }

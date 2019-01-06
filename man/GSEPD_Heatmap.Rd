\name{GSEPD_Heatmap}
\alias{GSEPD_Heatmap}
\title{
GSEPD_Heatmap
}
\description{
Plots the heatmap to the standard display. Uses heatmap.2 from gplots to display selected genes' expression level.
}
\usage{
GSEPD_Heatmap(G,genes,cap_range=3,cellnote="log10")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{G}{The GSEPD parameter object. Must be post Process.}
  \item{genes}{rownames of finalCounts, usually isoform ID#s.}
  \item{cap_range}{z-score of most extreme color}
  \item{cellnote}{display the log10 values in each cell. No other options are supported.}
  
}

\details{
Will use GSEPD$COLORFUNCTION scaled between samples of type GSEPD$Conditions in GSEPD$sampleMeta, including others in the mix. The heatmap's dendrograms (margin trees) are computed by the heatmap.2() function's default method hclust() on the supplied data, resulting in complete linkage heirarchical clustering. Because the magnitude of gene expression varies across a wide range, and we're interested in patterns more than scale, we first normalize each gene(row) by subtracting the mean, dividing by the standard deviation, and capping the min and max to the parameter cap_range=3. The heatmap function is run with no further scaling, ensuring genes with similar differential expression profiles are clustered together. The numbers written in each cell of the heatmap are simply the normalized counts directly from DESeq2::varianceStabilizingTransformation. 
}

\value{
No return value: generates a figure.
}

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
  G <- GSEPD_Process( G ) #have to have processed results to plot them
  
  GSEPD_Heatmap(G, genes=sample(rownames(G$finalCounts),8) )
  
}

\keyword{ heatmap }
\keyword{ plot }



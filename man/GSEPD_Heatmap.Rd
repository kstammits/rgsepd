\name{GSEPD_Heatmap}
\alias{GSEPD_Heatmap}
\title{
GSEPD_Heatmap
}
\description{
Plots the heatmap to the standard display. Uses heatmap.2 from gplots to display selected genes' expression level.
}
\usage{
GSEPD_Heatmap(G,genes,cap.range=3,cellnote="log10")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{G}{The GSEPD parameter object. Must be post Process.}
  \item{genes}{rownames of finalCounts, usually isoform ID#s.}
  \item{cap.range}{z-score of most extreme color}
  \item{cellnote}{display the log10 values in each cell. no other options are supported yet.}
  
}
\details{
Will use GSEPD$COLORFUNCTION scaled between samples of type GSEPD$Conditions in GSEPD$sampleMeta, including others in the mix.
}
\value{
No return value, generates a figure.
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



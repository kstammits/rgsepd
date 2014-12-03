\name{GSEPD_PullDEG}
\alias{GSEPD_PullDEG}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Pull Differentially Expressed Genes}
\description{
After processing, if you want to easily access the differentially expressed transcript listing, this function will read in the default generated files, and apply filters as specified by the GSEPD master object (default p-values). }
\usage{
GSEPD_PullDEG(GSEPD, PTHRESH)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{GSEPD}{
The master object should have been processed already such that differentially expressed genes are readily available.}
  \item{PTHRESH}{
Specify the degree of stringency.}
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
  
  Significant_Genes <- GSEPD_PullDEG(G, PTHRESH=0.0250)
  #then do more with these identifiers:
  print(Significant_Genes)
  # GSEPD_Heatmap(G, genes= Significant_Genes )
  
}
\value{
Returns a vector of ID#, suitable to row-subsetting of the finalCounts table.
}


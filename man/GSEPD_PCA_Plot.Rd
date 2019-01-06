\name{GSEPD_PCA_Plot}
\alias{GSEPD_PCA_Plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Principle Components Analysis figure generation}
\description{After processing the pipeline, users may want to have further PCA figures generated. This function takes a completed GSEPD object and generates informative figures, based on the differentially expressed genes.}
\usage{
GSEPD_PCA_Plot(GSEPD, customColors=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{GSEPD}{
The master object, it should have already been run through GSEPD_Process().}
  \item{customColors}{
a boolean value, when FALSE, default behavior is to color points to match the test conditions. When TRUE, use sampleMeta$CustomColor column for a sample-by-sample user specification. This behavior disables the built-in legend. }
}
\value{ No return value. Generates files.}


\examples{
  data("IlluminaBodymap")
  data("IlluminaBodymapMeta")
  set.seed(1000) #fixed randomness
  isoform_ids <- Name_to_RefSeq(c("HIF1A","EGFR","MYH7","CD33","BRCA2"))
  rows_of_interest <- unique( c( isoform_ids ,
                                 sample(rownames(IlluminaBodymap),
                                        size=2000,replace=FALSE)))
  G <- GSEPD_INIT(Output_Folder="OUT",
                finalCounts=round(IlluminaBodymap[rows_of_interest , ]),
                sampleMeta=IlluminaBodymapMeta,
                COLORS=c("green","black","red"))
  G <- GSEPD_ChangeConditions( G, c("A","B")) #set testing groups first!           
  G <- GSEPD_Process( G ) #have to have processed results to plot them
  GSEPD_PCA_Plot(G)
  
}

\seealso{
GSEPD_PCA_Spec}

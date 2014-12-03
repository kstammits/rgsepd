\name{GSEPD_ChangeOutput}
\alias{GSEPD_ChangeOutput}

\title{GSEPD_ChangeOutput}
\description{
Update the stored output folder designation, and create it if necessary. This is useful if you want to change some LIMIT parameters and re-run the pipeline. Don't forget to GSEPD_Process() after changing settings.}
\usage{
GSEPD_ChangeOutput(GSEPD, newFolder)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{GSEPD}{
The initial GSEPD parameter object to update the output folder of.}
  \item{newFolder}{
The new output folder to be created.}
}
\value{
Returns the updated GSEPD parameter object.
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

  G<- GSEPD_ChangeOutput(G, "Output2")
  #G <- GSEPD_Process( G ) #would output to folder Output2
  #now tweak some settings and re-do
  G$LIMIT$LFC <- 0.25 #lower than default log-fold-change limit
  G<- GSEPD_ChangeOutput(G, "Output-Low")
  #G <- GSEPD_Process( G ) #would output to folder Output-Low

}


\author{
karl.stamm@gmail.com
}

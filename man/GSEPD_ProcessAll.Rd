\name{GSEPD_ProcessAll}
\alias{GSEPD_ProcessAll}

\title{GSEPD_ProcessAll}

\description{
 Runs each pairing within GSEPD$sampleMeta$Conditions. }
\usage{
GSEPD_ProcessAll(G)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{G}{ The GSEPD object from GSEPD_INIT() }
}
\details{ Set your GSEPD$LIMIT before running each pairwise comparison.}
\value{ Returns the last GSEPD object.}


\examples{
  data("IlluminaBodymap")
  data("IlluminaBodymapMeta")
  head(IlluminaBodymap)
  set.seed(1000) #fixed randomness
  isoform_ids <- Name_to_RefSeq(c("HIF1A","EGFR","MYH7","CD33","BRCA2"))
  rows_of_interest <- unique( c( isoform_ids ,
                                 sample(rownames(IlluminaBodymap),
                                        size=2000,replace=FALSE)))
  G <- GSEPD_INIT(Output_Folder="OUT",
                finalCounts=round(IlluminaBodymap[rows_of_interest , ]),
                sampleMeta=IlluminaBodymapMeta,
                COLORS=c("green","black","red"))    
                
#  G <- GSEPD_ProcessAll( G ) #would run across all pairs of G$Condition

}


\seealso{
GSEPD_Process

}

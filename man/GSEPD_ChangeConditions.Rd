\name{GSEPD_ChangeConditions}
\alias{GSEPD_ChangeConditions}

\title{GSEPD_ChangeConditions}
\description{
This function is an interface to set which samples are the test conditions. Don't forget to GSEPD_Process() after changing settings. If you want to systematically try each condition pairing, try GSEPD_ProcessAll() }
\usage{
GSEPD_ChangeConditions(GSEPD, newConditions)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{GSEPD}{Parameters object.}
  \item{newConditions}{ a two-item vector matching some of your sampleMeta$Conditions}
}
\details{ Interface will check if the conditions are known, then set the C2T value.}
\value{
  Returns the GSEPD parameter object with its mode set via the C2T and Conditions element of the named list. These tell later steps which sample conditions you intend on comparing.
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
        ConditionsToTest <- c("A","B")
  G <- GSEPD_ChangeConditions( G, ConditionsToTest  )
  #G <- GSEPD_Process( G ) #would test samples A vs samples B
  G <- GSEPD_ChangeConditions( G, c("A","C"))
  #G <- GSEPD_Process( G ) #would test samples A vs samples C
}

\seealso{ GSEPD_Example }

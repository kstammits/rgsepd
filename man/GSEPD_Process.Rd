\name{GSEPD_Process}
\alias{GSEPD_Process}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Processing}
\description{Primary interface, use this function to kick off the pipeline.}
\usage{
GSEPD_Process(GSEPD)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{GSEPD}{
The initialized GSEPD master object to operate on.}
}
\details{
Runs the pipeline. If any files are already present matching the generated filenames, they will be reused. If you changed a parameter that would alter the generated filenames, new ones are created. If a customization parameter is not part of the filename (like a p-value cutoff), you should change the output folder to keep new files separate.}
\value{
Returns the GSEPD object post-processed, for use in further plotting functions. Optional.
}
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
  #G <- GSEPD_Process( G ) #would run DESeq2 and GOSeq and GSEPD comparing conditions A and B
 
}

\seealso{GSEPD_INIT}

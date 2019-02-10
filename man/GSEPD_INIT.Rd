\name{GSEPD_INIT}
\alias{GSEPD_INIT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Initialization}
\description{Initializes the system, here you will pass in the count dataset and the sample metadata, before any GSEPD processing. Return value is a named list holding configurable parameters. }
\usage{
GSEPD_INIT(Output_Folder = "OUT", finalCounts = NULL, sampleMeta = NULL,
DESeqDataSet = NULL, renormalize = TRUE, vstBlind=TRUE,
  COLORS = c("green", "gray", "red"),
  C2T = "x" )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Output_Folder}{ Specify the subdirectory to hold output/generated files. Defaults to "OUT".}
  \item{finalCounts}{This must be a matrix of count data, rows are transcript IDs and columns are samples.}
  \item{sampleMeta}{The sampleMeta matrix must be passed here. It is a data frame with a row for each sample in the finalCounts matrix. Some required columns are SHORTNAME= sample nicknames; Condition= treatment group for differential expression; and Sample are the column names of finalCounts. Other columns are permitted to facilitate subsetting (not automatically supported).}
  \item{DESeqDataSet}{Data may also be included in the format of a DESeqDataSet object, this is mutually exclusive of the finalCounts/sampleMeta scheme.}
  \item{renormalize}{Boolean performance flag. Default is TRUE, which causes a normalized counts table to be computed from your given raw reads 'finalCounts'. If you set this to FALSE, then the normCounts table is preloaded with the given finalCounts input matrix, short circuiting the built-in DESeq VST, and allowing the user to specify some sort of pre-normalized dataset.}
  \item{vstBlind}{Exposes the option from DESeq2 to change the way varianceStabilizingTransformation works. According to the DESeq manual: blind=TRUE should be used for comparing samples in an manner unbiased by prior information on samples ... If many of genes have large differences in counts due to the experimental design, it is important to set blind=FALSE for downstream analysis.}
  \item{COLORS}{A three element vector of colors to make the heatmaps, the first element is the under-expressed genes, and the third element is the over-expressed genes. Defaults to green-red through gray.}
  \item{C2T}{This symbol is used in the filenames to delimit sample groups.}
  
}
\details{This function sets up the master parameter object, and therefore must be called first. This object includes all configurable parameters you can change before running the pipeline. Count data should be provided in the finalCounts matrix, with phenotype and sample data in the sampleMeta matrix. Optionally, these data may be packages in a DESeqDataSet instead. Rows with no expression are dropped at the point of loading.}
\value{Returns the GSEPD named list master object, to be used in subsequent function calls.}

\examples{
  data("IlluminaBodymap")
  data("IlluminaBodymapMeta")
  isoform_ids <- Name_to_RefSeq(c("HIF1A","EGFR","MYH7","CD33","BRCA2"))
  rows_of_interest <- unique( c( isoform_ids ,
                                 sample(rownames(IlluminaBodymap),
                                        size=2000,replace=FALSE)))
  G <- GSEPD_INIT(Output_Folder="OUT",
                finalCounts=round(IlluminaBodymap[rows_of_interest , ]),
                sampleMeta=IlluminaBodymapMeta,
                COLORS=c("green","black","red"))   
  #now ready to run:
  # G<-GSEPD_ProcessAll(G);
  
}

\seealso{GSEPD_Process}


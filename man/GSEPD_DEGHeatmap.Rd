\name{GSEPD_DEGHeatmap}
\alias{GSEPD_DEGHeatmap}

\title{
 Differentially Expressed Genes Heatmap 
}
\description{
 Generates a gene-by-subject heatmap plot of differentially expressed genes.
 }
\usage{
 GSEPD_DEGHeatmap(G)
}
\arguments{
  \item{G}{
  The GSEPD master object carries sample information and gene expression data. It should have already run Process() to be eligible.  Parameters regarding differential expression limits are set within the G$LIMIT list object. 
 }
}
\details{
 After GSEPD_Process() has created differential expression tables with known filenames, this function can read those tables and make heatmap plots for a subset of genes.  We use the N most significant genes, specified by the \code{MAX_Genes_for_Heatmap} parameter of the passed GSEPD object.
}
\value{
 This function doesn't return anything. If successful, four PDF files are created. HM and HM- are all subjects from sampleMeta and finalCounts, HMS and HMS- are only those in the test groups. The hyphen indicates a smaller unlabeled figure. In each case the data is manipulated as in \code{GSEPD_Heatmap()} such that complete linkage clustering is performed on z-score normalized genes using the normalized counts directly from DESeq2::varianceStabilizingTransformation, which are displayed in the labeled figures. 
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
  GSEPD_DEGHeatmap(G) # all parameters automatic
  
}
\keyword{ plot }

\name{GSEPD_PCA_Spec}
\alias{GSEPD_PCA_Spec}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Specialized PCA Plot}
\description{After processing the pipeline, users may want to have further PCA figures generated. This function takes a completed GSEPD object and generates informative figures. This function includes parameters to specify a particular GO-Term of interest.}
\usage{
GSEPD_PCA_Spec(GSEPD, GOT, MDATA = NULL, customColors=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{GSEPD}{
The master GSEPD object, post-processed.}
  \item{GOT}{
The GO-Term you'd like to specifically analyse. It should be found in the .MERGE file.}
  \item{MDATA}{
Optionally, pass in the .MERGE dataset, if missing, we'll try to read the already-processed file from the output directory. This option exists because reading that file repeatedly is quite slow, so you're recommended to read it in once in advance if you intend on making more than a couple GO-Term specific plots.}
  \item{customColors}{
a boolean value, when FALSE, default behavior is to color points to match the test conditions. When TRUE, use sampleMeta$CustomColor column for a sample-by-sample user specification. This behavior disables the built-in legend. }
}
\value{ No return value. Generates files.}
\note{
This function uses either princomp() or prcomp() as neccesary, depending on sample count vs gene count.}
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
  
  GOT <- "GO:0012345" # specify a GO Term you'd like to review
  
  #it should be present in the MERGE file.
  MergeFile <- list.files(G$Output_Folder, pattern="MERGE")[1]
  MDATA<-read.csv(sprintf("\%s\%s\%s", G$Output_Folder, .Platform$file.sep, MergeFile),
      as.is=TRUE,header=TRUE)
      
  GOT=MDATA$category[1] #choose a GO term that is definitely in the output data.
  
  GSEPD_PCA_Spec(G, GOT,MDATA=MDATA) 
}
\seealso{
GSEPD_PCA_Plot}


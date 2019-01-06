\name{rgsepd-package}
\alias{rgsepd-package}
\alias{rgsepd}
\docType{package}
\title{
R/Gene Set Enrichment and Projection Displays
}
\description{
R/GSEPD is a Bioinformatics package for the R programming environment that helps you disambiguate transcriptome samples (Human RNA-Seq at RefSeq IDs) by automating differential expression (DESeq), then gene set enrichment (GOSeq hg19), and finally a N-dimensional projection to quantify in which ways each sample is like either treatment group. Many exploratory tables and plots are generated for you to browse the behavior of your samples in various gene-sets (defined by GO). Sets which significantly segregate your sample conditions by bootstrapped k-means are further explored.

See the Vignette for usage examples, and minimal examples within each function's reference.
}
\details{
\tabular{ll}{
Package: \tab rgsepd\cr
Type: \tab Package\cr
Version: \tab 1.15\cr
Date: \tab 2019-01-05\cr
License: \tab GPL-3\cr
}

}
\author{
Karl D. Stamm PhD <karl.stamm@gmail.com>
}
\references{
Initially a wrapper for DESeq2 and GOSeq, late-stage processing continues with a unique geneset based sample clustering. }
\keyword{ package }
\seealso{

Example data is public human RNA-Seq from Illumina Human Bodymap2 project, aligned to GRCh37 by Ensembl and downloaded from ftp://ftp.ensembl.org/pub/release-70/bam/homo_sapiens/genebuild/ then read counts are collected by coverageBed using RefSeq.GTF. We downloaded Adipose,Blood,Heart and Skeletal Muscle, and downsampled each to one third to create artificial replicates. 

See the Vignette for usage walkthrough and results summaries.

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
  G <- GSEPD_ChangeConditions( G, c("A","B"))
 # G <- GSEPD_Process( G ) #would run DESeq2 and GOSeq and GSEPD comparing conditions A and B
 
}

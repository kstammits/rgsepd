\name{IlluminaBodymap}
\alias{IlluminaBodymap}
\docType{data}
\title{
 Sample RNA-Seq Counts data
}
\description{
 A collection of counts datasets from Illumina Human Bodymap 2.0, one sample each for adipose, blood, heart and skeletal_muscle. Four technical replicates are created by downsampling the original Illumina data. Alignment was performed by Ensembl, so the source of this dataset is ftp://ftp.ensembl.org/pub/release-70/bam/homo_sapiens/genebuild . Each of the four Human Bodymap samples are downsampled four times. Read counts are collected with Bedtools CoverageBed and a RefSeq exon annotation.    
}
\usage{data(IlluminaBodymap)}
\format{
  A data frame with 37653 observations on the following 16 variables.
  \describe{
    \item{\code{adipose.1}}{Illumina Human Bodymap 2 'Adipose' sample, downsampled to one-third. }
    \item{\code{adipose.2}}{Illumina Human Bodymap 2 'Adipose' sample, downsampled to one-third. }
    \item{\code{adipose.3}}{Illumina Human Bodymap 2 'Adipose' sample, downsampled to one-third. }
    \item{\code{adipose.4}}{Illumina Human Bodymap 2 'Adipose' sample, downsampled to one-third. }
    \item{\code{blood.1}}{Illumina Human Bodymap 2 'Blood' sample, downsampled to one-third. }
    \item{\code{blood.2}}{Illumina Human Bodymap 2 'Blood' sample, downsampled to one-third. }
    \item{\code{blood.3}}{Illumina Human Bodymap 2 'Blood' sample, downsampled to one-third. }
    \item{\code{blood.4}}{Illumina Human Bodymap 2 'Blood' sample, downsampled to one-third. }
    \item{\code{heart.1}}{Illumina Human Bodymap 2 'Heart' sample, downsampled to one-third. }
    \item{\code{heart.2}}{Illumina Human Bodymap 2 'Heart' sample, downsampled to one-third. }
    \item{\code{heart.3}}{Illumina Human Bodymap 2 'Heart' sample, downsampled to one-third. }
    \item{\code{heart.4}}{Illumina Human Bodymap 2 'Heart' sample, downsampled to one-third. }
    \item{\code{skeletal_muscle.1}}{Illumina Human Bodymap 2 'Skeletal Muscle' sample, downsampled to one-third. }
    \item{\code{skeletal_muscle.2}}{Illumina Human Bodymap 2 'Skeletal Muscle' sample, downsampled to one-third.}
    \item{\code{skeletal_muscle.3}}{Illumina Human Bodymap 2 'Skeletal Muscle' sample, downsampled to one-third.}
    \item{\code{skeletal_muscle.4}}{Illumina Human Bodymap 2 'Skeletal Muscle' sample, downsampled to one-third.}
  }
}
\value{
 A numeric matrix of read-counts from RNA-Seq, measured at transcripts by coverageBed.
}
\source{
 http://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-513/
}
\references{
 Illumina Human Bodymap 2.0. Ensembl etc.
 }
\examples{
data(IlluminaBodymap)
head(IlluminaBodymap,30)

}
\keyword{datasets}

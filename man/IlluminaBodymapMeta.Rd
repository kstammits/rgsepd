\name{IlluminaBodymapMeta}
\alias{IlluminaBodymapMeta}
\docType{data}
\title{
Metadata table for the included sample data
}
\description{
  The metadata table required to inform GSEPD of the sample/condition and abbreviated names for each column of the included 'counts' dataset. You should mirror this table's structure for your dataset. 
}
\usage{data(IlluminaBodymapMeta)}
\format{
  A data frame with 16 observations on the following 3 variables.
  \describe{
    \item{\code{Sample}}{A vector of the column names in your counts table, for the included sample data, it's four tissue types repeated four times each. For your data this must correspond to the column labels in the counts table.}
    \item{\code{Condition}}{The sample categorizations for use in differential expression, this should also be a vector the same length as the number of columns in your counts table. Here we have 'A' for each Adipose, 'B' for each muscle type, and 'C' for the blood samples.}
    \item{\code{SHORTNAME}}{Abbreviated names for each sample to appear in plots.}
  }
}
\value{
 A dataframe of sample identifiers for the rgsepd::IlluminaBodymap matrix.
}
\examples{
data(IlluminaBodymapMeta)
str(IlluminaBodymapMeta) 
}
\keyword{datasets}

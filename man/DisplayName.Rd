\name{DisplayName}
\alias{DisplayName}

\title{
DisplayName}
\description{
Convert a NM## id number to the corresponding HGNC Gene name, where available.}
\usage{
DisplayName(NM)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{NM}{
The NM_### id number, or a vector thereof.}
}

\value{
The gene's human-readable name.
}

\note{
Uses org.Hs.eg.db and pulls the first Entrez Gene ID, then that ID's associated HGNC.}

\examples{
  DisplayName("NM_005228")
}

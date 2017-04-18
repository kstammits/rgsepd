\name{DisplayName}
\alias{DisplayName}

\title{
DisplayName}
\description{
Convert a transcript id number to the corresponding gene name, where available.}
\usage{
DisplayName(txid)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{txid}{
The transcript id number, or a vector thereof.}
}

\value{
The gene's human-readable name.
}

\note{
Uses org.Hs.eg.db and pulls the first Entrez Gene ID, then that ID's associated HGNC.}

\examples{
  DisplayName("NM_005228")
}

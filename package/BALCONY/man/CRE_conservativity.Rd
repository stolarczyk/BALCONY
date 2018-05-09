\name{CRE_conservativity}
\alias{CRE_conservativity}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate cumulative relative entropy score
}
\description{
This function calculates cumulative relative entropy score according to: \href{https://www.sciencedirect.com/science/article/pii/S0022283600940361?via\%3Dihub}{Hannenhalli and Russell (2000)}.

}
\usage{
CRE_conservativity(alignment, hmmbuild_path=NULL, pairwiseAlignemnt_scores=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{alignment}{An alignment object read with \code{\link[seqinr]{read.alignment}} function}
  \item{hmmbuild_path}{(optional if running under UNIX) The aboslute path to the hmmbuild binary}
  \item{pairwiseAlignemnt_scores}{(optional) A matrix with pairwise alignment scores. For example created by \code{\link[Biostrings]{pairwiseAlignment}}. If the matrix is not provideded by the user it is calculated automatically by the function (time consuming). The sequences are extracted from the alignemnt object.}
  }

\value{
\item{score}{A vector of length equal to the length of aligned sequences}
}

\details{

\strong{PSEUDO-ALGORITHM} (According to \href{https://www.sciencedirect.com/science/article/pii/S0022283600940361?via\%3Dihub}{Hannenhalli and Russell (2000)}):

\enumerate{

    \item (If score matrix is not provided) Run pairwise alignments for all available sequences in the input MSA and save scores to a matrix

    \item (If score matrix is not provided) Calculate a distance matrix based off of the alignment scores one

    \item Perform hierarchical clustering on the distance matrix (UPGMA method)

    \item Get the sequence clusters

    \item Divide the alignment into \code{sub_groups} which are the clusters

    \item Run hmmbuild for \code{whole_alignment} without \code{sub-group} and \code{sub_group}

    \item Calculate relative entropy using these two as indicated in the Reference and repeat for each \code{sub_group}

    \item Calculate the cumulative relative entropy
}

\strong{hmmbuild program}:


This function uses hmmbuild program of \href{http://www.hmmer.org/}{HMMER} suite for HMM profile generation for MSA.

We recommend downloading and installing HMMER by following the instructions and steps in the \href{http://hmmer.org/download.html}{ HMMER installation website }.
}
\references{
Hannenhalli, S. S. & Russell, R. B. Analysis and prediction of functional sub-types from protein sequence alignments11Edited by J. Thornton. Journal of Molecular Biology 303, 61–76 (2000).

}
\author{
Michal Stolarczyk & Alicja Pluciennik
}

\seealso{
\code{\link{consensus}}, \code{\link{cons2seqs_ident}}, \code{\link[seqinr]{read.alignment}}

}
\examples{
#No example due to external software requirements
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ cumulative relative entropy }% use one of  RShowDoc("KEYWORDS")
\keyword{ conservation }% __ONLY ONE__ keyword per line

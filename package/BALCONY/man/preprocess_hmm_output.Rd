\name{preprocess_hmm_output}
\alias{preprocess_hmm_output}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
preprocess HMM output
}
\description{
Preprocessing od HMMER output file to calculate CRE.
}
\usage{
preprocess_hmm_output(hmm_out)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hmm_out}{ path to ouptut file}
}

\value{
Returns list of
 \item{probabilities }{probabilities extracted from file}
 \item{alignment_positions }{index of each alignment position}
%% ...
}
\author{
Michal Stolarczyk & Alicja Pluciennik
}

\seealso{
CRE_conservativity()
}
\examples{
#No example due to external software requirements

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hmmer }% __ONLY ONE__ keyword per line

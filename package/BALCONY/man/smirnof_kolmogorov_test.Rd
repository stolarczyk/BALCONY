\name{smirnof_kolmogorov_test}
\alias{smirnof_kolmogorov_test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Perform Smirnof-Kolmogorov test for structural data
}
\description{
This function facilitates the comparion of conservativity of structure of interest with the rest of the protein. For example comparison of tunnel conservativity with overall protein conservativity.
}
\usage{
smirnof_kolmogorov_test(protein_cons, structure_cons, alternative, pdb_name, range = NULL, make_plot = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{protein_cons}{A list of calculated entropy scores (vectors of numeric values). The output of get_get_protein_entropy() function}
  \item{structure_cons}{
%%     ~~Describe \code{structure_cons} here~~
}
  \item{alternative}{
%%     ~~Describe \code{alternative} here~~
}
  \item{pdb_name}{
%%     ~~Describe \code{pdb_name} here~~
}
  \item{range}{
%%     ~~Describe \code{range} here~~
}
  \item{make_plot}{
%%     ~~Describe \code{make_plot} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (protein_cons, structure_cons, alternative, pdb_name,
    range = NULL, make_plot = NULL)
{
    if (is.null(make_plot)) {
        make_plot <- T
    }
    metrics_count = length(protein_cons)
    structures_count = length(structure_cons)
    structure_names = c()
    for (i in seq(1, structures_count)) {
        structure_names[i] = paste("stru", i)
    }
    alt_hip = c("two.sided", "less", "greater")[alternative]
    pval_mtx = matrix(NA, nrow = metrics_count, ncol = structures_count)
    rownames(pval_mtx) <- names(protein_cons)
    colnames(pval_mtx) <- structure_names
    for (i in seq(1, metrics_count)) {
        for (j in seq(1, structures_count)) {
            reference = protein_cons[[i]][-c(range, structure_cons[[j]][[2]])]
            temp = structure_cons[[j]][[1]][i, ]
            pval_mtx[i, j] = ks.test(reference, temp, alternative = alt_hip)$p.value
        }
    }
    if (make_plot == T) {
        colors = rainbow(structures_count)
        for (i in seq(1, metrics_count)) {
            par(new = F)
            for (j in seq(1, structures_count)) {
                reference = protein_cons[[i]][-c(range, structure_cons[[j]][[2]])]
                cumulative_distribution = ecdf(reference)
                temp = ecdf(structure_cons[[j]][[1]][i, ])
                plot(cumulative_distribution, xlim = c(0, 1),
                  ylim = c(0, 1), xlab = "entropy", ylab = "cumulative probability",
                  col = alpha("slategray", 0.7), main = paste("CDF of",
                    names(protein_cons)[i], " for", pdb_name,
                    "structures"))
                par(new = T)
                plot(temp, xlim = c(0, 1), ylim = c(0, 1), col = alpha(colors[j],
                  0.7), ylab = "", xlab = "", main = "")
                legend("bottomright", c(pdb_name, structure_names),
                  lty = c(1, 1, 1, 1), lwd = c(2.5, 2.5), col = alpha(c("slategray",
                    colors), 0.8))
            }
        }
    }
    return(pval_mtx)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

\name{simGeno}
\alias{simGeno}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
simulate genotypes with admixture
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
simGeno(z1, z2, a1, a2, f)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{z1}{
%%     ~~Describe \code{z1} here~~
}
  \item{z2}{
%%     ~~Describe \code{z2} here~~
}
  \item{a1}{
%%     ~~Describe \code{a1} here~~
}
  \item{a2}{
%%     ~~Describe \code{a2} here~~
}
  \item{f}{
%%     ~~Describe \code{f} here~~
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
function (z1, z2, a1, a2, f) 
{
    npop <- ncol(f)
    N <- nrow(f)
    f1 <- f \%*\% (norm(a1 - z2 - z1/2))
    f2 <- f \%*\% (norm(a2 - z2 - z1/2))
    hap11 <- rbinom(N, 1, f1)
    hap12 <- rbinom(N, 1, f1)
    hap21 <- rbinom(N, 1, f2)
    hap22 <- rbinom(N, 1, f2)
    states <- sample(0:(2 * npop), N, replace = T, prob = c(1 - 
        sum(z1) - sum(z2), z1, z2))
    for (p in 1:npop) {
        s1 <- rbinom(N, 1, f[, p])
        hap11[states == p] <- s1[states == p]
        hap21[states == p] <- s1[states == p]
        s21 <- rbinom(N, 1, f[, p])
        s22 <- rbinom(N, 1, f[, p])
        hap11[states == p + npop] <- s21[states == p + npop]
        hap21[states == p + npop] <- s21[states == p + npop]
        hap12[states == p + npop] <- s22[states == p + npop]
        hap22[states == p + npop] <- s22[states == p + npop]
    }
    return(cbind(hap11 + hap12, hap21 + hap22))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

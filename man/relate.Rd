\name{relate}
\alias{relate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
function for estimating relatedness for admixted individuals using genomic markers
}
\description{
Package for estimating relatedness for admixted individuals using genomic markers.
}
\usage{
relate(geno1, geno2, a1, a2, freq, maxIter = 300, useSq = 1, tol = 0.002, tolStop = 1e-06,start=c(0.7,0.2,0.1))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{geno1}{
 genotypes for individual 1. An array of genotypes coded as 0,1,2 (e.g. number of minor alleles)
}
  \item{geno2}{
 genotypes for individual 2. An array of genotypes coded as 0,1,2 (e.g. number of minor alleles)
}
  \item{a1}{
admixture proportions for individual 1. An array of length K (number of
populations) that sums to 1. 
}
  \item{a2}{
admixture proportions for individual 2. An array of length K (number of
populations) that sums to 1. 
}
\item{freq}{
  population specific allele frequencies of the minor allele. A matrix
  with K columns and N   rows. Each entry is the allele frequency for a sites in a specific population
}
  \item{maxIter}{
maximum number of iterations
}
  \item{useSq}{
use accelerated EM
}
  \item{tol}{
The mimimum admixture proportions that is counted as non zero (used for speed)
}
  \item{tolStop}{
    Stopping criteria for the algorithm
  }
  \item{start}{
    starting point for the optimization
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  returns a list with the estimated relatedness coefficents and the
  number of iterations used. 

}
\references{
Ida Moltke and Anders Albrechtsen
}
\author{
Ida Moltke and Anders Albrechtsen
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
#number of SNPs
Nsnp<-100000
#frequencies for 3 populations
f1<-runif(Nsnp)
f2<-runif(Nsnp)
f3<-runif(Nsnp)
f<-cbind(f1,f2,f3)
#relatedness coefficients in 3 populations
z2=c(0.1,0.1,0)
z1=c(0.3,0.15,0.1)

cat("Simulated relatedness c(k0,k1,k2)\n")
print(trueK<-c(1-sum(z1)-sum(z2),sum(z1),sum(z2)))

## admixture proportions
a1<-z2+z1/2+c(0.3,0.025,0.2)
a2<-z2+z1/2+c(0.3,0.125,0.1)
cat("Simulated admixture proportions\n")
print(a1)
print(a2)

geno<-simGeno(z1,z2,a1,a2,f)

r<-relate(geno[,1],geno[,2],a1,a2,f)
cat("Estimated relatedness\n")
print(r)

}
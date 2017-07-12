
"relate" <-
function (geno1,geno2,a1,a2,freq,maxIter=300,useSq=1,tol=0.002,tolStop=0.000001,start=c(0.7,0.2,0.1))
{


N<-length(geno1)
K=length(a1)
if(nrow(freq)!=N)
  stop("number of columns in freq does not match the number of sites")

if(ncol(freq)!=K)
  stop("number of rows in freq does not match the number of populations")

if(any(freq<0|freq>1))
  stop("freq must be between 0 and 1")
if(!all(geno1%in%c(0:2)))
  stop("geno1 but only contain 0,1,2")
if(!all(geno2%in%c(0:2)))
  stop("geno2 but only contain 0,1,2")







res<-.Call("relateC",
          freq,as.double(start),as.integer(geno1),as.integer(geno2),
	as.double(a1),as.double(a2),as.integer(useSq),as.double(tolStop),
	as.integer(maxIter),as.integer(K),as.double(tol),PACKAGE="relateAdmix")


return(res)
}


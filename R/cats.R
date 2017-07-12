"cats" <-
function (freq=0.5,freq2=-1,ncases=500,ncontrols=500,ncases2=500,ncontrols2=500,risk=1.5,risk2=-1,pisamples=-1,prevalence=0.1,prevalence2=-1,additive=0,recessive=0,dominant=0,multiplicative=1,alpha=0.0000001,pimarkers=0.00316)
{


model<-c(additive,recessive,dominant,multiplicative)

if(sum(model==1)!=1)
  stop("chose only one model. e.i. one model must be 1 the others 0")
if(sum(model==0)!=3)
  stop("chose only one model. e.i. one model must be 1 the others 0")


if(freq<0|freq>1)
  stop("freq must be between 0 and 1")
if((freq2<0|freq2>1)&freq2!=-1)
  stop("freq2 must be between 0 and 1 (or undefined as -1)")
if((pisamples<0|pisamples>1)&pisamples!=-1)
  stop("pisamples must be between 0 and 1")
if((prevalence2<0|prevalence2>1)&prevalence2!=-1)
  stop("prevalence2 must be between 0 and 1 (or undefined as -1)")
if(alpha<0|alpha>1)
  stop("alpha must be between 0 and 1")
if(prevalence<0|prevalence>1)
  stop("prevalence must be between 0 and 1")
if(pimarkers<0|pimarkers>1)
  stop("pimarkers must be between 0 and 1")
if(ncases!=as.integer(ncases)|ncases<0)
  stop("ncases must be a positive integer")
if(ncases2!=as.integer(ncases2)|ncases2<0)
  stop("ncases2 must be a positive integer")
if(ncontrols!=as.integer(ncontrols)|ncontrols<0)
  stop("ncontrols must be a positive integer")
if(ncontrols2!=as.integer(ncontrols2)|ncontrols2<0)
  stop("ncontrols2 must be a positive integer")
if(risk<0)
  stop("risk must be positive")
if(risk2<0&risk2!=-1)
  stop("risk2 must be positive(or undefined as -1)")



res<-.Call("cats",
          as.double(freq),as.double(freq2),as.integer(ncases),as.integer(ncontrols),
	as.integer(ncases2),as.integer(ncontrols2),as.double(risk),as.double(risk2),
	as.double(pisamples),as.double(prevalence),as.double(prevalence2),
	as.integer(additive),as.integer(recessive),as.integer(dominant),
	as.integer(multiplicative),as.double(alpha),as.double(pimarkers),PACKAGE="relateAdmix")



options<-cbind(freq,freq2,ncases=ncases,ncontrols=ncontrols,ncases2=ncases2,ncontrols2=ncontrols2,risk,risk2,pisamples,prevalence,prevalence2,additive,recessive,dominant,multiplicative,alpha,pimarkers)

result<-list(P.one.study=res[1,1],P.first.stage=res[2,1],P.rep.study=res[3,1],P.joint.min=res[4,1],P.joint=res[5,1],pi=res[6,1],T.one.study=res[7,1],T.first.stage=res[8,1],T.second.stage.rep=res[9,1],T.second.stage.joint=res[10,1],E.Disease.freq.cases1=res[11,1],E.Disease.freq.controls1=res[12,1],E.Disease.freq.cases2=res[13,1],E.Disease.freq.controls2=res[14,1],options=options)
class(result)<-"CATS"
return(result)
}


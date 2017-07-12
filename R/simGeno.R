simGeno <-
function(z1,z2,a1,a2,f){
  npop<-ncol(f)
  N<-nrow(f)
  f1<-f%*%((a1-z2-z1/2)/sum(a1-z2-z1/2)) #none IBD frequencies
  f2<-f%*%((a2-z2-z1/2)/sum(a2-z2-z1/2)) #none IBD frequencies

  #two unrelated individuals
  hap11<-rbinom(N,1,f1)
  hap12<-rbinom(N,1,f1)
  hap21<-rbinom(N,1,f2)
  hap22<-rbinom(N,1,f2)

  #samples state (0 is unrelated)
  states<-sample(0:(2*npop),N,replace=T,prob=c(1-sum(z1)-sum(z2),z1,z2))

  for(p in 1:npop){
    #ibd=1
    s1<-rbinom(N,1,f[,p])
    hap11[states==p] <-s1[states==p]
    hap21[states==p] <-s1[states==p]
    ##ibd=2
    s21<-rbinom(N,1,f[,p])
    s22<-rbinom(N,1,f[,p])
    hap11[states==p+npop] <-s21[states==p+npop]
    hap21[states==p+npop] <-s21[states==p+npop]
    hap12[states==p+npop] <-s22[states==p+npop]
    hap22[states==p+npop] <-s22[states==p+npop]
#    print(mean(states==p))
  }

  return(cbind(hap11+hap12,hap21+hap22))
}

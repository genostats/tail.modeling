###Calcul des estimations ----

p_estimate<- function(Nsim=150,queue=500,Nperm=1e6,Z,method){
  
  source("C:/Users/Marion/Documents/Stage/packageBC/BC_functions.R")
  
  P<-rep(0,Nsim)
  Zsim<-list() #liste qui contient les Nsim échantillons permutés Nperm fois de l'echantillon ini Z
  for (i in 1:Nsim){
    set.seed(i)
    perm<-sample(seq(length(Z)),length(Z),replace=FALSE)
    Zsim<-c(Zsim,list(Z[perm]))
  }
  
  if (method == "BC") {
    lbda<-rep(0,Nsim)
    for(i in 1:Nsim){
      pval<-calcul_p(Zsim[[i]],queue)
      P[i]<-pval$Pbc_z
      lbda[i]<-pval$lbda
    }
  }
  if (method == "GPD") {
    k<-rep(0,Nsim)
    a<-rep(0,Nsim)
    for(i in 1:Nsim){
      pval<-calcul_p(Zsim[[i]],queue)
      P[i]<-pval$Pgpd
      a[i]<-pval$a
      k[i]<-pval$k
      
    }
  }
  
  b<-boxplot( -log10(P),main=paste0("Estimate P-value with ",method))
  print("Quantiles")
  print(b$stats)  
  
  plot(-log10(P),
       col="red",
       type="l",
       ylab="",
       main=paste0("Estimate probability with ",method, " (-log10(P))")
  )
}
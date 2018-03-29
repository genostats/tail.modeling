##script de calcul des p valeurs pour différentes Pperm et loi de distribution pour Zobs

### Fonctions ----
source(file="p_value_estimate_funct.R") #charger les fonctions

###Calcul ----

Pg<-rep(0,length(Nsim))
Pe<-rep(0,length(Nsim))
Pli<-rep(0,length(Nsim))
Pb<-rep(0,length(Nsim))
k<-rep(0,length(Nsim))

for(i in 1:length(Nsim)){
  pval<-calcul_p(Zsim[[i]],queue)
  Pg[i]<-pval$Pgpd
  Pe[i]<-pval$Pecdf
  Pli[i]<-pval$Plin
  Pb[i]<-pval$Pbc
  k[i]<-pval$k
}

##Graphiques ----
# boxplot( -log10(Pe), -log10(Pg),-log10(Pli) ,-log10(Pb),
#          names=c("ECDF","GPD","Reg Lin","BoxCox"))
boxplot( -log10(Pli) ,-log10(Pb),
         names=c("Reg Lin","BoxCox"))
abline(h=-log10(Pperm),col="red")
abline(h=-0.9*log10(Pperm),col="red",lty=3)
abline(h=-1.1*log10(Pperm),col="red",lty=3)

plot(Nsim,-log10(Pli),
     col="red",
     type="l",
     ylab="",
     main=paste0("Loi ",loi),
     ylim=c(-log10(Pperm)-1,-log10(Pperm)+1)
     )
lines(Nsim,-log10(Pb),col="black")
# lines(Nsim,-log10(Pe),col="green")
# lines(Nsim,-log10(Pg),col="coral3")
#lines(Nsim,rep(-log10(Pperm),length(Nsim)),col="blue")
abline(h=-0.9*log10(Pperm),col="blue",lty=3)
abline(h=-1.1*log10(Pperm),col="blue",lty=3)
legend("topright",
       # legend=c("reg lin","BoxCox","ECDF","GPD",paste0("Pperm = ",Pperm)),
       # fill=c("red","black","green","coral3","blue"),horiz = T)
       legend=c("reg lin","BoxCox","Pperm"),
       fill=c("red","black","blue"),horiz = T)
grid(lty = 1,ny=NA,nx=15)







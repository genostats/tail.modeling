##script de calcul des p valeurs pour différentes Pperm et loi de distribution pour Zobs

### Fonctions ----
source(file="p_value_estimate_funct.R") #charger les fonctions

###Calcul ----

Pg<-rep(0,length(Nsim))
Pe<-rep(0,length(Nsim))
Pli<-rep(0,length(Nsim))
Pb<-rep(0,length(Nsim))
Pbz<-rep(0,length(Nsim))
k<-rep(0,length(Nsim))
lbda<-rep(0,length((Nsim)))

for(i in 1:length(Nsim)){
  pval<-calcul_p(Zsim[[i]],queue)
  Pg[i]<-pval$Pgpd
  Pe[i]<-pval$Pecdf
  Pli[i]<-pval$Plin
  Pb[i]<-pval$Pbc
  Pbz[i]<-pval$Pbc_z
  k[i]<-pval$k
  lbda[i]<-pval$lbda
}

##Graphiques ----

jpeg(paste0("C:/Users/Marion/Documents/Stage/Boxplots_EMV/boxplot_loi_",loi,"_",Pperm, ".jpeg"),res = 450, height = 12, width = 16, units = 'cm')
b<-boxplot( -log10(Pli) ,-log10(Pg),-log10(Pbz),
         names=c("Reg Lin","GPD","Box-Cox Z"),main=paste0("Boxplots pour la loi ",loi," et la pvaleur ",Pperm))
abline(h=-log10(Pperm),col="red")
abline(h=-0.9*log10(Pperm),col="red",lty=3)
abline(h=-1.1*log10(Pperm),col="red",lty=3)
dev.off()

sink(paste0("C:/Users/Marion/Documents/Stage/Boxplots_EMV/tableau_boxplot_loi_",loi,"_",Pperm, ".txt"))
print("Quantiles et valeurs extrêmes pour les méthodes ML (1), GPD (2) et BC (3)")
print(b$stats)
sink()



plot(-log10(Pli),
     col="red",
     type="l",
     ylab="",
     main=paste0("Probabilite estimee avec 10^6 permutations")
     )
lines(-log10(Pg),col="black")
lines(-log10(Pbz),col="forestgreen")

# abline(h=-0.9*log10(Pperm),col="blue",lty=3)
# abline(h=-1.1*log10(Pperm),col="blue",lty=3)
legend("topright",
       legend=c("reg lin","GPD","Box-Cox Z"),
       fill=c("red","black","forestgreen"),horiz = T)
#grid(lty = 1,ny=NA,nx=15)

#plot(Nsim,lbda,col="blueviolet",main="Box-Cox parameter lambda")

print("Queue lourde ?")
print(critereLourd(Zsim[[nbsim]],end))
print(paste0("Mean lambda: ",mean(lbda)))




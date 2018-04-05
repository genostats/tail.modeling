#### Implementation de la methode GDP 
#### en utilisant la methode PWM pour estimer les parametres de Pareto a et k
#### pour differentes valeur reelles de la pvaleur Pperm = la vraie Pvaleur 
#### pour differentes distributions



### Variables ----
#p valeur visee
Pperm <- 1e-9

# loi simulee

##critere de queue lourde: differences entre les quartiles

d<-diff(qchisq( c(0.97,0.98, 0.99,0.999),df=3))
R<-(d[3]-d[2])/(d[2]-d[1])
if (R>5.6)
  print("Queue de distribution lourde, prendre le log des statistiques de test")

Zobs<-log(qcauchy(Pperm,lower.tail = F))
#pf(exp(Zobs),df1=5,df2=10,lower.tail=F)

# Parametres de simulations par permutation
start<-1.4e3
end<-1e6
nbsim<-300
queue<-500

Nsim<-seq(start,end,length.out=nbsim)
Zsim<-list() #liste qui contient les nbsim simulations contenant Nsim[i] permutations
for (i in 1:length(Nsim)){
  Zsim<-c(Zsim,list(log(rcauchy(Nsim[i]))))
}

### Fichier de calcul ----
source(file="p_value_estimate_script.R") #charger le script de calcul

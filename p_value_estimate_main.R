#### Implementation de la methode GDP 
#### en utilisant la methode PWM pour estimer les parametres de Pareto a et k
#### pour differentes valeur reelles de la pvaleur Pperm = la vraie Pvaleur 
#### pour differentes distributions



### Variables ----
#p valeur visee
Pperm <- 1e-5

# loi simulee

Zobs<-qf(Pperm,df1=5,df2=10,lower.tail = F)
#pf(exp(Zobs),df1=5,df2=10,lower.tail=F)

# Parametres de simulations par permutation
start<-1.4e3
end<-5e5
nbsim<-300
queue<-500

Nsim<-seq(start,end,length.out=nbsim)
Zsim<-list() #liste qui contient les nbsim simulations contenant Nsim[i] permutations
for (i in 1:length(Nsim)){
  Zsim<-c(Zsim,list(rf(Nsim[i],df1=5,df2=10)))
}

### Fichier de calcul ----
source(file="p_value_estimate_script.R") #charger le script de calcul

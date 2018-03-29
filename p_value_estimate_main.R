#### Implementation de la methode GDP 
#### en utilisant la methode PWM pour estimer les parametres de Pareto a et k
#### pour differentes valeur reelles de la pvaleur Pperm = la vraie Pvaleur 
#### pour differentes distributions



### Variables ----
#p valeur visee
Pperm <- 1e-9

# loi simulee
loi<-""
Zobs<-qlnorm(Pperm,sdlog=2,lower.tail = FALSE) 
#pnorm(Zobs,lower.tail=F)

# Parametres de simulations par permutation
start<-1.4e3
end<-1e6
nbsim<-300
queue<-500

Nsim<-seq(start,end,length.out=nbsim)
Zsim<-list() #liste qui contient les nbsim simulations contenant Nsim[i] permutations
for (i in 1:length(Nsim)){
  Zsim<-c(Zsim,list(rlnorm(Nsim[i],sdlog=2)))
}

### Fichier de calcul ----
source(file="p_value_estimate_script.R") #charger le script de calcul

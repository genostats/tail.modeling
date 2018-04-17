#### Implementation de la methode GDP 
#### en utilisant la methode PWM pour estimer les parametres de Pareto a et k
#### pour differentes valeur reelles de la pvaleur Pperm = la vraie Pvaleur 
#### pour differentes distributions



### Variables ----
#p valeur visee
loi <- "F"
P <- c(1e-4,1e-5,1e-6,1e-7,1e-8)

for (Pperm in P) {
# loi simulee
print(Pperm)

Zobs<-qf(Pperm,df1=5,df2=10,lower.tail = F)
#pf(exp(Zobs),df1=5,df2=10,lower.tail=F)

# Parametres de simulations par permutation
# start<-1.4e3
end<-1e6  #useful for critere lourd

nbsim<-150
queue<-500

#Nsim<-seq(start,end,length.out=nbsim)
Nsim<-rep(1e6,nbsim)
Zsim<-list() #liste qui contient les nbsim simulations contenant Nsim[i] permutations
for (i in 1:length(Nsim)){
  Zsim<-c(Zsim,list(rf(Nsim[i],df1=5,df2=10)))
}

# critereLourd<-function(z,N){
#   z<-sort(z)
#   p<-seq(1,N)/N #proba empiriques p(X<Z)
#   d<-diff(c(z[sum(p<=0.97)],z[sum(p<=0.98)], z[sum(p<=0.99)],z[sum(p<=0.999)])) #diff des ecarts des quantiles empiriques
#   R<-(d[3]-d[2])/(d[2]-d[1])
#   return(R)}
# 
# R<-rep(0,length(Nsim))
# for (i in 1:length(Nsim))
#   R[i]<-critereLourd(Zsim[[i]],end)
# boxplot(R)  
  
  
### Fichier de calcul ----
source(file="p_value_estimate_script.R") #charger le script de calcul
}
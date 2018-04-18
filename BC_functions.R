####Fonctions utiles pour estimer la p-valeur 

#calcul du nombre de stats simulees superieures a Zobs
nb_exc<-function(x0,z){  ##pas besoin d'une boucle en fait
  M<-length(z[z>=x0])
  return(M)
}

#GPD ----
# estimateur PWM pour les parametres de la GPD
PWM <- function(z) {
  n <- length(z)
  p <- (seq.int(n)-0.35)/n
  t <- mean( (1-p)*z )
  mu <- mean(z)
  k <- mu/(mu-2*t)-2
  a <- 2*mu*t/(mu-2*t)
  list(a = a, k = k)
}

# estimateur du max de vraisemblance (via une vraisemblance profilée)
k_func <- function(u, z) -mean(log(1-u*z))

Lp <- function(u,z) {
  n <- length(z)
  if(u == 0) 
    return( -n*(log(mean(z)) + 1) )
  k <- k_func(u,z)
  n*(log(u/k) + k - 1)
}

EMV <- function(z) {
  u <- optimize(function(u) Lp(u,z), c(-10, min(1/z)), maximum=TRUE )$maximum
  if(u == 0) 
    return(list( a = mean(z), k = 0))
  k <- k_func(u, z)
  list(a = k/u, k = k)
}


#Fonction de repartition de Pareto
FGPD<-function(z,zexc,estim){
  if (estim == "EMV")
    coeff<-EMV(zexc)
  else
    coeff<-PWM(zexc)
  a<-coeff$a
  k<-coeff$k
  if (k!= 0)
    return(1-(1-k*z/a)^(1/k))
  else
    return(1-exp(-z/a))
}

#calcul de Pgpd
PGPD<-function(x0,zexc,y,seuil,estim){
  M<-nb_exc(x0,y)
  if(M >= 10)
    return(M/length(y))
  else
    return(length(zexc)/length(y)*(1-FGPD(x0-seuil,zexc,estim)))
}
##Box-Cox transformation ----

BoxCox<-function(z,lambda){
  if (lambda != 0)
    return((z^lambda - 1)/lambda)
  else
    return(log(z))
}

#BC applique a la statistique z

BoxCox_lm <- function(Y,X) {  #estimation de lambda par les moindres carrés
  f <- function(lambda) {
    if(lambda == 0)
      X <- log(X)
    else
      X <- (X**lambda - 1)/lambda
    
    b <-  cov(X,Y) / var(X)
    a <- mean(Y) - b*mean(X)
    sum( (Y - a - b*X)**2 )
  }
  optimize(f,c(0,1))$minimum
}

PBC_Z<-function(z,Ntot,N,param=NA,Zobs){  
  p<-seq(N,1)/Ntot
  if (is.na(param))
    lambda <- BoxCox_lm(log(-log(p)),log(z))
  else
    lambda <- param
  coeffs <- lm(  log(-log(p)) ~ BoxCox(log(z),lambda) )$coefficients
  return(list(p = exp(-exp(sum( coeffs*c(1, BoxCox(log(Zobs),lambda))))), 
              interc = coeffs[[1]], 
              pente = coeffs[[2]],
              lbda = lambda))
  
}

## regroupement des p valeurs des differentes methodes ----

calcul_p<-function(zsim,Ntail,estim,Zobs,param){  
  # les Ntail plus grandes valeurs (les dernieres)
  z1 <- tail( sort(zsim) , Ntail )
  #seuil pour la GDP
  t<-(z1[1] + head(sort(zsim,decreasing =TRUE),Ntail+1)[Ntail+1])/2
  #calcul des excedents de la GDP, ceux qui lui sont sup?rieurs
  zgpd<-z1-t
  zgpd<-zgpd[zgpd>0] #uniquement ceux superieurs au seuil
  
  if (estim == "PWM")
    return(list(Pgpd=PGPD(Zobs,zgpd,zsim,t,estim),
                Pbc_z = PBC_Z(z1,length(zsim),Ntail,param,Zobs)$p,
                lbda = PBC_Z(z1,length(zsim),Ntail,param,Zobs)$lbda,
                a=PWM(zgpd)$a,
                k=PWM(zgpd)$k))
  else
    return(list(Pgpd=PGPD(Zobs,zgpd,zsim,t,estim),
                Pbc_z = PBC_Z(z1,length(zsim),Ntail,param,Zobs)$p,
                lbda = PBC_Z(z1,length(zsim),Ntail,param,Zobs)$lbda,
                a=EMV(zgpd)$a,
                k=EMV(zgpd)$k))
    
}

###Calcul des estimations ----

p_estimate<- function(Nsim=150,queue=500,Nperm=1e6,Z,method,Zobs,estim,param){
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
      pval<-calcul_p(Zsim[[i]],queue,estim,Zobs,param)
      P[i]<-pval$Pbc_z
      lbda[i]<-pval$lbda
    }
  }
  if (method == "GPD") {
    k<-rep(0,Nsim)
    a<-rep(0,Nsim)
    for(i in 1:Nsim){
      pval<-calcul_p(Zsim[[i]],queue,estim,Zobs,param)
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


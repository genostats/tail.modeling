####Fonctions utile pour estimer la p-valeur 

### Librairie a charger ----
library(MASS)

### Fonctions ----

#calcul du nombre de stats simulees superieures a Zobs
nb_exc<-function(x0,z){  ##pas besoin d'une boucle en fait
  M<-length(z[z>=x0])
  return(M)
}
#calcul de Pecdf----
PECDF<-function(x0,y){
  return((1+nb_exc(x0,y))/length(y))
  
}

#GPD ----
# estimateur PWM pour les param?tres de la GPD
PWM <- function(z) {
  n <- length(z)
  p <- (seq.int(n)-0.35)/n
  t <- mean( (1-p)*z )
  mu <- mean(z)
  k <- mu/(mu-2*t)-2
  a <- 2*mu*t/(mu-2*t)
  list(a = a, k = k)
}

#Fonction de repartition de Pareto
FGPD<-function(z,zexc){
  coeff<-PWM(zexc)
  a<-coeff$a
  k<-coeff$k
  if (k!= 0)
    return(1-(1-k*z/a)^(1/k))
  else
    return(1-exp(-z/a))
}

#calcul de Pgpd
PGPD<-function(x0,zexc,y,seuil){
  M<-nb_exc(x0,y)
  if(M >= 10)
    return(M/length(y))
  else
    return(length(zexc)/length(y)*(1-FGPD(x0-seuil,zexc)))
}

#le modele lineaire ----
PML<-function(z,Ntot,N){
  p <- seq(N,1)/Ntot
  # le mod?le ~ lin?aire 
  # estimation de P( Z > Zobs) par r?gression lin?aire
  coeffs <- lm(  log(-log(p)) ~ log(z) )$coefficients
  return(exp(-exp( sum( coeffs*c(1, log(Zobs)) ) )))
}

#Box-Cox transformation ----

BoxCox<-function(z,lambda){
  if (lambda != 0)
    return((z^lambda - 1)/lambda)
  else
    return(log(z))
}

PBC<-function(z,Ntot,N){
  p<-seq(N,1)/Ntot
  b<-boxcox( log(-log(p)) ~ log(z),plotit=F) ##which lambda to choose
  lambda <- b$x[which.max(b$y)]
  coeffs <- lm(  BoxCox(log(-log(p)),lambda) ~ log(z) )$coefficients
  Bp<-sum( coeffs*c(1, log(Zobs)))
  if (lambda != 0)
    return(exp(-exp((Bp*lambda+1)^(1/lambda))))
  else
    return(exp(-exp(exp(Bp))))
}


## regroupement des p valeurs des differentes methodes ----

calcul_p<-function(zsim,Ntail){  
  # les Ntail plus grandes valeurs (les derni?res)
  z1 <- tail( sort(zsim) , Ntail )
  #seuil pour la GDP
  t<-(z1[1] + head(sort(zsim,decreasing =TRUE),Ntail+1)[Ntail+1])/2
  #calcul des excedents de la GDP, ceux qui lui sont sup?rieurs
  zgpd<-z1-t
  zgpd<-zgpd[zgpd>0] #uniquement ceux sup?rieurs au seuil
  
  return(list(Pecdf = PECDF(Zobs,zsim),
              Pgpd=PGPD(Zobs,zgpd,zsim,t),
              Plin=PML(z1,length(zsim),Ntail),
              Pbc = PBC(z1,length(zsim),Ntail),
              a=PWM(zgpd)$a,
              k=PWM(zgpd)$k))
}


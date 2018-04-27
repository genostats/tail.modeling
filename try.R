library(dplyr)
library(PvalueEsti)

### Normal distribution ----
data<-read.csv2(file = "C:/Users/Marion/Documents/s057.csv") #tailles de chaussures
data2<-data %>%
        group_by(taille) %>%
        count()
        
plot(data2$taille, data2$n)
mean(data$taille)
sd(data$taille)

hist(rnorm(nrow(data),mean=46,sd=4))
pnorm(55,mean=46,sd=4,lower.tail=F)

# p_estimate(Z=data$taille,
#            Zobs=55,
#            Nsim=1,
#            method="GPD",
#            estim="PWM")
# p_estimate(Z=data$taille,
#            Zobs=55,
#            Nsim=1,
#            method="GPD",
#            estim="EMV")
# p_estimate(Z=data$taille,
#            Zobs=55,
#            Nsim=1,
#            method="BC")

calcul_p(zsim=data$taille,
         Zobs=55,
         method="GPD",
         estim="PWM")

calcul_p(zsim=data$taille,
         Zobs=55,
         method="GPD",
         estim="EMV")

calcul_p(zsim=data$taille,
         Zobs=55,
         method="BC")

### Exponential distribution ----
div<-read.csv(file = "C:/Users/Marion/Documents/divorce.csv",sep=";")[,c(2,5)] #divorces aux usa (survival data)
colnames(div)<-c("group","years")
div$averageY<-ceiling(div$years)

div2<-div %>%
  group_by(averageY)%>%
  count()

plot(div2$averageY,div2$n)
mean(div$averageY)
sd(div$averageY)

hist(rexp(nrow(div),rate=0.07))
pexp(60,rate=0.07,lower.tail=F)

calcul_p(zsim=div$averageY,
         Zobs=60,
         method="GPD",
         estim="PWM")

calcul_p(zsim=div$averageY,
         Zobs=60,
         method="GPD",
         estim="EMV")

calcul_p(zsim=div$averageY,
         Zobs=60,
         method="BC")

### Normal distribution ---
iris<-read.csv(file = "C:/Users/Marion/Documents/iris.csv",sep=";") # Fisher's iris data set
colnames(iris)<-c("sep.length","sep.width","pet.length","pet.width","species")

iris2<-iris%>%
      group_by(sep.width)%>%
      count()

plot(iris2$sep.width,iris2$n) 
mean(iris$sep.width)
sd(iris$sep.width)

hist(rnorm(nrow(iris),mean=3,sd=0.5))  #en fait normal
pnorm(4,mean=3,sd=0.5,lower.tail=F)

calcul_p(zsim=iris$sep.width,
         Zobs=4,
         method="BC")

calcul_p(zsim=iris$sep.width,
         Zobs=4,
         method="GPD",
         estim="PWM")

calcul_p(zsim=iris$sep.width,
         Zobs=4,
         method="GPD",
         estim="EMV")

### Amazone and Coca Cola Data set ---- 

amz<-read.csv(file = "C:/Users/Marion/Documents/AMZN-KO.csv",sep=";",dec=",") 
plot(amz$AMZN.Closing.Price)
plot(amz$AMZN.Daily.Percent.Return)

amz2<-amz %>%
  group_by(ceiling(AMZN.Closing.Price)) %>%
  count()

colnames(amz2)[1]<-"Price"
#amz2$freq<-amz2$n/sum(amz2$n)

plot(amz2$Price,amz2$n)
mean(amz$AMZN.Closing.Price)
sd(amz$AMZN.Closing.Price)

calcul_p(zsim=amz$AMZN.Closing.Price,
         Zobs=400,
         method="BC")

calcul_p(zsim=amz$AMZN.Closing.Price,
         Zobs=400,
         method="GPD",
         estim="PWM")

calcul_p(zsim=amz$AMZN.Closing.Price,
         Zobs=400,
         method="GPD",
         estim="EMV")

hist(rexp(nrow(amz),r=0.015))
pexp(400,rate=0.015,lower.tail=F)


### Croissance mondiale
grw<-read.csv(file = "C:/Users/Marion/Documents/growth.csv",sep=";",dec=",")
grw[is.na(grw)]<-0

## on regroupe toutes les annees dans une seule colonne 
grw2<-rep(0,nrow(grw)*(ncol(grw)-1))
for (i in 2:length(grw)) {
  j<-(i-2)*nrow(grw) + 1
  k<-(i-1)*nrow(grw)
  grw2[j:k]<-grw[,i]
}

grw2<-subset(data.frame(grw2),grw2!=0)

grw3 <- grw2 %>%
  group_by(round(grw2,digits=1)) %>%
  count()

plot(grw3$`round(grw2, digits = 1)`,grw3$n)  
mean(grw2$grw2)
sd(grw2$grw2)

hist(rpois(nrow(grw2),lambda = 8))
ppois(20,lambda=8,lower.tail=F)

calcul_p(zsim=grw2$grw2,
         Zobs=20,
         method="BC")

calcul_p(zsim=grw2$grw2,
         Zobs=20,
         method="GPD",
         estim="PWM")

calcul_p(zsim=grw2$grw2,
         Zobs=20,
         method="GPD",
         estim="EMV")
## Tracer les courbes de regression avec le Nc obtenu

Nc<-end

z<-rf(Nc,df1=5,df2=10)
z1 <- tail( sort(z) , queue )

##vraie p-valeur
pvraie<-pf(z1,df1=5,df2=10,lower.tail = F) 
plot(z1,-log(pvraie),type="l",col="red",log="xy")


##modele lineaire
p <- seq(queue,1)/Nc
ml<-PML(z1,Nc,queue)
if (critereLourd(z,Nc)==FALSE)
  lines(z1,exp(ml$inter + ml$pente*log(z1)),col="blue")
if (critereLourd(z,Nc)==TRUE)
  lines(z1,exp(ml$inter + ml$pente*log(log(z1))),col="blue")

##BOX-COX
bc<-PBC(z1,Nc,queue)
lambda<-bc$lbda
if (lambda != 0){
  if (critereLourd(z,Nc)==FALSE)
    lines(z1,exp((lambda*(bc$inter + bc$pente*log(z1)) + 1)^(1/lambda)),col="darkgreen")
  if (critereLourd(z,Nc)==TRUE)
    lines(z1,exp((lambda*(bc$inter + bc$pente*log(log(z1))) + 1)^(1/lambda)),col="darkgreen")
}
if (lambda == 0){
  if (critereLourd(z,Nc)==FALSE)
    lines(z1,exp(exp(bc$inter + bc$pente*log(z1))),col="darkgreen")
  if (critereLourd(z,Nc)==TRUE)
    lines(z1,exp(exp(bc$inter + bc$pente*log(log(z1)))),col="darkgreen")
}
legend("topleft",
       legend=c("Modele lineaire","Box-Cox","Vraie p-valeur"),
       fill=c("blue","darkgreen","red"))

## Tracer les courbes de regression avec le Nc obtenu

Nc<-1e6

z<-log(rcauchy(Nc))
z1 <- tail( sort(z) , queue )

##vraie p-valeur
pvraie<-pcauchy(exp(z1),lower.tail = F)
plot(z1,-log(pvraie),type="l",col="red",log="xy")

##modele lineaire
p <- seq(queue,1)/Nc
ml<-PML(z1,Nc,queue)
lines(z1,exp(ml$inter + ml$pente*log(z1)),col="blue")


##BOX-COX
bc<-PBC(z1,Nc,queue)
lambda<-bc$lbda
if (lambda != 0)
  lines(z1,exp((lambda*(bc$inter + bc$pente*log(z1)) + 1)^(1/lambda)),col="darkgreen")
if (lambda == 0)
  lines(z1,exp(exp(bc$inter + bc$pente*log(z1))),col="darkgreen")

legend("topleft",
       legend=c("Modele lineaire","Box-Cox","Vraie p-valeur"),
       fill=c("blue","darkgreen","red"))

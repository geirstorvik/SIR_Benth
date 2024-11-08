#Function that calculates residual sum for a given set of (beta,gamma) parameters
sir = function(beta,gamma=1,I=c(10,11),N=6000000,S=N-100)
{
  Iseq = rep(0,25)
  Sseq = rep(0,25)
  Iseq[1] = I[1]
  Sseq[1] = S[1]
  beta2 = beta/(N*24)
  gamma2 = gamma/24
  for(i in 1:24)
  {
    Sseq[i+1]=(1-beta2*Iseq[i])*Sseq[i]
    Iseq[i+1]=(1+beta2*Sseq[i]-gamma2)*Iseq[i]
  }
  res = (I[2]-Iseq[25])^2
  #show(c(beta,gamma,I[2],Iseq[25],res))
  res
}

#General settings
N = 5305127
I = 100
R = 100
S = N-I-R
M = 100
gamma.hat = seq(0.15,2,length=M)

#I(t+1)=110
beta.hat = rep(NA,M)
obj = rep(NA,M)
for(m in 1:M)
{
 foo = optimize(sir,c(0,10),gamma=gamma.hat[m],I=c(100,110),S=S,N=N)
 beta.hat[m] = foo$minimum
 obj[m] = foo$objective
}
R.hat = beta.hat/gamma.hat

#I(t+1)=90
beta.hat2 = rep(NA,M)
obj2 = rep(NA,M)
for(m in 1:M)
{
  foo = optimize(sir,c(0,10),gamma=gamma.hat[m],I=c(100,90),S=S,N=N)
  beta.hat2[m] = foo$minimum
  obj2[m] = foo$objective
}
R.hat2 = beta.hat2/gamma.hat

#Plotting results
par(mfrow=c(1,2))
matplot(cbind(gamma.hat,gamma.hat),cbind(R.hat,R.hat2),type="l",lwd=2,lty=1,col=1:2,ylim=c(0,2),
        xlab=expression(gamma[t]),ylab=expression(R[t]))
points(gamma.hat,sqrt(obj),cex=0.5)
points(gamma.hat,sqrt(obj2),cex=0.5,col=2)
abline(v=0.5,col=3)
abline(v=0.2,col=3)
#matplot(cbind(gamma.hat,gamma.hat),cbind(obj,obj2),pch=16,
#        xlab=expression(gamma),ylab="residual",col=1:2,cex=0.5)
matplot(cbind(gamma.hat,gamma.hat),cbind(beta.hat,beta.hat2),type="l",lwd=2,lty=1,col=1:2,ylim=c(0,2),
        xlab=expression(gamma[t]),ylab=expression(beta[t]))

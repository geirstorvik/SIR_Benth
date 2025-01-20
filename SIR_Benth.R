library(EpiEstim)
library(zoo)
library(ggplot2)


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
  res = abs(I[2]-Iseq[25])
  Snext <<-Sseq[25]
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


#Plotting results (figure 1)
par(mfrow=c(1,2))
matplot(cbind(gamma.hat,gamma.hat),cbind(R.hat,R.hat2),type="l",lwd=2,lty=1,col=1:2,ylim=c(0,2),
        xlab=expression(gamma[t]),ylab=expression(R[t]))
abline(v=0.5,col=3)
abline(v=0.2,col=3)
matplot(cbind(gamma.hat,gamma.hat),cbind(beta.hat,beta.hat2),type="l",lwd=2,lty=1,col=1:2,ylim=c(0,2),
        xlab=expression(gamma[t]),ylab=expression(beta[t]))


##Figure 3
library(data.table)
x = read.csv("data_covid19_lab_by_time_latest.csv")
library(ggplot2)
# Extract relevant period
xsub = x[x$date>as.Date("2021-01-31") & x$date<as.Date("2021-07-15"),]
xsub$date = as.Date(xsub$date)
xsub$ntests =  frollmean(xsub$n_pos+xsub$n_neg,7)
xsub$date = as.Date(xsub$date)
Iobs = xsub$n_pos
Imean = c(NA,NA,NA,rollmean(x$n_pos,7),NA,NA,NA)
Imean = Imean[(x$date>as.Date("2021-01-31")) & (x$date<as.Date("2021-07-15"))]

##Estimated R-curves for whole period
# First EpiEstim, two different "gamma" values
config = make_config(list(mean_si = 5.0,std_si = 0.1))
r <- EpiEstim::estimate_R(xsub$n_pos,method="parametric_si",config=config)  
xsub$rR = NA
med = apply(r$R[,1:2],1,median)
xsub$rR[med] = r$R[,3]

config = make_config(list(mean_si = 2.0,std_si = 0.1))
r <- EpiEstim::estimate_R(xsub$n_pos,method="parametric_si",config=config)  
med = apply(r$R[,1:2],1,median)
xsub$rR2 = NA
xsub$rR2[med] = r$R[,3]

# SIR-BBN with gamma=0.2
nT = length(Iobs)
It2 = Imean
St2 = rep(NA,length(Iobs))
St2[1] = N-It2[1]
R1.hat = rep(NA,nT)
SSQ1 = rep(NA,nT)
for(j in 2:nT)
{
  foo1 = optimize(sir,c(0,10),gamma=0.2,I=It2[(j-1):j],S=St2[j-1],N=N,
                  tol=.Machine$double.eps^0.5)
  R1.hat[j] = foo1$minimum/0.2
  St2[j] = Snext
  SSQ1[j] = foo1$objective
}
xsub$R1 = R1.hat

# SIR-BBN with gamma=0.5
R2.hat = rep(NA,nT)
SSQ2 = rep(NA,nT)
for(j in 2:nT)
{
  foo2 = optimize(sir,c(0,10),gamma=0.5,I=It2[(j-1):j],S=St2[j-1],N=N,
                  tol=.Machine$double.eps^0.5)
  R2.hat[j] = foo2$minimum/0.5
  St2[j] = Snext
  SSQ2[j] = foo2$objective
}
xsub$R2 = R2.hat

## Weekly averages
xsub$R1mean = rollmean(R1.hat,7,fill="extend")
xsub$R2mean = rollmean(R2.hat,7,fill="extend")
xsub$rRmean = rollmean(xsub$rR,7,fill="extend")
xsub$rR2mean = rollmean(xsub$rR2,7,fill="extend")

df = data.frame(date=c(xsub$date,xsub$date,xsub$date,xsub$date),
                R=c(xsub$R2mean,xsub$R1mean,xsub$rRmean,xsub$rR2mean),
                gamma=as.factor(c(rep(0.5,nT),rep(0.2,nT),
                                  rep("E0.2",nT),rep("E0.5",nT))))
ggplot(data=df,aes(x=date,y=R,group=gamma,colour=gamma)) + 
  geom_line(size=1.25) + ylim(0.5,1.6) +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14)) +
  labs(x="2021",y="Estimated R") +
  scale_color_manual(values = c("red","black","green","blue")) 


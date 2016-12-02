######### Authors: Alicja Wolny-Dominiak 
######### University of Economics in Katowice, mailto: woali@ue.katowice.pl
######### R version 3.2.5

# CHAPTER 4

#The list of packages
library(ChainLadder)
library(hglm)
library(statmod)
library(tweedie)
library(StatMatch) ##dummy

#Example - ata
library(ChainLadder)
cum.t <-as.triangle(auto$CommercialAutoPaid/1000)
cum.t

# Figure
par(mfrow = c(1,2))
plot(c(1:10), cum.t[1,1:10], xlab="Development year", xlim=c(1,10), ylim=c(1,130),ylab="cumulative loss", type='b')
for (j in 0:9) {
points(c(1:(10-j)), cum.t[(j+1),1:(10-j)], type='b', col=j)
}
legend("bottomright",col=c(1,2,3,4,5,6,7,8,10),lwd=2,
c("i=1", "i=2","i=3","i=4","i=5","i=6","i=7","i=8","i=9","i=10"), cex=0.8)
plot(cum.t)


ata <-ata(cum.t)
ata
d=c(1,2,3,4,5,6,7,8,9)

# Figure
plot(d,ata[1,], xlab="development", ylab="ata", col=1, type='b')
for (i in 2:9) {	
points(d, ata[i,], type='b', col=i)}
points(d,attr(ata,"vwtd"), type='l', lwd=2, col="black")
legend("topright",col=c(1,2,3,4,5,6,7,8,9),lwd=2,
c("i=1", "i=2","i=3","i=4","i=5","i=6","i=7","i=8","i=9","vwtd"), cex=0.8)

#Example - Chain-ladder Loss Reserve
library(ChainLadder)

Mortgage
MackChainLadder(Mortgage, mse.method="Mack")

#Example - loss reserving using GLM
library(ChainLadder)
cum.t <-as.triangle(auto$CommercialAutoPaid/1000)

glm.reserve <-glmReserve(cum.t, var.power = 1, link.power = 0,
cum = TRUE, mse.method = c("formula"))
glm.reserve

glm.reserve.boot <-glmReserve(cum.t, var.power = 1, link.power = 0,
cum = TRUE, mse.method = c("bootstrap"), nsim = 1000)

glm.reserve <-glmReserve(cum.t, var.power = 1, link.power = 0, cum = TRUE, mse.method = c("formula"))
glm.reserve

glm.reserve.boot <-glmReserve(cum.t, var.power = 1, link.power = 0, cum = TRUE, mse.method = c("bootstrap"), nsim = 1000)
glm.reserve.boot$coefficients[10,5]
		
#Example - Claim Reserving Using HGLM
# HGLM
library(hglm)
library(ChainLadder)
library(statmod)
library(tweedie)
library(StatMatch) ##dummy
set.seed(456)

cum.t <-as.triangle(auto$CommercialAutoPaid/1000)
cum.t

tr <-cum2incr(cum.t)

loss <- as.data.frame(tr, origin = names(dimnames(tr))[1], dev = names(dimnames(tr))[2])

triangle_upper <- subset(loss, !is.na(loss$value))
triangle_lower <- subset(loss, is.na(loss$value))

X.upper = model.matrix(~as.factor(dev),data=triangle_upper) 
Z.upper = fact2dummy(as.factor(triangle_upper$origin))
X.lower = as.matrix(cbind(rep(1,length(triangle_lower$dev)),fact2dummy(as.factor(triangle_lower$dev)))) 
Z.lower = fact2dummy(as.factor(triangle_lower$origin))

claim <- triangle_upper$value
p <-1 
beta_hglm=hglm(fixed=claim~as.factor(triangle_upper$dev),random=~1|as.factor(triangle_upper$origin), 
family=tweedie(var.power=p,link.power=0), rand.family=Gamma(log))

beta=beta_hglm$fixef
u=beta_hglm$ranef
v=log(beta_hglm$ranef)
phi <- beta_hglm$varFix

Y.u.lower <- exp(X.lower%*%beta+Z.lower%*%v[2:10])
Y.u.upper <- exp(X.upper%*%beta+Z.upper%*%v) 

Ri_hglm <- tapply(Y.u.lower, triangle_lower$origin, sum)
R_hglm <- sum(Ri_hglm)
R_hglm

summary <-rbind(as.matrix(Ri_hglm), R_hglm)
colnames(summary) <-c("HGLM Loss reserve")
summary

#Example - Bootstrap RMSEP and QAPE for HGLM Loss Reserve
boot.error <-function(claim, Y.u.upper, X.lower, Z.lower, nsim, p) 
{
Ri_hglmB <- Ri_hglmBB <- R_hglmB <- R_hglmBB <- NULL
resid=(claim-Y.u.upper)/sqrt(Y.u.upper^p) 

for (i in 1:nsim) { 
residB <- sample(resid, nrow(triangle_upper), replace = TRUE)
Y.u.upperB <- abs(residB * sqrt(Y.u.upper^p) + Y.u.upper) 
			
beta_hglmB <- hglm(fixed=Y.u.upperB~as.factor(dev),random=~1|as.factor(origin), family=tweedie(var.power=p,link.power=0), 
rand.family=Gamma(log), data=triangle_upper, maxit=10)

betaB=beta_hglmB$fixef
uB=beta_hglmB$ranef
vB=log(uB)
phiB <- beta_hglmB$varFix

Y.u.lowerB <- exp(X.lower%*%betaB+Z.lower%*%vB[2:10])
uBB = rtweedie(length(u), mu = c(u), phi =  beta_hglm$varRanef, power = 2)
vBB = log(uBB)
Y.u.lowerBB <- rtweedie(length(Y.u.lowerB), mu = c(exp(X.lower%*%beta+Z.lower%*%vB[2:10])), phi = phiB, power = p)

Ri_hglmB <- rbind(Ri_hglmB, tapply(Y.u.lowerB, triangle_lower$origin, sum))
Ri_hglmBB <- rbind(Ri_hglmBB, tapply(Y.u.lowerBB, triangle_lower$origin, sum))
R_hglmB=rbind(R_hglmB, sum(Ri_hglmB[i,]))
R_hglmBB=rbind(R_hglmBB, sum(Ri_hglmBB[i,]))
}

RMSEPi_hglm <- QAPEi_hglm <- NULL
for (k in 1:9){
RMSEPi_hglm <-rbind(RMSEPi_hglm, sqrt(sum((Ri_hglmBB[,k] - Ri_hglmB[,k])^2)/nsim))
QAPEi_hglm <- rbind(QAPEi_hglm, quantile(abs(Ri_hglmBB[,k] - Ri_hglmB[,k]), probs=c(0.5,0.75,0.9, 0.95)))
}
RMSEP_hglm <- sqrt(sum((R_hglmBB - R_hglmB)^2)/nsim)
QAPE_hglm <-quantile(abs(R_hglmBB - R_hglmB), probs=c(0.5,0.75,0.9, 0.95))
	
summary <-rbind(cbind(RMSEPi_hglm, QAPEi_hglm),c(RMSEP_hglm, QAPE_hglm))
rownames(summary) <-c("i=2", "i=3","i=4", "i=5", "i=6", "i=7", "i=8", "i=9", "i=10", "Total")
colnames(summary) <-c("RMSEP", "QAPE0.5",  "QAPE0.75", "QAPE0.9", "QAPE0.95") 
summary 
return(summary)
}

boot.er <-boot.error(claim, Y.u.upper, X.lower, Z.lower, 1000, 1) 


#Figure 
EP <-as.data.frame(boot.er[1:9,1:5])

par(mfrow=c(1,2))
plot(EP[,1],xaxt="n", type="b",ylim=c(0,20), lty=2,pch=0, 
ylab="",xlab="", main="RMSEPi_HGLM and Qi_HGLM for origin")
points(EP[,2],xaxt="n",type="b",lty=2,pch=1, ylab="",xlab="")
points(EP[,3],xaxt="n",type="b", lty=2,pch=2, ylab="",xlab="")
points(EP[,4],xaxt="n",type="b", lty=2,pch=3, ylab="",xlab="")
legend("topleft",c("RMSEPi","Qi_0.5","Qi_0.75","Qi_0.9"), 
lty=c(2,2,2,2),pch = c(0,1,2,3), cex=0.75, inset=.05, bty = "n")
labele <-c("i=2","i=3","i=4","i=5","i=6","i=7","i=8", "i=9", "i=10")
mtext(labele,at=1:9,side=1)

plot(c(RMSEP_hglm, QAPE_hglm),xaxt="n",type="b",ylab="",
xlab="", main="RMSEP and QAPE _HGLM")
labele = c("RMSEP","Q_0.5","Q_0.75","Q_0.9")
mtext(labele,at=1:4,side=1)



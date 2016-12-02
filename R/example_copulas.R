######### Authors: Stanis³aw Wanat, Alicja Wolny-Dominiak 
######### University of Economics in Katowice, mailto: woali@ue.katowice.pl
######### R version 3.2.5

# CHAPTER 3

#The list of packages
library(copula)
library(statmod)
library(fitdistrplus)
library(CopulaRegression)
library(ggplot2)
library(plotrix)
library(cplm)

#3.1.1.1 Example
library(copula)
sigma <-0.3
kopula.norm <- normalCopula(param = sigma, dim = 2)
par(mfrow=c(2,2))
contour(kopula.norm, dCopula,  nlevels=50)
contour(kopula.norm, pCopula,  nlevels=20)
persp(kopula.norm,dCopula, xlab="u1", ylab="u2", zlab="c(u1, u2)")
persp(kopula.norm,pCopula, xlab="u1", ylab="u2", zlab="C(u1, u2)")

#3.1.1.2 Example
library(statmod)
mvd<- mvdc(kopula.norm, c("invgauss", "norm"),
list(list(mean=350,dispersion=4.5), list(mean = 0, sd =0.2)))

par(mfrow=c(2,2))
persp(mvd, dMvdc, col=4, xlim = c(0, 2), ylim=c(-0.5, 0.5),
xlab="x1", ylab="x2", zlab="f(x1, x2)")
persp(mvd, pMvdc, col=4, xlim = c(0, 2), ylim=c(-0.5, 0.5),
xlab="x1", ylab="x2", zlab="f(x1, x2)")
contour(mvd, dMvdc, col=4, xlim = c(0, 2), ylim=c(-0.5, 0.5), 
xlab="x1", ylab="x2",nlevels=30)
contour(mvd, pMvdc, col=4, xlim = c(0, 2), ylim=c(-0.5, 0.5), 
xlab="x1", ylab="x2",nlevels=30)

#3.1.2.1 Example
library(fitdistrplus)

sim <- rMvdc(1000,mvd)

fit.X1<-fitdist(sim[,1], dinvgauss, start= list(mean=350,dispersion=4.5), method="mle")
(mean.invgauss<-coef(fit.X1)[1]) 
(dispersion.invgauss<-coef(fit.X1)[2])
#plot(fit.X1)

fit.X2<-fitdist(sim[,2], "norm", method="mle")
(mean.norm<-coef(fit.X2)[1]) 
(sd.norm<-coef(fit.X2)[2])

u.1<-pinvgauss(sim[,1],mean=mean.invgauss,dispersion=dispersion.invgauss)
u.2<-pnorm(sim[,2], mean=mean.norm, sd=sd.norm)
U<-cbind(u.1, u.2)

#Figure
par(mfrow=c(1,2))
plot(sim, xlim = c(0, 2), ylim=c(-0.5, 0.5),xlab=expression(x[1]), ylab=expression(x[2]))
plot(U, xlim = c(0, 1), ylim=c(0, 1) ,xlab=expression(u[1]), ylab=expression(u[2]))

fit.normal<- fitCopula(normalCopula(dim = 2),U, method="ml")
summary(fit.normal)
(sigma.fit<-coef(fit.normal))

gofCopula(normalCopula(sigma.fit), sim)

#3.2.1.1 Example - Fitting Frequency-severity Parameters
library(CopulaRegression)
library(ggplot2)

auto3 <- read.csv2("http://web.ue.katowice.pl/woali/auto_data3.csv", header=TRUE, sep=";")
attach(auto3)

#Figure
par(mfrow=c(1,2))
hist(N, xlab="",main="Number of claims")
plot(density(Y.N/1000), xlab="",main="Average value of claims (in thou.)", axes = F)
axis(1)
axis(2) 

#Figure 
desc.stat <- function(x)
{out <- c(min(x), mean(x) - sd(x), mean(x), mean(x) + sd(x), max(x))
 names(out) <- c("ymin", "lower", "middle", "upper", "ymax")
 out}
 
plot.1 <- ggplot(aes(y = Y.N, x = as.factor(X1)), data = auto3)
plot.11 <- plot.1 + stat_summary(fun.data = desc.stat, geom = "boxplot")+ 
geom_jitter(position=position_jitter(width=.2), size=0.5)+ ggtitle("")+ xlab("X1") + ylab("")+ theme_bw()

plot.2 <- ggplot(aes(y = Y.N, x = as.factor(X2)), data = auto3)
plot.12 <- plot.2 + stat_summary(fun.data = desc.stat, geom = "boxplot")+ 
geom_jitter(position=position_jitter(width=.2), size=0.5)+ ggtitle("")+ xlab("X2") + ylab("")+ theme_bw()

plot.3 <- ggplot(aes(y = Y.N, x = as.factor(X3)), data = auto3)
plot.13 <- plot.3 + stat_summary(fun.data = desc.stat, geom = "boxplot")+ 
geom_jitter(position=position_jitter(width=.2), size=0.5)+ ggtitle("")+ xlab("X3") + ylab("")+ theme_bw()

plot.4 <- ggplot(aes(y = Y.N, x = as.factor(X4)), data = auto3)
plot.14 <- plot.4 + stat_summary(fun.data = desc.stat, geom = "boxplot")+ 
geom_jitter(position=position_jitter(width=.2), size=0.5)+ ggtitle("")+ xlab("X4") + ylab("")+ theme_bw()


f <- grid.arrange(plot.11, plot.12, plot.13, plot.14, nrow=2, ncol=2)
ggsave("boxplotX1.jpg", f)

#Model
R=model.matrix(~as.factor(X1)+X2+X3+as.factor(X4), auto3)
S.bez=model.matrix(~1,auto3)
m1=copreg(Y.N,N,R,S.bez,family=1, E, sd.error=TRUE, zt=TRUE)
m3=copreg(Y.N,N,R,S.bez,family=3, E, sd.error=TRUE, zt=TRUE)
m4=copreg(Y.N,N,R,S.bez,family=4, E, sd.error=TRUE, zt=TRUE)
m5=copreg(Y.N,N,R,S.bez,family=5, E, sd.error=TRUE, zt=TRUE)

summary <- cbind(rbind(as.matrix(m1$alpha), m1$delta, m1$beta, m1$theta), 
rbind(as.matrix(m3$alpha), m3$delta, m3$beta, m3$theta),
rbind(as.matrix(m4$alpha), m4$delta, m4$beta, m4$theta),
rbind(as.matrix(m5$alpha), m5$delta, m5$beta, m5$theta))
colnames(summary) <-c("Gauss", "Clayton", "Gumbel", "Frank")
rownames(summary) <-c("R(Intercept)", "Ras.factor(X1)1", "RX224-27", 
"RX228-44", "RX245-57", "RX258-75", "RX276+", "RX3DIE", "Ras.factor(X4)1", "delta", "beta", "theta")
summary

m1$loglik
m3$loglik
m4$loglik
m5$loglik
m3$theta
tau <-BiCopPar2Tau(par = m3$theta, family = 3)
tau

hat.Y.N <-exp(R%*%m3$alpha)
hat.Y.N[1:10]
hat.lambda <-E*exp(S.bez%*%m3$beta)
hat.N<-hat.lambda/(1-exp(-hat.lambda))
hat.N[1:10]

predicted<-predict(m3,R,S.bez, exposure=E)
predicted$x.pred[1:10]
predicted$y.pred[1:10]

#Example - the joint and conditional densities of the variable $(Y_i,N_i)$

x<-seq(from=20, to =10000, by=20)
fxy1<-fxy2<-fxy3<-rep(0,500)
for (i in 1:500)
{fxy1[i]<-density_joint(x[i],1,mean(Y.N),m3$delta,hat.lambda[1],
m3$theta,family=3, zt=TRUE)/dztp(1, hat.lambda[1])
fxy2[i]<-density_joint(x[i],2,mean(Y.N),m3$delta,hat.lambda[1],
m3$theta,family=3, zt=TRUE)/dztp(2, hat.lambda[1])
fxy3[i]<-density_joint(x[i],3,mean(Y.N),m3$delta,hat.lambda[1],
m3$theta,family=3, zt=TRUE)/dztp(3, hat.lambda[1])
}

fxy1
fxy2
fxy3

#Figure
plot(x,fxy1, col=1, ylab="f(Y|N=1)", type="l",
xlim = c(0,2000), ylim=c(0, 0.0003), lty=1, lwd=2, axes = F)
axis(1)
axis(2)
lines(x,fxy2, col=2, ylab="f(Y|N=2)", type="l",
xlim = c(0,2000), ylim=c(0, 0.0003), lty=1, lwd=2)
lines(x,fxy3, col=3, ylab="f(Y|N=3)", type="l",
xlim = c(0,2000), ylim=c(0, 0.0003), lty=1, lwd=2)
legend("topright",col=c(1,2,3),lwd=2,c("f(Y|N=1)",
"f(Y|N=2)","f(Y|N=3)"), lty=1:3)

#3.2.1.3 Example - Copula-based distribution of the aggregate value of claims for the whole portfolio
library(CopulaRegression)
library(plotrix)

auto3 <- read.csv2("http://web.ue.katowice.pl/woali/auto_data3.csv", header=TRUE, sep=";", dec=",")
attach(auto3)

S.bez=model.matrix(~1,auto3)
fit.0<- function(cop)
{model<-copreg(Y.N,N,S.bez,S.bez,family=cop,E, sd.error=TRUE, joint=TRUE, zt=TRUE)}

fit.Gauss<-fit.0(1)
fit.Clayton<-fit.0(3)
fit.Gumbel<-fit.0(4)
fit.Frank<-fit.0(5)

# Model parameters for different copulas
param<-function(model)
{out<-as.matrix(c(model$alpha[[1]],model$sd.alpha[[1]],model$beta[[1]], model$sd.beta[[1]],model$delta[[1]],model$theta[[1]],
	model$tau[[1]],-2*model$loglik[[1]]+2*model$npar))
	out}
	
summary<-cbind(param(fit.Gauss), param(fit.Clayton),param(fit.Gumbel), param(fit.Frank))
par<-c("alpha.0", "sd.alpha.0","beta.0","sd.beta.0", "phi", "theta", "tau", "AIC")
cop<-c("Gauss", "Clayton", "Gumbel", "Frank")
rownames(summary)<-par
colnames(summary)<-cop
round(summary,4)

total.loss.sim<-function(model)
	{mu<-exp(model$alpha)
	lambda<-exp(model$beta)
	rand <-simulate_joint(322, mu, model$delta, lambda, model$theta, family=model$family, max.y = 500, eps = 1e-05)
	loss<-rand[,1]*rand[,2]
	loss<-sum(loss)
	loss}

n.boot<-500
total.Gauss<-replicate(n.boot, total.loss.sim(fit.Gauss))
total.Clayton<-replicate(n.boot, total.loss.sim(fit.Clayton))
total.Gumbel<-replicate(n.boot, total.loss.sim(fit.Gumbel))
total.Frank<-replicate(n.boot, total.loss.sim(fit.Frank))

#Figure
par(mfrow=c(1,1))
plot(density(total.Gauss), lty=1, col=1, ylim=c(0, 4.0e-06),main="Total loss densinty")
lines(density(total.Clayton), lty=2,col=2)
lines(density(total.Gumbel), lty=3, col=3)
lines(density(total.Frank), lty=4,col=4)
legend("topright",lty=1:4,c("Gauss", "Clayton", "Gumbel", "Frank"), col=1:4)


#3.2.1.2 Example - Plotting Conditional Probability Mass Function of the Number of Claims
library(CopulaRegression)
library(plotrix)
library(copula)

condit.prob <- function(tau)
{theta<-iTau(normalCopula(), tau)
out<-density_conditional(y=1:10,x=1000,mu=900,delta=0.1, lambda=3, theta=theta,family=1)}

prob<-as.data.frame(rbind(condit.prob(0),condit.prob(0.1),condit.prob(0.3), condit.prob(0.5)))
names(prob)=1:10
tau<-c(expression(paste(tau,"=0.0")),expression(paste(tau,"=0.1")),expression(paste(tau,"=0.3")),expression(paste(tau,"=0.5")))

#Figure
barp(prob, col=c("grey30","grey50","grey70","grey90"),names.arg=names(prob),legend.lab=tau, legend.pos=list(x=9,y=0.25))


theta<-rep(0,5)
theta[1]<-iTau(normalCopula(), 0.3)
theta[3]<-iTau(claytonCopula(), 0.3)
theta[4]<-iTau(gumbelCopula(), 0.3)
theta[5]<-iTau(frankCopula(), 0.3)
condit.prob <- function(fam){
out<-density_conditional(y=1:10,x=1000,mu=900,delta=0.1, lambda=3,theta=theta[fam],family=fam)}

density <-as.data.frame(rbind(condit.prob(1),condit.prob(3),condit.prob(4), condit.prob(5)))
names(density)=1:10
copula <-c("Gauss","Clayton", "Gumbel","Frank")

#Figure
barp(density, col=c("grey30","grey50","grey70","grey90"), names.arg=names(density),legend.lab=copula,legend.pos=list(x=9,y=0.25))

#3.3.1.1 Example - Pure Premium Calculation Using Compound Poisson Generalized Linear Models
library(cplm)
library(ggplot2)

data <- read.csv2("http://web.ue.katowice.pl/woali/brazdata.csv", header=TRUE, sep=";", dec=",")
attach(data)
data[11:20, ]

group<-interaction(Gender,DrivAge,VehAge )
g.1<-g.2<-group
for (j in 1:30) {g.2<-ifelse(g.1==levels(group)[j],j,g.2)}
group.nb<-as.factor(g.2)
legend<-data.frame(Number=levels(group.nb),Group=levels(group))
legend
lob<-data.frame(Group=group.nb,lob.1=ClaimAmountRob,lob.2=ClaimAmountPartColl,
lob.3=ClaimAmountTotColl, lob.4=ClaimAmountOther)
lob[11:20,]

### Figures
par(mfrow=c(2,2))
hist(lob$lob.1[lob$lob.1<1000], xlim=c(0,1000), main="LOB 1")
hist(lob$lob.2[lob$lob.2<1000], xlim=c(0,1000), main="LOB 2")
hist(lob$lob.3[lob$lob.1<1000], xlim=c(0,1000), main="LOB 3")
hist(lob$lob.4[lob$lob.1<100], xlim=c(0,1000), main="LOB 4")

par(mfrow=c(3,2))
plot(lob$lob.1, lob$lob.2, main="LOB 1 and LOB 2")
plot(lob$lob.1, lob$lob.3, main="LOB 1 and LOB 3")
plot(lob$lob.1, lob$lob.4, main="LOB 1 and LOB 4")
plot(lob$lob.2, lob$lob.3, main="LOB 2 and LOB 3")
plot(lob$lob.2, lob$lob.4, main="LOB 2 and LOB 4")
plot(lob$lob.3, lob$lob.4, main="LOB 3 and LOB 4")

avr.dev.max <- function(x) {
  wyn <- c(0, 0, mean(x), mean(x) + sd(x), max(x))
  names(wyn) <- c("ymin", "lower", "middle", "upper", "ymax")
  wyn}


lob.1.p1 <- ggplot(aes(y = lob.1, x = Group), data = lob)
lob.1.p12 <- lob.1.p1 + stat_summary(fun.data =avr.dev.max, geom = "boxplot",col=2)+ geom_jitter(position=position_jitter(width=.4), size=0.5)+
ggtitle("LOB 1")+xlab("Group") + ylab("")
lob.1.p12

lob.2.p1 <- ggplot(aes(y = lob.2, x = Group), data = lob)
lob.2.p12 <- lob.2.p1 + stat_summary(fun.data =avr.dev.max, geom = "boxplot",col=2)+ geom_jitter(position=position_jitter(width=.4), size=0.5)+
ggtitle("LOB 2")+xlab("Group") + ylab("")
lob.2.p12

lob.3.p1 <- ggplot(aes(y = lob.3, x = Group), data = lob)
lob.3.p12 <- lob.3.p1 + stat_summary(fun.data =avr.dev.max, geom = "boxplot",col=2)+ geom_jitter(position=position_jitter(width=.4), size=0.5)+
ggtitle("LOB 3")+xlab("Group") + ylab("")
lob.3.p12

lob.4.p1 <- ggplot(aes(y = lob.4, x = Group), data = lob)
lob.4.p12 <- lob.4.p1 + stat_summary(fun.data =avr.dev.max, geom = "boxplot",col=2)+ geom_jitter(position=position_jitter(width=.4), size=0.5)+
ggtitle("LOB 4")+xlab("Group") + ylab("")
lob.4.p12

### Estimation of marginal model parameters
model.lob.1<-cpglm(lob$lob.1~Gender+DrivAge+VehAge, offset=NULL, link="log")
summary(model.lob.1)

model.lob.2<-cpglm(lob$lob.2~Gender+DrivAge+VehAge, offset=NULL, link="log")
summary(model.lob.2)

model.lob.3<-cpglm(lob$lob.3~Gender+DrivAge+VehAge, offset=NULL, link="log")
summary(model.lob.3)

model.lob.4<-cpglm(lob$lob.4~Gender+DrivAge+VehAge, offset=NULL, link="log",optimizer="bobyqa" )
summary(model.lob.4)

mu.lob<-as.matrix(cbind(fitted(model.lob.1),fitted(model.lob.2),fitted(model.lob.3),fitted(model.lob.4)))
mu.lob[1:10,]

phi.lob<-c(model.lob.1$phi,model.lob.2$phi,model.lob.3$phi,model.lob.4$phi)
phi.lob

p.lob<-c(model.lob.1$p,model.lob.2$p,model.lob.3$p,model.lob.4$p)
p.lob

### Pure premium 
pure.prem.lob<-sapply(1:4, function(i) mu.lob[,i] )

par(mfrow=c(2,2))
hist(pure.prem.lob[,1],  xlab=" ",main="LOB1")
hist(pure.prem.lob[,2],  xlab=" ",main="LOB2")
hist(pure.prem.lob[,3],  xlab=" ",main="LOB3")
hist(pure.prem.lob[,4],  xlab=" ",main="LOB4")
summary(pure.prem.lob)

#Figure
par(mfrow=c(2,2))
plot(pure.prem.lob[,1]~lob$Group, main="LOB1")
plot(pure.prem.lob[,2]~lob$Group, main="LOB2")
plot(pure.prem.lob[,3]~lob$Group, main="LOB3")
plot(pure.prem.lob[,4]~lob$Group, main="LOB4")

n.policy<-length(lob$lob.1)
pure.prem<-sapply(1:n.policy, function(i) sum(pure.prem.lob[i,]))

#Figure
par(mfrow=c(1,2))
hist(pure.prem,  xlab=" ",main="Total LOB")
plot(pure.prem~lob$Group,xlab="Group", ylab="Pure premium", main="Total LOB")

#3.3.2.1 Example - Quantile premium in groups - comonotonic
library(tweedie)
quant.prem.lob<-function(lob){
	quant.prem<-matrix(0, nrow=n.policy, ncol=3)
	for (i in 1:n.policy){quant.prem[i,]<-qtweedie(c(0.95, 0.98, 0.99), power=p.lob[lob], mu=mu.lob[i,lob], phi=phi.lob[lob])}
	quant.prem}
	
quant.prem.lob.1<-quant.prem.lob(1)
quant.prem.lob.2<-quant.prem.lob(2)
quant.prem.lob.3<-quant.prem.lob(3)
quant.prem.lob.4<-quant.prem.lob(4)

quant.prem.comon<-quant.prem.lob.1+quant.prem.lob.2+quant.prem.lob.3+quant.prem.lob.4
quant.prem.comon.regr<-data.frame("Q(0.95)"=quant.prem.comon[,1], "Q(0.98)"=quant.prem.comon[,2], "Q(0.99)"=quant.prem.comon[,3],Gender=Gender,DrivAge=DrivAge,VehAge=VehAge)
quant.prem.comon.regr[1:10,]

sumary.quant.prem.comon<-sapply(1:3, function(i) summary(quant.prem.comon[,i]))
colnames(sumary.quant.prem.comon)<- c("Q(95)","Q(98)","Q(99)")
sumary.quant.prem.comon

#Figure
par(mfrow=c(3,1))
boxplot(quant.prem.comon.regr$Q.0.95.~quant.prem.comon.regr$Group, ylim=c(0,10000), main="Q(0.95)")
boxplot(quant.prem.comon.regr$Q.0.98.~quant.prem.comon.regr$Group, ylim=c(0,10000), main="Q(0.98)")
boxplot(quant.prem.comon.regr$Q.0.99.~quant.prem.comon.regr$Group, ylim=c(0,10000), main="Q(0.99)")

#3.3.3.1 Example - Copula-based Quantile Premium for Dependent LOBs
library(cplm)
library(tweedie)
library(copula)
library(ggplot2)
library(plotrix)

data <- read.csv2("http://web.ue.katowice.pl/woali/brazdata.csv", header=TRUE, sep=";")
out<-data[data$ClaimAmountOther>2000,][1]
data<-data[-out[,1],]
attach(data)
lob<-data.frame(lob.1=ClaimAmountRob,lob.2=ClaimAmountPartColl,lob.3=ClaimAmountTotColl,lob.4=ClaimAmountOther)
attach(lob)
lob[11:20,]

#Margin parameter estimation
model.lob.1<-cpglm(lob.1~1)
model.lob.2<-cpglm(lob.2~1)
model.lob.3<-cpglm(lob.3~1, optimizer="bobyqa")
model.lob.4<-cpglm(lob.4~1)

phi.lob <-c(model.lob.1$phi,model.lob.2$phi, model.lob.3$phi,model.lob.4$phi)
p.lob <-c(model.lob.1$p,model.lob.2$p,model.lob.3$p,model.lob.4$p)
mu.lob<-c(exp(model.lob.1$coefficients),exp(model.lob.2$coefficients), exp(model.lob.3$coefficients),exp(model.lob.4$coefficients))

phi.lob
p.lob
mu.lob

# Copula parameter estimation
U <-sapply(1:4, function(i) ptweedie(lob[,i], power=p.lob[i],mu=mu.lob[i], phi=phi.lob[i]))

#Gauss copula
#Start parameters
X<-pobs(lob)
normal.cop<-normalCopula(c(0,0,0,0,0,0),dim = 4, dispstr = "un")
fit.normal.mpl<- fitCopula(normal.cop, X, method="mpl")
rho.start<-sapply(1:6, function(i) coefficients(fit.normal.mpl)[i])
rho.start

# IFM estimation
fit.normal<- fitCopula(normal.cop, U, start = rho.start,method="ml" )
summary(fit.normal)
par.normal<-sapply(1:6, function(i) coefficients(fit.normal)[i])

# Student copula
# Start parameters
student.cop<-tCopula(c(0,0,0,0,0,0),dim = 4, dispstr = "un")
fit.student.mpl<- fitCopula(student.cop, X, method="mpl")
start<-sapply(1:7, function(i) coefficients(fit.student.mpl)[i])
start

# IFM estimation
fit.student<- fitCopula(student.cop, U, start=start, method="ml")
summary(fit.student)
par.student<-sapply(1:6, function(i) coefficients(fit.student)[i])
df<-coefficients(fit.student)[7]

###Simulation
sim <-function(n.sim, copula)
{rand <- rCopula(n.sim, copula)
rand.lob<-sapply(1:4, function(i) qtweedie(rand[,i], power=p.lob[i],mu=mu.lob[i], phi=phi.lob[i]))
sim <-rand.lob[,1]+rand.lob[,2]+rand.lob[,3]+rand.lob[,4]
return(as.vector(sim))
}
copula.normal<-normalCopula(par.normal,dim = 4, dispstr = "un")
copula.student<-tCopula(par.student,dim = 4,df=df, dispstr = "un")
n.sim<-1000
n.boot<-200
sim.total.normal<-replicate(n.boot, sim(n.sim, copula.normal))
sim.total.student<-replicate(n.boot, sim(n.sim, copula.student))

#Simulation Results
#sim.total.normal<-read.table("http://web.ue.katowice.pl/woali/sim_Gauss")
#sim.total.student <- read.table("http://web.ue.katowice.pl/woali/sim_Student")

#Quantile premium
quant.normal <-apply(sim.total.normal,2, quantile,probs = c(0.75, 0.90, 0.95, 0.99))
quant.student <-apply(sim.total.student,2, quantile,probs = c(0.75, 0.90, 0.95, 0.99))
(quant.prem.normal <-apply(quant.normal,1,mean))
(quant.prem.student <-apply(quant.student,1,mean))

#Figure
par(mfrow=c(2,2))
plot(density(quant.student[1,]), lty=1,col=1,main="Quantile premium - Q(0.75)")
lines(density(quant.normal[1,]), lty=2, col=2)
legend("topleft",lty=1:2,c("Student", "Gauss"), col=1:2)
plot(density(quant.student[2,]), lty=1,col=1,main="Quantile premium - Q(0.90)")
lines(density(quant.normal[2,]), lty=2, col=2)
legend("topleft",lty=1:2,c("Student", "Gauss"), col=1:2)
plot(density(quant.student[3,]), lty=1,col=1,main="Quantile premium - Q(0.95)")
lines(density(quant.normal[3,]), lty=2, col=2)
legend("topleft",lty=1:2,c("Student", "Gauss"), col=1:2)
plot(density(quant.student[4,]), lty=1,col=1,main="Quantile premium - Q(0.99)")
lines(density(quant.normal[4,]), lty=2, col=2)
legend("topleft",lty=1:2,c("Student", "Gauss"), col=1:2)

#Figure
par(mfrow=c(1,1))
dist.S.normal<-apply(sim.total.normal,1,mean)
dist.S.student<-apply(sim.total.student,1,mean)
plot(density(dist.S.student), lty=1,col=1,ylim=c(0,0.0045),main="Distribution of S")
lines(density(dist.S.normal), lty=2, col=2)
legend("topleft",lty=1:2,c("Student", "Gauss"), col=1:2)

#Comparison 
quant.prem.co.lob <-sapply(1:4, function(i) qtweedie(c(0.75, 0.90, 0.95, 0.99), power=p.lob[i],mu=mu.lob[i], phi=phi.lob[i]))
quant.prem.co <-apply(quant.prem.co.lob,1,sum)
compar<-rbind(quant.prem.co,quant.prem.normal, quant.prem.student)
colnames(compar)<- c("Prem.Q(75)","Prem.Q(90)","Prem.Q(95)","Prem.Q(99)")
rownames(compar)<- c("Comonotonicity","Gauss copula","Student copula")
compar

#Fugure
par(mfrow=c(1,1))
barp(premium, col=c("grey10","grey50","grey90"), names.arg=c("Q(75)", "Q(90)", "Q(95)","Q(99)"), 
legend.lab=rownames(premium),   legend.pos=list(x=1,y=7000), main="Quantile premium")


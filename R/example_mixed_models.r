######### Author: Alicja Wolny-Dominiak, Daniel Sobiecki 
######### University of Economics in Katowice, mailto: woali@ue.katowice.pl
######### R version 3.2.5

####CHAPTER 2
##List of packages 
library(lme4)
library(MASS)
library(glmmML)
library(faraway)
library(insuranceData)
library(MixedPoisson)
library(hglm)
library(StatMatch)
library(lattice)
library(pls)
library(stats)
library(cplm)
library(pscl)
library(graphics)
library(gridExtra)

#Example - fitting the GLMM 
library(lme4)
library(MASS)
library(glmmML)
library(faraway)
attach(motorins)

glmm1 <-glmer(Claims~Kilometres+(1|Bonus), family=poisson(link="log"), offset=log(Insured), data=motorins)
glmm2 <-glmmPQL(fixed=Claims~Kilometres+offset(log(Insured)), random=~1|Bonus, family=poisson(link="log"), data=motorins)
glmm3 <-glmmML(Claims~Kilometres, cluster=Bonus, offset=log(Insured), method = c("Laplace"), family=poisson(link="log"), data=motorins)

##Summary
summ <-rbind(cbind(glmm1@beta, glmm2$coefficients$fixed, glmm3$coefficients),
cbind(glmm1@theta, glmm2$sigma, glmm3$sigma))
col <-c("AG-H quadrature","PQL", "Laplace")
row <-c("Intercept", "K2", "K3", "K4", "K5", "Dispersion")
rownames(summ)<-row
colnames(summ)<-col
summ

##random effects
cbind("glmer"=glmm1@u, "glmmPQL"=glmm2$coefficients$random$Bonus)

#Example - Fitting HGLM Inverse-Gaussian
library(hglm)
library(insuranceData)
data(dataCar)
attach(dataCar)
names(dataCar)

Y <-claimcst0[claimcst0>0]
X <-as.factor(veh_age)[claimcst0>0]
Z <-as.factor(veh_body)[claimcst0>0]

hglm.ig <- hglm(fixed=Y~X, 
random=~1|veh_body[claimcst0>0], vcovmat = TRUE,
family=inverse.gaussian(log), disp.family=Gamma(log),
data=dataCar)

print(summary(hglm.ig), print.ranef = TRUE)


#Example - Estimation of Parameters of NB Distribution
library(insuranceData)
library(MASS)
data(SingaporeAuto)

claims<-SingaporeAuto$Clm_Count
fitdistr(claims, "Negative Binomial") #MLEs with BFGS method

library(fitdistrplus) #requires R >= 3.2.0
fitdist(claims, "nbinom",method = "mle",optim.method="Nelder-Mead")
fitdist(claims, "nbinom",method = "mle",optim.method="BFGS")
fitdist(claims, "nbinom",method = "mle",optim.method="SANN")
fitdist(claims, "nbinom",method = "mme") # Moment matching estimation

library(vcd)
out <-goodfit(claims, type =  "nbinomial", method =  "MinChisq") #ML or Minchisq method
summary(out)
out$par
out

#Example - Estimation of Parameters of Mixed Poisson Distribution
library(MixedPoisson2)
library(Rmpfr)
library(gaussquad)

variable <-rpois(322,2)
X=as.matrix(rep(1, length(variable)))

EM_lognormal <-function(variable, X, n.iter.max){
lambda1=c(); nu1=c(); beta1=c()
lambda.old = c(); nuold = c()
lambda.new = c(); nu.new = c()
lambda.old = lambda_start(variable, X)$lambda
nu.old = 1
n.iter=0

for (i in 1:n.iter.max){
### E-Step
t = pseudo_values(variable, mixing=c("lognorm"), 
lambda=lambda.old, nu=nu.old, n=100)

### M-Step
lam <-lambda_m_step(variable, X, t$pseudo_values)
lambda.new <-lam$lambda
nu.new = est.nu(t=t$pseudo_values)$nu
nu.old = nu.new 
lambda.old = lambda.new
n.iter = n.iter + 1
lambda1=c(lambda1,lambda.new[1])
nu1=c(nu1,nu.new)
beta1=c(beta1,lam$beta)
}

tab <-rbind(lambda1, beta1, nu1)
rownames(tab)=c("lambda", "beta_0", "nu")
colnames(tab)=NULL
tab

outlist = list(tab)
return(tab)

}

EM_lognormal(variable, X, n.iter.max=4)

###Mixed Models in Ratemaking

#Example - Frequency-severity HGLM and Pure Risk Premium
library(hglm)
library(StatMatch)
library(lattice)
library(pls)
library(stats)
auto5 <- read.csv2("http://web.ue.katowice.pl/woali/auto5.csv", header=TRUE, sep=";", dec=",")

attach(auto5)
nz <-CLAIM_AMOUNT>0

meanY=data.frame(CLAIM_AMOUNT[nz], GENDER[nz], ENGINE[nz])
meanN=data.frame(EXPOSURE[nz]*CLAIM_COUNT[nz], GENDER[nz], ENGINE[nz])

par(mfrow=c(2,2))
plot.design(meanY, fun=mean, xlab="The average value of claims", ylab="Mean", cex.lab=1.5, cex.axis=1.3)
plot.design(meanN, fun=mean, xlab="The number of claims", ylab="Mean", cex.lab=1.5, cex.axis=1.3) 
plot.design(meanY, fun=sd, xlab="The average value of claims",ylab="Standard deviation", cex.lab=1.5, cex.axis=1.3) 
plot.design(meanN, fun=sd, xlab="The number of claims", ylab="Standard deviation", cex.lab=1.5, cex.axis=1.3) 

group <-interaction(GENDER,ENGINE)
data.frame(levels(group))

##cv.error
S=CLAIM_AMOUNT[nz]
X1= GENDER[nz]
X2= ENGINE[nz]
N= CLAIM_COUNT[nz]

dataset <-data.frame(S, X1, X2, N)
attach(dataset)

	mse.cv.glm.gamma=function(K){		
	cvseg=c()
	set.seed(113) 
	cvseg=cvsegments(nrow(dataset), K) 
	ModelMSE=c()
	
		for (i in 1:K) {
			validset=NULL
			validset=eval(parse(text=paste("cvseg$V", i, sep="")))
			datasetTrainCV= NULL; datasetValidCV= NULL
			datasetTrainCV= dataset[-validset,]
			datasetValidCV= dataset[validset,]
			Model=NULL
			Model=glm(S~X1+X2, weight=N, family=Gamma(log),  data=datasetTrainCV)
			pred.valid=NULL
			pred.valid=exp(predict(Model, datasetValidCV))
			MSE.valid=NULL  
			MSE.valid=sum((datasetValidCV$S-pred.valid)^2)/length(pred.valid)
			ModelMSE=c(ModelMSE, MSE.valid)
			
			}
	MSE.CV=NULL
	MSE.CV=mean(ModelMSE)
	
outlist <- list(ModelMSE=ModelMSE, ModelRMSE=sqrt(ModelMSE), MSE.CV=MSE.CV, RMSE.CV=sqrt(MSE.CV))
return(outlist)
	}

	mse.cv.glm.ig=function(K){		
	cvseg=c()
	set.seed(113) 
	cvseg=cvsegments(nrow(dataset), K) 
	ModelMSE=c()
	
		for (i in 1:K) {
			validset=NULL
			validset=eval(parse(text=paste("cvseg$V", i, sep="")))
			datasetTrainCV= NULL; datasetValidCV= NULL
			datasetTrainCV= dataset[-validset,]
			datasetValidCV= dataset[validset,]
			Model=NULL
			Model=glm(S~X1+X2, weight=N, family=inverse.gaussian(log),  data=datasetTrainCV)
			pred.valid=NULL
			pred.valid=exp(predict(Model, datasetValidCV))
			MSE.valid=NULL  
			MSE.valid=sum((datasetValidCV$S-pred.valid)^2)/length(pred.valid)
			ModelMSE=c(ModelMSE, MSE.valid)
			
			}
	MSE.CV=NULL
	MSE.CV=mean(ModelMSE)
	
outlist <- list(ModelMSE=ModelMSE, ModelRMSE=sqrt(ModelMSE), MSE.CV=MSE.CV, RMSE.CV=sqrt(MSE.CV))
return(outlist)
	}

mse.cv.glm.gamma(K=5)
mse.cv.glm.ig(K=5)

##HGLM Gamma
hglm.gamma.gamma <-hglm(fixed=CLAIM_AMOUNT[nz]~GENDER[nz]+ENGINE[nz], family=Gamma(log), random=~1|CAR_MAKE[nz], rand.family=Gamma(log), data=auto5) 
#summary(hglm.gamma.gamma)
hglm.poisson <-hglm(fixed=CLAIM_COUNT~GENDER+ENGINE, family=poisson(log), random=~1|CAR_MAKE, rand.family=Gamma(log), data=auto5) 
#summary(hglm.poisson)

beta.Y <-data.frame("beta.Y"=hglm.gamma.gamma$fixef, "Tarrif rates"=exp(hglm.gamma.gamma$fixef))
phi.Y <-hglm.gamma.gamma$varFix
phiu.Y <-hglm.gamma.gamma$varRanef
beta.N <-data.frame("beta.N"=hglm.poisson$fixef, "Tarrif rates"=exp(hglm.poisson$fixef))
phi.N <-hglm.poisson$varFix
phiu.N <-hglm.poisson$varRanef

##summary
cbind(round(beta.Y,2), round(beta.N,2))
data.frame(phi.Y, phiu.Y, phi.N,phiu.N)

ranef.Y <-data.frame("vu.Y"=exp(hglm.gamma.gamma$ranef), "u.Y"=hglm.gamma.gamma$ranef)
ranef.N <-data.frame("vu.N"=exp(hglm.poisson$ranef), "u.N"=hglm.poisson$ranef)
cbind(round(ranef.Y,2),round(ranef.N,2))


##Model matrix
X <-model.matrix(~GENDER+ENGINE, auto5)
Z <-fact2dummy(CAR_MAKE)

mu.Y <-exp(X%*%beta.Y[,1])
u.Y <-Z%*%ranef.Y[,2]
E.Y <-mu.Y*u.Y
severity <-unique(data.frame(E.Y[nz], "groups"=interaction(GENDER[nz],ENGINE[nz]), CAR_MAKE[nz]))
severity[1:10,]

mu.N <-exp(X%*%beta.N[,1])
u.N <-Z%*%ranef.N[,1]
E.N <-mu.N*u.N
frequency <-unique(data.frame(E.N, "groups"=interaction(GENDER,ENGINE), CAR_MAKE))
frequency[1:10,]
##Figure
bwtheme <- standard.theme("pdf", color=FALSE)
densityplot(~severity[,1]|severity[,3],xlab="HGLM claim severity", par.settings=bwtheme)

histogram(~frequency[,1]|frequency[,3], xlab="HGLM claim frequency", par.settings=bwtheme)


##Pure Premium for Clustered Portfolio
pi <-EXPOSURE*E.Y*E.N

histogram(~pi|GENDER[nz]+ENGINE[nz],  panel = function(...) {
        panel.histogram(...)
		panel.abline(v=mean(pi), lwd=3, lty=2)
       }, cex.axis=2, cex.lab=1.5,  xlab="HGLM pure premium", par.settings=bwtheme)
       
	   
#Example - Modeling the Aggregate Value of Claims with GLMM-TCPoisson
library(cplm)

glmm.tcpoisson <-cpglmm(CLAIM_AMOUNT~GENDER+ENGINE+(1|CAR_MAKE), 
offset=log(EXPOSURE), link="log", data=auto5)

glmm.tcpoisson@ghw
glmm.tcpoisson@p

data.frame("fixed.effects"=glmm.tcpoisson@fixef, "exp.fixed.effects"=exp(glmm.tcpoisson@fixef))
data.frame("random.effects"=glmm.tcpoisson@ranef,"exp.random.effects"=exp(glmm.tcpoisson@ranef))


pi1 <-glmm.tcpoisson@mu

histogram(~pi1|GENDER[nz]+ENGINE[nz],  panel = function(...) {
        panel.histogram(...)
		panel.abline(v=mean(pi1), lwd=3, lty=2)
       }, cex.axis=2, cex.lab=1.5,  xlab="HGLM TCPoisson pure premium", par.settings=bwtheme)
      
par(mfrow=c(1,2))
hist(pi, breaks=30, xlab="HGLM frequency-severity pure premium", main="")
abline(v=mean(pi), col="red", lty=2, lwd=3)
hist(pi1, breaks=30, xlab="HGLM TCPoisson pure premium", main="")
abline(v=mean(pi1), col="red", lty=2, lwd=3)

##A posteriori Ratemaking

#Example - NB Regression
library(MASS)
library(pscl)

auto1 <- read.csv2("http://web.ue.katowice.pl/woali/auto_data.csv", header=TRUE, sep=";", dec=",")
summary(auto1)

glm.nb<-glm.nb(CLAIM_COUNT~PREMIUM_SPLIT+SEX+CLIENT_AGE+CAR_AGE+POWER+CAPACITY, offset(EXPOSURE), link=log, data=auto1) 
summary(glm.nb)

glm.pois<-glm(CLAIM_COUNT~PREMIUM_SPLIT+SEX+CLIENT_AGE+CAR_AGE+POWER+CAPACITY, offset=EXPOSURE, family=poisson(log), data=auto1) 
vuong(glm.pois,glm.nb)

X <- 2 * (logLik(glm.nb) - logLik(glm.pois))
X
pchisq(X, df = 1, lower.tail=FALSE)

glm.nb2 <- update(glm.nb, . ~ . - CAPACITY)
anova(glm.nb, glm.nb2)

est <- cbind(Estimate = coef(glm.nb), confint(glm.nb))
round(est, digits = 2)
round(exp(est), digits = 2)


sample <- data.frame(
  EXPOSURE = c(1,1,0.5),
  PREMIUM_SPLIT = c(0,1,0),
  SEX=c(0,1,1),
  CLIENT_AGE=c('28-44','45-57','58-75'),
  CAR_AGE=c('0','1-16','17+'),
  POWER=c('67-124','66-','125+'),
  CAPACITY=c('900-','900-','901-2500')
)
sample $PHAT <- predict(glm.nb, sample , type = "response")
sample 


#Example - A Posteriori Premium - Cross-sectional Data
library(MixedPoisson)

auto3 <- read.csv2("http://web.ue.katowice.pl/woali/auto_data3.csv")

attach(auto3)
variable <-N
X <-model.matrix(~X1+X2)

lambda.new = c(); gamma.par.new = c()
lambda.new1 = c(); nu.new = c()
lambda.new2 = c(); delta.new = c()

lambda.matrix = matrix(nrow = length(variable), ncol = n.iter.max)
gamma.par.vec = vector(length = n.iter.max)
lambda.matrix1 = matrix(nrow = length(variable), ncol = n.iter.max)
nu.vec = vector(length = n.iter.max)
lambda.matrix2 = matrix(nrow = length(variable), ncol = n.iter.max)
delta.vec = vector(length = n.iter.max)

lambda.old = lambda_start(variable, X)$lambda
lambda.old1 = lambda.old2 = lambda.old
gamma.par.old = nu.old = delta.old =1

n.iter.max=4
for (i in 1:n.iter.max){
### E-Step
t = pseudo_values(variable, mixing=c("Gamma"), 
lambda=lambda.old, gamma.par=gamma.par.old, n=100)

t1 = pseudo_values(variable, mixing=c("lognorm"), 
lambda=lambda.old, nu=nu.old, n=100)

t2 = pseudo_values(variable, mixing=c("invGauss"), 
lambda=lambda.old, delta=delta.old, n=100)

### M-Step
lambda.new <-lambda_m_step(variable, X, t$pseudo_values)$lambda
lambda.new1 <-lambda_m_step(variable, X, t1$pseudo_values)$lambda
lambda.new2 <-lambda_m_step(variable, X, t2$pseudo_values)$lambda

gamma.par.new = est.gamma(t=t$pseudo_values)$gamma.par
nu.new = est.nu(t=t1$pseudo_values)$nu
delta.new = est.delta(t=t2$pseudo_values)$delta

gamma.par.old = gamma.par.new 
nu.old = nu.new 
delta.old = delta.new 

lambda.old = lambda.new
lambda.old1 = lambda.new1
lambda.old2 = lambda.new2

#Results
lambda.matrix[,i] <-lambda.new
gamma.par.vec[i] <-gamma.par.new
lambda.matrix1[,i] <-lambda.new1
nu.vec[i] <-nu.new
lambda.matrix2[,i] <-lambda.new2
delta.vec[i] <-delta.new
}

group<-interaction(X1, X2)
tab.gamma <-cbind(unique(data.frame(group)), unique(lambda.matrix))
colnames(tab.gamma) <-c("Group GENDER.AGE", "lambda i=1", "lambda i=2", "lambda i=3", "lambda i=4")
tab.gamma
gamma.par.vec

tab.lognorm <-cbind(unique(data.frame(group)), unique(lambda.matrix1))
colnames(tab.lognorm) <-c("Group GENDER.AGE", "lambda i=1", "lambda i=2", "lambda i=3", "lambda i=4")
tab.lognorm
nu.vec

tab.invGauss <-cbind(unique(data.frame(group)), unique(lambda.matrix2))
colnames(tab.invGauss) <-c("Group GENDER.AGE", "lambda i=1", "lambda i=2", "lambda i=3", "lambda i=4")
tab.invGauss
delta.vec

##Figure
par(mfrow=c(3,2))
boxplot(lambda.matrix[,4]~X1)
boxplot(lambda.matrix[,4]~X2)
boxplot(lambda.matrix1[,4]~X1)
boxplot(lambda.matrix1[,4]~X2)
boxplot(lambda.matrix2[,4]~X1)
boxplot(lambda.matrix2[,4]~X2)

#Example - Claim Frequency Mixed Model with Risk Profile for Longitudinal Data
library(hglm)
library(graphics)
library(lattice)
library(gridExtra)
library(StatMatch)

panel <-read.csv2("http://web.ue.katowice.pl/woali/panel.sp.csv", header=TRUE, sep=";", dec=",")
attach(panel)

##Figure
bwtheme <- standard.theme("pdf", color=FALSE)
histogram(CLAIM_COUNT~as.factor(BEGIN_YEAR)|GENDER:ENGINE:KIND_OF_PAYMENT, par.settings=bwtheme)

table(CLAIM_COUNT,BEGIN_YEAR)

table(CLAIM_COUNT,GENDER:ENGINE:KIND_OF_PAYMENT)

N.p.d=hglm(fixed=CLAIM_COUNT~ENGINE+GENDER+KIND_OF_PAYMENT, 
random=~1|CLIENT_NUMBER,rand.family=Gamma(log), disp=~GENDER,
offset=log(panel$EXPOSURE), family=poisson(log), data=panel)
summary(N.p.d)
print(N.p.d, print.ranef = TRUE)

N=hglm(fixed=CLAIM_COUNT~ENGINE+GENDER+KIND_OF_PAYMENT, 
random=~1|CLIENT_NUMBER,rand.family=Gamma(log), disp=~GENDER, rand.disp=~ENGINE,
offset=log(panel$EXPOSURE), family=poisson(log), data=panel)
summary(N.p.d)
print(N.p.d, print.ranef = TRUE)


##Model parameters
beta.N <-data.frame("Fixed effects"=N.p.d$fixef, "Tarrif rates"=exp(N.p.d$fixef))
beta.N
gamma.phi <- N.p.d$SummVC1
phi <-tapply(N.p.d$disp.fv, GENDER, unique)
phi
gamma.phi
phiu <-N.p.d$varRanef
phiu

v.N <-data.frame("v(u)"=exp(N.p.d$ranef), "Risk profile u"=N.p.d$ranef)
v.N[1:10,]

##Aposteriori premium
X <-model.matrix(~ENGINE+GENDER+KIND_OF_PAYMENT, panel)
Z <-fact2dummy(as.factor(CLIENT_NUMBER))

mu.N <-exp(X%*%beta.N[,1])
u.N <-Z%*%v.N[,2]
E.N <-mu.N*u.N
pi=EXPOSURE*E.N

summ <-data.frame(panel, "mu.N"=mu.N, "phi"=phi,"u.N"=u.N, "pi"=pi)

##Figure
par(mfrow=c(2,2))
hist(pi[BEGIN_YEAR=="2010"], xlab="A posteriori premium 2010", main="")
abline(v=mean(CLAIM_COUNT), lty=2, lwd=3)
hist(EXPOSURE*mu.N[BEGIN_YEAR=="2010"], xlab="A priori premium 2010", main="")
abline(v=mean(CLAIM_COUNT), lty=2, lwd=3)
plot(quantile(pi[BEGIN_YEAR=="2010"], probs=seq(0,1,0.005)), lwd=1, type="b", ylab="Quantiles 2010")
lines(quantile(EXPOSURE*mu.N[BEGIN_YEAR=="2010"], probs=seq(0,1,0.005)), lwd=3, type="b", xlab="")
legend("topleft", lwd=c(1,3), c("A posteriori premium 2010", "A priori premium 2010"))
boxplot(pi~as.factor(BEGIN_YEAR), xlab="", ylab="A posteriori premium 2007-2010")

##Figure
grid.arrange(
histogram(unique(summ$u.N), breaks=30, main="Random effects histogram", xlab="", par.settings=bwtheme),
bwplot(summ$u.N~GENDER:ENGINE:KIND_OF_PAYMENT, main="Random effects in groups", ylab="", par.settings=bwtheme))

hist(pi, xlab="Aposteriori premium", main="")
lines(density(pi))


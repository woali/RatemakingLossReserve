######### Authors: Alicja Wolny-Dominiak, Daniel Sobiecki 
######### University of Economics in Katowice, mailto: woali@ue.katowice.pl
######### R version 3.2.5

####CHAPTER 1

##List of packages 
library(stats4)
library(MASS) 
library(insuranceData)
library(tweedie)
library(statmod)
library(dglm)
library(faraway)
library(gamlss.dist)
library(pscl)
library(pscl)
library(insuranceData)
library(cplm)
library(pls) 
library(fitdistrplus)
library(gridExtra)
library(lattice)
library(ggplot2)
library(extRemes)
library(gPdtest)

options(digits=2)

#Example - the estimation of $(\lambda, \alpha)'$ in Gamma distribution
y.gam <- rgamma(100, 2)

#1st method:

l.gamma <- function(y,par) {n<-length(y)
-n*par[2]*log(par[1])+n*log(gamma(par[2]))-(par[2]-1)*sum(log(y))+par[1]*sum(y)} # -log-likelihood function			
result <- optim( par = c(1,2), l.gamma, y = y.gam, method = "Nelder-Mead")
result

#2nd method:

library(stats4) 
l.gamma<-function(lambda,alfa) {y<-y.gam
n<-length(y)
-n*alfa*log(lambda)+n*log(gamma(alfa))-(alfa-1)*sum(log(y))+lambda*sum(y)} # -log-likelihood function	
est<-mle(minuslog=l.gamma, start=list(lambda=1,alfa=2))
summary(est)

#3rd method:

library(MASS) 
fitdistr(y.gam,"gamma", method = "Nelder-Mead") # fitting gamma pdf parameters

##Example - fitting GLM Inverse-Gaussian
library(insuranceData)

data(dataCar)
attach(dataCar)

glm.ig <- glm(claimcst0[dataCar$claimcst0>0]~as.factor(veh_age[dataCar$claimcst0>0]), 
family=inverse.gaussian(link="log"), data=dataCar)
summary(glm.ig)

mean.group <-tapply(claimcst0[dataCar$claimcst0>0], as.factor(veh_age[dataCar$claimcst0>0]), mean)
data.frame("mean"=mean.group, "fitted values"=unique(glm.ig$fitted.values))

plot(glm.ig, which=1:4)

##Example - density of Tweedie distribution 
library(tweedie)
y=seq(0,4,length=1000)

d1 <-dtweedie( y, power=1.5, mu=1, phi=1)
d2 <-dtweedie( y, power=2, mu=1, phi=1)
d3 <-dtweedie( y, power=2.5, mu=1, phi=1)
d4 <-dtweedie( y, power=5, mu=1, phi=1)

## Figure
split.screen(c(2,2))
screen(1,new=TRUE)
plot(y, d1, type="l", lwd=2, col=1)
screen(2,new=TRUE)
plot(y, d2, type="l", lwd=2, col=2)
screen(3,new=TRUE)
plot(y, d3, type="l", lwd=2, col=3)
screen(4,new=TRUE)
plot(y, d4, type="l", lwd=2, col=4)


##Example - the estimation of $p$ 
library(tweedie)
y <- rnorm(100, 2)
p.fit <- tweedie.profile(y ~ 1, p.vec=seq(1.5, 3, by=0.1))

####Example -fitting GLM Tweeedie
library(tweedie)
library(statmod)
y <- rinvgauss(100,mean=1.5,dispersion=1)
x <- rpois(100,1)

p.fit <- tweedie.profile(y ~ 1, p.vec=seq(2, 3.5, by=0.1))
tweedie.reg <-glm(y~as.factor(x), 
family=tweedie(var.power=p.fit$p.max, link.power=1))
summary(tweedie.reg)

##Figure
par(mfrow = c(2,2)) 
plot(tweedie.reg, which = 1:4)

##Example - fitting DGLM Poisson 
library(dglm)
library(faraway)

attach(motorins)
Kilom <-as.factor(Kilometres)
Bonus <-as.factor(Bonus)
Make <-as.factor(Make)

dglm.claims <-dglm(Claims~Kilom+Bonus+Make, offset=log(motorins$Insured), 
~Make, family=poisson(link="log"), data=motorins)
summary(dglm.claims)

##Figure
par(mfrow=c(2,2))
plot(dglm.claims)

glm.claims <-glm(Claims~Kilom+Bonus+Make, offset=log(motorins$Insured),
family=poisson(link="log"), data=motorins)

plot(quantile(dglm.claims$fitted.values, probs=seq(0,1,0.00001)), 
quantile(glm.claims$fitted.values, probs=seq(0,1,0.00001)),
xlab="DGLM", 
ylab="GLM" )
abline(0,1,lwd=3,col="gray")



##Example - the estimation of $\lambda$ of Poisson distribution
k <- rpois(20,4)		
l.lambda <- function(lambda) {log(prod(((lambda^k)*exp(-lambda))/factorial(k)))}		
optimize(f=l.lambda, interval = c(0,10), maximum=TRUE)

##Example - simulation random sample from ZINB and ZIP
library(gamlss.dist)
library(pscl)

count.zinb <- rZINBI(1000, mu=5, sigma=1, nu=0.5) #random sample from ZINB distribution
count.zip <- rZIP(1000, mu=4, sigma=0.5) #random sample from ZIP distribution

## Figure
par(mfrow=c(1,2))
hist(count.zinb, main="")
hist(count.zip, main="")


##Example - estimation of parameters in ZI distribution
library(pscl)
library(insuranceData)
data(dataOhlsson)
number.of.claims <- dataOhlsson$antskad

table(number.of.claims)

zip.pois <-zeroinfl(number.of.claims~1, dist = "poisson")  
zip.nb <-zeroinfl(number.of.claims~1, dist = "negbin")  

lambda.pois <-exp(zip.pois$coefficients$count)
logit.tau.pois <-zip.pois$coefficients$zero
tau.pois <-exp(logit.tau.pois)/(1+exp(logit.tau.pois)) 

lambda.nb <-exp(zip.nb$coefficients$count)
logit.tau.nb <-zip.nb$coefficients$zero
tau.nb <-exp(logit.tau.nb)/(1+exp(logit.tau.nb)) 

##summary
cbind(lambda.pois, tau.pois. lambda.nb, tau.nb)

##Example - the probability mass function of ZTPois
dztpois <- function (y, lambda) 
{
    n <- length(y)
    dztp <- vector(length=n)
    dztp[y>0] <- dpois(y[y>0], lambda)/(1 - dpois(0, lambda))
    return(dztp)
}
k <- rpois(5,lambda=5)
round(dztpois(k[k>0], lambda=5), digits=4)

k.plot <- 1:10  
dk <- length(10) 
for (i in 1:10) {
dk[i] <- dztpois(i, lambda=5)
}
dk

##Figure
plot(k.plot, dk, type="h", xlab="k", ylab="dztpois", axes = FALSE)
axis(1, 1:10, 1:10)
axis(2, round(dk, digits=2))
par(new = T)
plot(k.plot, dk, type="p", xlab="", ylab="")
#box()

##Example - the estimation of $\lambda$ in ZTPois
k <- rpois(200,3)
k.nz <- as.data.frame(k[k>0])
X <- model.matrix(~1, k.nz)

beta <- coef(glm(k[k>0] ~ 1, family = poisson(link=log)))

ll.ztp <- function(beta){
lambda <- exp(X %*% beta)	
ll <-  sum(k[k>0]*log(lambda) - log(exp(lambda)-1)-log(factorial(k[k>0])))
return(-ll)
}
opt <- optim(beta, ll.ztp, hessian = TRUE, method = "BFGS")
log.likelihood <- (-opt$value)
opt$par
log.likelihood
opt$hessian

##Example - simulation Tweedie's compound Poisson distribution
library(tweedie)
library(statmod)

lambda=2
iota=30
varsigma=.2

mu=lambda*iota*varsigma
p=(iota+2)/(iota+1)
phi=lambda^(1-p)*(iota*varsigma)^(2-p)/(2-p)

n=1000
S=NULL
for (i in 1:n){
S[i] <-sum(rgamma(rpois(n,lambda),shape=iota,scale=varsigma))
S=cbind(S,S[i])
 }

##Figure
hist(S, breaks=30, main="")
abline(v=mean(S), lty=2, lwd=3)

##Example - estimation of Tweedie's compound Poisson parameters
library(insuranceData)
library(cplm)
data(IndustryAuto)

cpglm <- cpglm(Claim~1, link="log", data=IndustryAuto)
summary(cpglm) 
cbind("Intercept"=cpglm$coefficients, "Power p"=cpglm$p, "Dispersion"=cpglm$phi)

##Figure
hist(IndustryAuto$Claim, main="")


##Example - Cross-validation error
library(pls) 
dataset <- read.csv2("http://web.ue.katowice.pl/woali/auto_data2.csv", header=TRUE, sep=";", dec=",")

	mse.cv.glm=function(K){		
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
			Model=glm(CLAIM_AMOUNT~PREMIUM_SPLIT+POWER+CLIENT_AGE, weight=CLAIM_COUNT, family=Gamma(log), data=datasetTrainCV)
			pred.valid=NULL
			pred.valid=exp(predict(Model, datasetValidCV))
			MSE.valid=NULL  
			MSE.valid=sum((datasetValidCV$CLAIM_AMOUNT-pred.valid)^2)/length(pred.valid)
			ModelMSE=c(ModelMSE, MSE.valid)
			
			}
	MSE.CV=NULL
	MSE.CV=mean(ModelMSE)
	
outlist <- list(ModelMSE=ModelMSE, ModelRMSE=sqrt(ModelMSE), MSE.CV=MSE.CV, RMSE.CV=sqrt(MSE.CV))
return(outlist)
	}

mse.cv.glm(K=10)

#Example - Bootstrap MSEP and Confidence Interval
library(fitdistrplus)
library(tweedie)

dataset <- read.csv2("http://web.ue.katowice.pl/woali/auto_data2.csv", header=TRUE, sep=";", dec=",")
attach(dataset)

model<-glm(CLAIM_AMOUNT~PREMIUM_SPLIT+GENDER+POWER, weight=CLAIM_COUNT, family=Gamma(log), data=dataset)
res<-residuals(model, type = "pearson")
mu.hat<-model$fitted.values
phi.hat<-summary(model)$dispersion
p <-2
n <-length(res)

B <-1
difb <-dif <-0
for (i in 1:B) {
res.boot <- sample(abs(res), replace = TRUE)
yb <- res.boot*sqrt(phi.hat*mu.hat^p) + mu.hat
modelb <-glm(yb~PREMIUM_SPLIT+GENDER+POWER, weight=CLAIM_COUNT, family=Gamma(log), data=dataset)    
mub.hat<-modelb$fitted.values
phib.hat<-summary(modelb)$dispersion
ybb <-rtweedie(length(yb), mu = mub.hat, phi = phib.hat, power = p)
dif[i] <-sum((yb-ybb)^2)/n
difb <- difb+dif[i]
}
msep.boot <-difb/B
msep.boot
sqrt(msep.boot)


z.star <- double(B)
for (i in 1:B) {
    res.star <- rnorm(n, mean = mu.hat, sd = sigma.hat)
    z.star[i] <- (mean(res.star) - mu.hat) / sd(res.star)
}
crit.val <- quantile(z.star, probs = c(0.975, 0.025))
mu.hat - crit.val * sigma.hat
	
###Fixed Effect Models in Ratemaking

##Example - fitting GLM Gamma
library(gridExtra)
library(lattice)
bwtheme <- standard.theme("pdf", color=FALSE)

auto2 <- read.csv2("http://web.ue.katowice.pl/woali/auto_data2.csv", header=TRUE, sep=";", dec=",")
attach(auto2)

data.frame(table(PREMIUM_SPLIT:GENDER:POWER))

mean <-tapply(CLAIM_AMOUNT,PREMIUM_SPLIT:GENDER:POWER, mean)
data.frame(mean)

##Figure
grid.arrange(
densityplot(~CLAIM_AMOUNT|PREMIUM_SPLIT, xlab="", main="PREMIUM_SPLIT",par.settings=bwtheme),
densityplot(~CLAIM_AMOUNT|GENDER, xlab="", main="GENDER", par.settings=bwtheme),
densityplot(~CLAIM_AMOUNT|POWER, xlab="", main="POWER", par.settings=bwtheme))

glm.gamma <-glm(CLAIM_AMOUNT~PREMIUM_SPLIT+GENDER+POWER, weight=CLAIM_COUNT, family=Gamma(log), data=auto2) 
summary(glm.gamma)

##Figure
par(mfrow=c(2,2))
plot(glm.gamma)

relevel(x,ref,...)

one.way<-aggregate(auto2$CLAIM_AMOUNT/auto2$CLAIM_COUNT,
list(auto2$CLIENT_AGE), mean)
fit.gamma<-aggregate(glm.gamma$fitted.values,
list(glm.gamma$data$CLIENT_AGE), mean)
hist.cl_age<-aggregate(auto2$CLAIM_AMOUNT/auto2$CLAIM_COUNT,
list(auto2$CLIENT_AGE), na.omit(length))
labels<-levels(as.factor(auto2$CLIENT_AGE))

##Figure
par(mfrow=c(1,2))
 plot(one.way, type="o", lty=1, lwd=1, xlab="CLIENT_AGE",
ylab="Average claim" )
lines(fit.gamma)
points(one.way, pch=16, cex=2)
points(fit.gamma, pch=21, bg="white", cex=2)
legend(3,9500,c("observed","fitted"), pch=c(16,21))

barplot(hist.cl_age$x, beside=FALSE, names.arg=labels ,
xlab="CLIENT_AGE", ylab="Observations")

##Example -  Fitting GLM Poisson
library(lattice)
library(gridExtra)
bwtheme <- standard.theme("pdf", color=FALSE)

auto1 <- read.csv2("http://web.ue.katowice.pl/woali/auto_data.csv", header=TRUE, sep=";", dec=",")
attach(auto1)

table(CLAIM_COUNT)

X1 <- as.factor(auto1$PREMIUM_SPLIT)
X2 <- as.factor(auto1$SEX)
X3 <- as.factor(auto1$CLIENT_AGE)
X4 <- as.factor(auto1$CAR_AGE)
X5 <- as.factor(auto1$POWER)
X6 <- as.factor(auto1$CAPACITY)

sum <-tapply(auto1$CLAIM_COUNT,X1:X2:X3:X4:X5:X6, sum)

##Figure
grid.arrange(histogram(~auto1$CLAIM_COUNT|X1, panel = function(...) {
        panel.histogram(...)
		panel.abline(v=mean(CLAIM_COUNT), lwd=3, lty=2)
       },xlab="PREMIUM_SPLIT", par.settings=bwtheme),
histogram(~auto1$CLAIM_COUNT|X2, panel = function(...) {
        panel.histogram(...)
		panel.abline(v=mean(CLAIM_COUNT), lwd=3, lty=2)
       },xlab="GENDER", par.settings=bwtheme),
histogram(~auto1$CLAIM_COUNT|X3, panel = function(...) {
        panel.histogram(...)
		panel.abline(v=mean(CLAIM_COUNT), lwd=3, lty=2)
       },xlab="CLIENT_AGE", par.settings=bwtheme),
histogram(~auto1$CLAIM_COUNT|X4, panel = function(...) {
        panel.histogram(...)
		panel.abline(v=mean(CLAIM_COUNT), lwd=3, lty=2)
       },xlab="CAR_AGE", par.settings=bwtheme),
histogram(~auto1$CLAIM_COUNT|X5, panel = function(...) {
        panel.histogram(...)
		panel.abline(v=mean(CLAIM_COUNT), lwd=3, lty=2)
       },xlab="POWER", par.settings=bwtheme),
histogram(~auto1$CLAIM_COUNT|X6, panel = function(...) {
        panel.histogram(...)
		panel.abline(v=mean(CLAIM_COUNT), lwd=3, lty=2)
       },xlab="CAPACITY", par.settings=bwtheme), ncol=2, nrow=3)

glm.pois=glm(CLAIM_COUNT~PREMIUM_SPLIT+SEX+CLIENT_AGE+CAR_AGE+POWER+CAPACITY, offset=EXPOSURE, family=poisson(log), data=auto1) 
summary(glm.pois)

##Figure
par(mfrow=c(2,2))
plot(glm.pois)
anova(glm.pois, test="Chisq")

with(glm.pois, cbind(res.deviance = deviance, df = df.residual,
p = pchisq(deviance, df.residual, lower.tail=FALSE)))

glmexp.est <- exp(coef(glm.pois))
glmexp.est

glm.pois2=glm(CLAIM_COUNT~PREMIUM_SPLIT+SEX+SEX:CLIENT_AGE+CAR_AGE+POWER+CAPACITY, offset=EXPOSURE, family=poisson(log), data=auto1)
summary(glm.pois2) 
anova(glm.pois, glm.pois2, test="Rao") #Lagrange multiplier test (comparison of 2 GLMs)

AIC(glm.pois)
AIC(glm.pois2)
BIC(glm.pois)
BIC(glm.pois2)
pois$deviance/glm.pois$df.residual
pois2$deviance/glm.pois2$df.residual

##Example - ZIP Claim Frequency Model
library(pscl)
library(insuranceData)
data(dataCar)
attach(dataCar)

summary( zip1 <- zeroinfl(numclaims ~ gender + as.factor(agecat) + 
as.factor(veh_age)+offset(log(exposure)) | area, data=dataCar, dist=c("poisson")))
	
poi1 <- glm(numclaims ~ gender + as.factor(agecat) + as.factor(veh_age)+offset(log(exposure)), family=poisson, data=dataCar)
vuong(poi1, zip1)

plotdata <- expand.grid(gender=c("F","M"), agecat=factor(1:6), veh_age=factor(1:4), exposure=1, area=c("A"))
plotdata$freq_hat <- predict(zip1, plotdata, type=c("response"))
plotdata[1:10,]

#Figure
library(ggplot2)
ggplot(plotdata, aes(x=agecat, y=freq_hat, color=veh_age, group=veh_age))+geom_point()+geom_line()+
facet_wrap(~gender)+labs(x="Age category", y="Predicted number of claims")

##Example - fitting GLM TCPoisson
library(cplm)
library(tweedie)
library(statmod)
library(lattice)
library(gridExtra)
bwtheme <- standard.theme("pdf", color=FALSE)

auto4 <- read.csv2("http://web.ue.katowice.pl/woali/auto_data4.csv", header=TRUE, sep=";", dec=",")
attach(auto4)

grid.arrange(densityplot(~CLAIM_AMOUNT|PREMIUM_SPLIT, xlab="PREMIUM_SPLIT", par.settings=bwtheme),
densityplot(~CLAIM_AMOUNT|ENGINE, xlab="ENGINE", par.settings=bwtheme),
densityplot(~CLAIM_AMOUNT|CAPACITY, xlab="CAPACITY", par.settings=bwtheme),ncol=1, nrow=3)

p.fit <- tweedie.profile(CLAIM_AMOUNT ~ 1, p.vec=seq(1, 1.8, by=0.1), data=auto4)
p.fit$p.max
p.fit$phi.max
tweedie1 <-glm(CLAIM_AMOUNT~PREMIUM_SPLIT+ENGINE+CAPACITY,  
family=tweedie(var.power=p.fit$p.max, link.power=0), offset=log(EXPOSURE))
summary(tweedie1)$coefficients

tweedie2 <-cpglm(CLAIM_AMOUNT~PREMIUM_SPLIT+ENGINE+CAPACITY, offset=log(EXPOSURE), data=auto4)
rbind(data.frame(tweedie2@coefficients), "power p"=tweedie2@p, "dispersion"=tweedie2@phi)

mean <-tapply(CLAIM_AMOUNT,PREMIUM_SPLIT:ENGINE:CAPACITY, mean)
mu1 <-tapply(tweedie1$fitted.values,PREMIUM_SPLIT:ENGINE:CAPACITY, mean)
mu2 <-tapply(tweedie2$fitted.values,PREMIUM_SPLIT:ENGINE:CAPACITY, mean)

cbind(mean, mu1, mu2)

##Example - Zero-inflated TCPoisson Model
library(insuranceData)
data(dataCar)
attach(dataCar)

group <-interaction(gender,vah_age,agecat)
data.frame(levels(group))

tapply(claimcst0,gender, mean)
tapply(claimcst0,veh_age, mean)
tapply(claimcst0,agecat, mean)

mean(claimcst0)
var(claimcst0)

tweedie.zi <-zcpglm(claimcst0~gender+as.factor(veh_age)+agecat + offset(exposure)||agecat, data=dataCar)
tweedie.zi

##Example - Pure Risk Premium Estimation
auto4 <- read.csv2("http://web.ue.katowice.pl/woali/auto_data4.csv", header=TRUE, sep=";", dec=",")
attach(auto4)
group <-interaction(PREMIUM_SPLIT,ENGINE,CAPACITY)
data.frame(levels(group))
 
nz <-CLAIM_AMOUNT>0

glm.gamma <-glm(CLAIM_AMOUNT[nz]~PREMIUM_SPLIT[nz]+ENGINE[nz]+CAPACITY[nz], weight=CLAIM_COUNT[nz], family=Gamma(log), data=auto4) 
glm.poisson <-glm(CLAIM_COUNT~PREMIUM_SPLIT+ENGINE+CAPACITY, family=poisson(log), data=auto4) 

## summary
data.frame("severity"=round(glm.gamma$coefficients,2), "exp(sev)"=round(exp(glm.gamma$coefficients),2), 
"frequency"=round(glm.poisson$coefficients,2),"exp(freq)"=round(exp(glm.poisson$coefficients),2) )

severity <-tapply(glm.gamma$fitted.values, PREMIUM_SPLIT[nz]:ENGINE[nz]:CAPACITY[nz], unique)
frequency <-tapply(glm.poisson$fitted.values, PREMIUM_SPLIT:ENGINE:CAPACITY, unique)
pi <- severity*frequency
summary <-cbind(severity, frequency, data.frame(pi))

pi.portfel <-EXPOSURE[nz]*glm.gamma$fitted.values*glm.poisson$fitted.values[nz]

##Figure
bwtheme <- standard.theme("pdf", color=FALSE)
histogram(~pi.portfel|PREMIUM_SPLIT:ENGINE:CAPACITY, par.settings=bwtheme)

##Figure
par(mfrow=c(2,2))
hist(pi.portfel, breaks=30, main="Pure premium histogram")
boxplot(pi.portfel~PREMIUM_SPLIT[nz], main="PREMIUM_SPLIT")
boxplot(pi.portfel~ENGINE[nz], main="ENGINE")
boxplot(pi.portfel~POWER[nz], main="POWER")

#Matrix notation
X.Y <- model.matrix(~PREMIUM_SPLIT[nz]+ENGINE[nz]+CAPACITY[nz],auto4)
X.N <- model.matrix(~PREMIUM_SPLIT[nz]+ENGINE[nz]+CAPACITY[nz],auto4)
beta.Y <-coef(glm.gamma) 
beta.N <-coef(glm.poisson) 

pi.matrix <- exp(X.Y%*%beta.Y)*exp(X.N%*%beta.N)
pi.group <-  data.frame(Level[nz], pi.matrix)
unique(pi.group)
pi.portfel.matrix <- EXPOSURE[nz]*pi.matrix

##Example - Mean excess plot
library(extRemes)

dataset <- read.csv2("auto_data2.csv", header=TRUE, sep=";", dec=",")
x<-dataset[,"CLAIM_AMOUNT"]
mrlplot(x, umax=200000) 

##Example - Bootstrap goodness-of-fit test for GPD
library(gPdtest)

a<-double(53)
for (i in 0:53) {
	x2<-x[x>=50000-i*500]
	a[i]<-gpd.test(x2)
}
a

x3<-x[x>=25500]
gpd.fit(x3,"amle")


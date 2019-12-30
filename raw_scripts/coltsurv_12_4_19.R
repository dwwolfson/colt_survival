library(KMsurv)
library(survival)
library(gtools)
colt.dat<-read.csv("raw_data/fates_table.csv")

########glm conditional on survival to 2.75 weeks
#GLM 
library(Epi)
library(coda)
library(splines)
library(arm)
library(boot)
#Create Lexis diagram
Colt.Lexis<-read.csv("raw_data/fates_table.csv", header=TRUE)
spl.dat<-Lexis(entry=list(age=Cap.age.days),exit=list(age=age_at_fate_days),exit.status = censor, data=Colt.Lexis)
plot(spl.dat)
age.days<-seq(0,365,1)
dayobs<-splitLexis(spl.dat,breaks=age.days,time.scale="age")
dayobs$Fail<-factor(dayobs$lex.Xst)
write.csv(dayobs,file = "data_product/dayobs.csv",row.names=TRUE)
options(na.action = "na.fail")

m1.age <- glm(Fail ~ ns(age, df=1), family = binomial(link=cloglog), data=dayobs)
m2.age <- glm(Fail ~ ns(age, df=2), family = binomial(link=cloglog), data=dayobs)
m3.age <- glm(Fail ~ ns(age, df=3), family = binomial(link=cloglog), data=dayobs)
m4.age <- glm(Fail ~ ns(age, df=4), family = binomial(link=cloglog), data=dayobs)
m5.age <- glm(Fail ~ ns(age, df=5), family = binomial(link=cloglog), data=dayobs)
AIC(m1.age,m2.age,m3.age,m4.age,m5.age)

pred.m.age<-1 - predict(m1.age, type = "resp",data.frame(age = seq(round(7*2.7), 365, 1))) # start predicting at 2.7 weeks (or, min left truncation time)
St.m.age <- c(1, cumprod(pred.m.age))
#m.age <- glm(Fail ~ ns(age),family = binomial(link=cloglog), data=dayobs)
#pred.m.age<-1 - predict(m.age, type = "resp",data.frame(age = seq(0, 365, 1)))
#St.m.age <- c(1, cumprod(pred.m.age))

par(mar=c(5,6,4,1)+0.1)
plot(St.m.age, type="l",lwd=2,ylim=c(0,1),xlab="Age (days)",
     xlim=c(0,150), ylab=expression(hat(italic(S)))) #italics S-hat
################BOOTSTRAP 95% CIs################################################
nsims = 500
dayobs$Fail<-factor(dayobs$lex.Xst)
predictions <- matrix(NA, nrow = 348, ncol = nsims) 
uid <- unique(dayobs$colt_ID) 
nID <- length(uid) 
for(ii in 1:nsims){
  bootIDs <- data.frame(colt_ID = sample(x = uid, size = nID, replace = TRUE))
  bootDat <- merge(bootIDs, dayobs)
  dailyS <- glm(Fail ~ ns(age, df=3), family = binomial(link=cloglog), data=bootDat)
  dailySP <- 1 - predict(dailyS, type = "resp",data.frame(age = seq(round(7*2.7), 365, 1)))
  cum.boot <- c(1, cumprod(dailySP))
  predictions[,ii]<-cum.boot
}  
CIs<-apply(as.matrix(predictions), 1, function(x){quantile(x,probs=c(0.025,0.975))})
lines(CIs[1,], lty=3)
lines(CIs[2,], lty=3)

colt <- survfit(Surv(cond.age, censor) ~ 1, data=colt.dat)
lines(colt, col="red", xlim = c(0, 280))

##NOT WORKING##########
################BOOTSTRAP 95% CIs for conditional daily mort################################################
set.seed(5)
nsims = 10000
dayobs$Fail<-factor(dayobs$lex.Xst)
predictions.m <- matrix(NA, nrow = 348, ncol = nsims) #matrix for estimates
uid.m <- unique(dayobs$colt_ID)
nID.m <- length(uid.m)
for(iii in 1:nsims){
  bootIDs.m <- data.frame(colt_ID = sample(x = uid.m, size = nID.m, replace = TRUE))
  bootDat.m <- merge(bootIDs.m, dayobs)
  table(bootDat.m$colt_ID)
  length(dayobs$colt_ID) # original data
  length(bootDat.m$colt_ID) # bootstrap sample
  dailyS.m <- glm(Fail ~ ns(age, df=3),family = binomial(link=cloglog), data=bootDat.m)
  dailySP.m <- 1 - predict(dailyS.m, type = "resp",data.frame(age = seq(round(7*2.7), 365, 1)))
  cum.boot.m <- dailySP.m
  predictions.m[,iii]<-cum.boot.m
}  

m.CIs<-apply(as.matrix(predictions.m), 1, function(x){quantile(x,probs=c(0.025,0.975))})
lines(m.CIs[1,], lty=3)
lines(m.CIs[2,], lty=3)


############2-panel fig.  All models and CIs must be run beforehand
par(mfrow=c(1,2))
par(pin=c(3.75,3))
par(xaxs="i", yaxs="i") 


#par(mar=c(5,6,4,1)+0.1)
plot(St.m.age, frame.plot=F,type="l",lwd=2,ylim=c(0,1),xlab="Age (days)",
     ylab=expression(hat(italic(S)))) #italics S-hat
abline(h=0)
axis(1,at=1,labels=c("0"))
lines(CIs[1,], lty=5,lwd=2)
lines(CIs[2,], lty=5,lwd=2)
lines(kmfit, col="red", lwd=1.5, xlim = c(0, 280), xlab="Age (days)", ylab="Survival")
legend("topright", legend=c("A"),bty="n")
plot(1-(pred.m.age),type="l", 
     frame.plot=F,lwd=2,ylab="Daily mortality hazard", xlab="Age (days)")
abline(h=0.0003)
abline(v=1)
axis(1,at=1,labels=c("0"))
axis(2,at=0.000306,labels=c("0.000"))
lines(m.CIs[1,], lty=5, lwd=2)
lines(m.CIs[2,], lty=5,lwd=2)
legend("topright", legend=c("B"),bty="n")

#MuHAZ
library("muhaz")
kpfit1<-kphaz.fit(time=colt.dat$age_at_fate_days,status=colt.dat$censor, method="product-limit")
kpfit.sm1<-muhaz(time=colt.dat$age_at_fate_days, delta=colt.dat$censor, bw.method = "g")
# Use "kphaz.fit" to generate a hazard estimate
mar.default<-c(5,1,4,2) + 0.1
par(ps = 20, cex = 1, cex.main = 1,mar = mar.default + c(0, 4, 0, 0))
kphaz.plot(kpfit1)
lines(kpfit.sm1,col=2,lwd=3)
plot(kpfit.sm1)

head(kpfit1)
head(kpfit.sm1)
summary(kpfit.sm1)
#hazard plots for KM survival
fit1<-muhaz(colt.dat$age_at_fate_days,colt.dat$censor,bw.method="knn")
plot(fit1,col="black",lty=1, lwd=2)
h.df<-data.frame(est=fit1$est.grid, h.orig=fit1$haz.est)

for (i in 1:10000){
  d.s.fixedKM<-colt.dat[sample(1:nrow(colt.dat), nrow(colt.dat), replace = T),]
  d.s.muhaz<-muhaz(d.s.fixedKM$age_at_fate_days, d.s.fixedKM$censor, min.time=0, max.time=800, bw.method="knn")
  h.df<-cbind(h.df, d.s.muhaz$haz.est)
}


h.df$upper.ci<-apply(h.df[,c(-1,-2)], 1,  FUN=function(x) quantile(x, probs = 0.975))
h.df$lower.ci<-apply(h.df[,c(-1,-2)], 1,  FUN=function(x) quantile(x, probs = 0.025))
plot(h.df$est, h.df$h.orig, type="l", xlim=c(19,200),ylim=c(0,0.015), lwd=3, 
     xlab= "Age (days)", ylab = "Hazard")
lines(h.df$est, h.df$upper.ci,  lty=3, lwd=3)
lines(h.df$est, h.df$lower.ci,  lty=3, lwd=3)



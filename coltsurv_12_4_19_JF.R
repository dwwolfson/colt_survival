#' ---
#' title: "Crane Survival Analysis" 
#' output: html_document
#' ---
#' 
#' Load libraries
#+ warning=FALSE, message=FALSE
library(KMsurv)
library(survival)
library(gtools)
library(Epi)
library(coda)
library(splines)
library(arm)
library(boot)
library(dplyr)

#' Read in data
colt.dat<-read.csv("raw_data/fates_table.csv")

#' Check fates to make sure that data are entered correctly.
#' Note: there were 4 individuals that died or were censored on the same day 
#' they were captured
colt.dat %>% filter(age_at_fate_days <= Cap.age.days)

#' Change the age_at_fate_days for these individuals so that they die mid-way
#' through the day
colt.dat$age_at_fate_days<-with(colt.dat, ifelse(age_at_fate_days <= Cap.age.days, 
                                  age_at_fate_days + 0.5, age_at_fate_days))

#' ## GLM survival model 
#' 
#' ### First, create daily records using the splitLexis function
Colt.Lexis<-read.csv("raw_data/fates_table.csv", header=TRUE)
spl.dat<-Lexis(entry=list(age=Cap.age.days),exit=list(age=age_at_fate_days),exit.status = censor, data=Colt.Lexis)
plot(spl.dat)
age.days<-seq(0,365,1)
dayobs<-splitLexis(spl.dat,breaks=age.days,time.scale="age")
dayobs$Fail<-factor(dayobs$lex.Xst)
write.csv(dayobs,file = "data_product/dayobs.csv",row.names=TRUE)
options(na.action = "na.fail")

#' ### Fit models.
#' 
#' Note, using more htan 1 df results in warnings that 
#' the models result in fitted probabilities with 0 or 1, indicative
#' of problems with complete separation. So, only use age as linear effect.
m1.age <- glm(Fail ~ age, family = binomial(link=cloglog), data=dayobs)
 
#' Estimate the probability survive each day, given alive at start of the day.
#' Start predicting at 2.7 weeks (or, min left truncation time)
age<-seq(round(7*2.7), 365, 1) 
pred.m.age<-1 - predict(m1.age, type = "resp",data.frame(age=age)) 
St.m.age <- c(1, cumprod(pred.m.age))
age2<-c(18, age) # to allow plot of survival to age 19, given alive on day 18

 
#' ### Use a bootstrap to calculate 95% CIs 
#' 
#' Set up information to hold the bootstrap results
nsims = 500 # number of bootstraps
dayobs$Fail<-factor(dayobs$lex.Xst)

#' Store results in Shat (for survival curve) and hhat for hazard curve
Shat <- matrix(NA, nrow = 348, ncol = nsims) 
hhat <- matrix(NA, nrow = 347, ncol = nsims) 
uid <- unique(dayobs$colt_ID) 
nID <- length(uid) 
for(ii in 1:nsims){
  bootIDs <- data.frame(colt_ID = sample(x = uid, size = nID, replace = TRUE))
  bootDat <- merge(bootIDs, dayobs)
  dailyS <- glm(Fail ~ age, family = binomial(link=cloglog), data=bootDat)
  dailySP <- 1 - predict(dailyS, type = "resp",data.frame(age = seq(round(7*2.7), 365, 1)))
  cum.boot <- c(1, cumprod(dailySP))
  Shat[,ii]<-cum.boot
  hhat[,ii]<-dailySP
}  

#' Calculate confidence intervals for Shat (CIs) and hhat (m.CIs)
CIs<-apply(as.matrix(Shat), 1, function(x){quantile(x,probs=c(0.025,0.975))})
m.CIs<-apply(as.matrix(hhat), 1, function(x){quantile(x,probs=c(0.025,0.975))})

#' ## KM survival
#' 
#'  Make sure all records are complete for entry and exit dates & censor indicator
colt.dat %>% select("Cap.age.days", "age_at_fate_days", "censor") %>% complete.cases

#' Inspect the one obs with incomplete data
ind<-colt.dat %>% select("Cap.age.days", "age_at_fate_days", "censor") %>% 
      complete.cases 
colt.dat[ind==FALSE,]

#' Drop this obs when fitting model
colt <- survfit(Surv(Cap.age.days, age_at_fate_days, censor) ~ 1, data=colt.dat[ind,])
(kmfit<-summary(colt))

#' Create data set for plotting starting at 2.7 weeks = 19 days
#+ fig.keep="last" 
kmdat<-data.frame(time=c(19, kmfit$time), survival=c(1, kmfit$surv), 
                  lowCI=c(1,kmfit$lower), upCI=c(1, kmfit$upper))
  
#' Plot in a single 2-panel fig.   
#+ fig.height=6.5, fig.width=10, fig.keep="last" 
par(mfrow=c(1,2), bty="L")
par(mar=c(5,6,4,1)+0.1)
plot(age2, St.m.age, type="l",lwd=2,ylim=c(0,1),xlab="Age (days)",
     ylab=expression(paste(hat(italic(S)),"(t) | ",italic(S),"(19)=1")), xlim=c(18, 100)) #italics S-hat
lines(age2, CIs[1,], lty=5,lwd=2)
lines(age2, CIs[2,], lty=5,lwd=2)
lines(kmdat$time, kmdat$survival, col="red", type="s", lwd=1.5)
lines(kmdat$time, kmdat$lowCI, col="red", type="s", lty=2)
lines(kmdat$time, kmdat$upCI, col="red", type="s", lty=2)
legend("topright", legend=c("A"),bty="n")
plot(age, 1-(pred.m.age),type="l", 
     lwd=2,ylab="Daily mortality hazard", xlab="Age (days)",
      xlim=c(19, 100))
lines(age, 1-m.CIs[1,], lty=3, lwd=2)
lines(age, 1-m.CIs[2,], lty=3, lwd=2)
legend("topright", legend=c("B"),bty="n")
 
###########################################################################
###########################################################################
############### STATS 255 : FINAL EXAM ####################################
############### Lars Hertel ###############################################
###########################################################################
###########################################################################


###########################################################################
##### Loading Data
library(survival)
library(xtable)
library(Hmisc)
source("http://www.ics.uci.edu/~dgillen/STAT255/Handouts/Stat255Functions.R")
usrds <- read.csv("http://www.ics.uci.edu./~dgillen/STAT255/Data/usrdsData.csv")
albumin <- read.csv("http://www.ics.uci.edu./~dgillen/STAT255/Data/LongitudinalAlbumin.csv")
head(usrds)
summary(usrds)

###########################################################################
##### Simple data properties
mean(usrds$hist.cvd[usrds$age>=60 & usrds$age<=80], na.rm=T)
mean(usrds$hist.cvd[usrds$age>=20 & usrds$age<=40], na.rm=T)
mean(usrds$hist.cvd[usrds$age>=40 & usrds$age<=60], na.rm=T)
mean(usrds$hist.cvd[usrds$age<=80], na.rm=T)

###########################################################################
##### Missing Data
# NAs
# Racegroup, Smokegrp, Hist CVD, Diabetes, Under nourished, BMI, Cholest,
# Trigly, SBP
describe(usrds)
na.patterns <- naclus(usrds)
naplot( na.patterns )
par(mar=c(3, 2, 2, 0) + 0.1); plot( na.patterns ) ; par(mar=c(5, 4, 4, 2) + 0.1);
pdf("./plots/hist_na.pdf"); plot( summary( is.na(hist.cvd) ~ albumin.0 + esrdtime + diabetes + 
                 racegrp + female + age + smokegrp, data=usrds)); dev.off()
plot( summary( is.na(smokegrp) ~ albumin.0 + esrdtime + diabetes + 
                 racegrp + female + age + hist.cvd, data=usrds))
plot( summary( is.na(undnour) ~ albumin.0 + esrdtime + diabetes + 
                 racegrp + female + age, data=usrds))


###########################################################################
##### Univariate Summary Statistics
par(mfrow=c(1,1))
summary(usrds)

# Continuous
# albumin.0, cholest, trigly, pst.sbp, bmi, esrdtime, age
hist(usrds$albumin.0, col="darkred", xlab="Serum Albumin (g/dL) at recruitment", freq = F, main="Histogram of Serum Albumin values")
par(mfrow=c(1,5)); boxplot(usrds$age, main="Age", col="darkred");
boxplot(usrds$esrdtime, main="ESRD Time", ylab="Years", col="darkred");
boxplot(usrds$bmi, main="BMI", col="darkred");
boxplot(usrds$cholest, main="Cholest", col="darkred");
#boxplot(usrds$trigly, main="Trigly", col="darkred")
boxplot(usrds$pst.sbp, main="SBP", col="darkred"); par(mfrow=c(1,1))



par(mfrow=c(1,2))
id.died <- which(usrds$death==1)
hist(usrds$esrdtime[id.died]*365.25 + usrds$tdeath[id.died], col="darkred", xlab="ESRD time till death (days)", main="Distribution of survival times on ESRD\n for subjects that died", freq=T)
abline(v=median(usrds$esrdtime[id.died]*365.25 + usrds$tdeath[id.died]), col="black", lty=2)
id.cens <- which(usrds$death==0)
hist(usrds$esrdtime[id.cens]*365.25 + usrds$tdeath[id.cens], col="darkred", xlab="ESRD time till censoring (days)", main="Distribution of times on ESRD\n for subjects that did not die", freq=T)
abline(v=median(usrds$esrdtime[id.cens]*365.25 + usrds$tdeath[id.cens]), col="black", lty=2)
par(mfrow=c(1,1))


### Table
cont.names <- names(usrds)[c(2, 4, 10,12, 13, 14, 15,16)]
cont.vars <- usrds[,cont.names]
n <- nrow(cont.vars); p <- ncol(cont.vars)
cont.table <- matrix(0, p, 2)
for(i in 1:p){
  cont.table[i,2] <- sum(is.na(cont.vars[,i]))
  #cont.table[i, 1] <- cont.names[i]
  cont.table[i,1] <- paste(round(median(cont.vars[,i],na.rm=T),2), "(",round(sd(cont.vars[,i], na.rm=T),2),")", sep="")
}
cont.table <- as.matrix(cont.table)
rownames(cont.table) <- cont.names
colnames(cont.table) <- c("median (SD)", "Number of NAs")
xtable(cont.table)

###########################################################################
##### Pairs
pdf("./plots/pairs.pdf")
pairs(cont.vars, pch=".")
dev.off()
  


# Categorical
cat.names <- names(usrds)[c(3, 5, 6, 7, 8, 9, 11)]
cat.vars <- usrds[,cat.names]
n <- nrow(cat.vars); p <- ncol(cat.vars)
cat.table <- matrix(0, p+6, 2)
levelcount <- 0
i=1
for(i in 1:p){
  var <- cat.vars[!is.na(cat.vars[,i]),i]
  no_na <- sum(is.na(cat.vars[,i]))
  if(max(var)>1){
    cat.table[i+levelcount,1] <- "---"
    cat.table[i+levelcount,2] <- no_na
    for(j in 1:length(unique(var))){
      cat.table[i+levelcount+j,1] <- paste(round(sum(var==unique(var)[j]),2), "(",round(mean(var==unique(var)[j])*100,2),"%)", sep="")
      cat.table[i+levelcount+j,2] <- "---"
    }
    levelcount <- levelcount + j
  }else{
    cat.table[i+levelcount,1] <- paste(round(sum(var),2), "(",round(mean(var),2)*100,"%)", sep="")
    cat.table[i+levelcount,2] <- no_na
  }
}
cat.table <- as.matrix(cat.table)
rownames(cat.table) <- c(cat.names[1], cat.names[2], cat.names[3], "caucasian", "african", "other", cat.names[4],
"never smoked", "former smoker", "current smoker", cat.names[5], cat.names[6], cat.names[7])
colnames(cat.table) <- c("count (proportion)", "Number NAs")
xtable(cat.table)



# Time variable
fit <- survfit(Surv(tdeath, death)~1, data=usrds)
plot(fit)

sum(usrds$death); mean(usrds$death)
par(mfrow=c(1,2)); hist(usrds$tdeath, freq=F,col="darkred", main = "Distribution of time in study", xlab="Time in study (in days)") # Most people censored ?
hist(usrds$tdeath[usrds$death==1], freq=F,col="darkred", main = "Distribution of time to death", xlab="Time to death (in days)"); par(mfrow=c(1,1));  # Most people who die, die at the beginning
mean(usrds$tdeath==427)
sum(usrds$tdeath==427 & usrds$death==1) # Only one person died at the end of study
# So ~ 75% of people get censored 

fit <- survfit(Surv(tdeath, death)~cut(albumin.0, quantile(albumin.0, c(0.,.33,.66,.1))), data=usrds)
par(mfrow=c(1,1)) ;plot(fit, lty=1:3, col="darkred", lwd=2, main="Kaplan-Meier Plot of Death by albumin level", ylab="Survival", xlab="Time in study (days)"); legend( 0.25, .5, lty=1:3, lwd=2, col="darkred", legend=c("low albumin","med albumin","high albumin"), bty="n" )
kmPlot(fit, groupLabels = c("low albumin","med albumin","high albumin")); par(mar=c(5, 4, 4, 2) + 0.1); legend( 0.25, .5, lty=1:3, legend=c("low albumin","med albumin","high albumin"), bty="n" )
pKM(fit, q=c(0.,0.5,1.))



# Bivariate statistics
# Association of serum albumin with predictors
par(mfrow=c(3,3))
par(mar=c(3, 2, 2, 0) + 0.1); 
boxplot(albumin.0~factor(death, labels=c("no death", "death")), usrds, col="darkred", ylab="Serum Albumin");
boxplot(albumin.0~factor(female, labels=c("male", "female")), usrds, col="darkred", ylab="Serum Albumin");
boxplot(albumin.0~factor(racegrp, labels=c("cauc", "African Am", "other")), main="Race", usrds, col="darkred", ylab="Serum Albumin");
boxplot(albumin.0~factor(smokegrp, labels=c("never", "former", "current")), main="Smoking" ,usrds, col="darkred", ylab="Serum Albumin");
boxplot(albumin.0~factor(hist.cvd, labels=c("No", "Yes")), main="CVD History", usrds, col="darkred", ylab="Serum Albumin");
boxplot(albumin.0~factor(diabetes, labels=c("No", "Yes")), main="Diabetes",usrds, col="darkred", ylab="Serum Albumin");
boxplot(albumin.0~factor(undnour, labels=c("No", "Yes")), main="Undernourished", usrds, col="darkred", ylab="Serum Albumin"); 
boxplot(usrds$albumin.0~cut(usrds$bmi, quantile(usrds$bmi, c(0.,.5,.1), na.rm=T)), main="BMI", col="darkred", ylab="Serum Albumin");
boxplot(usrds$albumin.0~cut(usrds$esrdtime, quantile(usrds$esrdtime, c(0.,.5,.1), na.rm=T)), main="ESRD Time (years)", col="darkred", ylab="Serum Albumin");
par(mfrow=c(1,1)); par(mar=c(5, 4, 4, 2) + 0.1); 

boxplot(usrds$albumin.0~as.character(as.numeric(usrds$undnour)), main="Undernourished",col="darkred", ylab="Serum Albumin"); 

boxplot(esrdtime~factor(undnour, labels=c("No", "Yes")), main="Undernourished", usrds, col="darkred", ylab="ESRD time"); 

par(mfrow=c(1,2)); plot(usrds$bmi, usrds$albumin.0); plot(usrds$esrdtime, usrds$albumin.0)

boxplot(usrds$albumin.0~cut(usrds$bmi, quantile(usrds$bmi, c(0.,.5,.1), na.rm=T)), main="BMI", col="darkred", ylab="Serum Albumin");
boxplot(usrds$albumin.0~cut(usrds$esrdtime, quantile(usrds$esrdtime, c(0.,.5,.1), na.rm=T)), main="ESRD Time (years)", col="darkred", ylab="Serum Albumin");




###########################################################################
##### Follow up measurements of albumin
head(albumin)
plot(albumin$measday, albumin$albumin, ylab="Serum Albumin (g/dL)", xlab="Days since recruitment", main="Development of Serum Albumin")

slope.1 <- NULL
count <- 1
for(i in 1:length(unique(albumin$usrds.id))){
  #i=1
  id <- unique(albumin$usrds.id)[i]
  idxs <- which(albumin$usrds.id==id)
  if(length(idxs)>1){
    slope.1[count] <- summary(lm(albumin~measday, data=albumin, subset=idxs))$coefficients[2,1]
    count <- count + 1
  }
}
slope.2 <- NULL
count <- 1
for(i in 1:length(unique(albumin$usrds.id))){
  #i=1
  id <- unique(albumin$usrds.id)[i]
  idxs <- which(albumin$usrds.id==id)
  if(length(idxs)>1 & usrds$death[which(usrds$usrds.id==id)]==1){
    slope.2[count] <- summary(lm(albumin~measday, data=albumin, subset=idxs))$coefficients[2,1]
    count <- count + 1
  }
}

par(mfrow=c(1,2)); hist(slope.1, col="darkred", freq=F, main="Distribution of slopes in LM fits to\n albumin values", xlab="Slope coefficient")
hist(slope.2, col="darkred", freq=F, main="Distribution of slopes in LM fits to\n albumin values for subjects that died", xlab="Slope coefficient"); par(mfrow=c(1,1)); 



###########################################################################
##### Ties
length(which(table(usrds$tdeath)!=1))
# We have a lot of ties
# Model uses the Efron approximation

###########################################################################
##### Print results
results.tab <- function(fit){
  results.fit <- cbind(paste(round(summary(fit)$conf.int[,1],2),"(", round(summary(fit)$conf.int[,3],2)," ,", round(summary(fit)$conf.int[,4],2),")", sep=""), round(summary(fit)$coefficients[,5],3))
  colnames(results.fit) <- c("Adjusted relative risk (95% CI)", "p-Value")
  return(xtable(results.fit))
}

###########################################################################
##### Fit model
fit <- coxph(Surv(tdeath, death)~albumin.0 + esrdtime + bmi + factor(smokegrp)
      + cholest + diabetes + pst.sbp + female + age + hist.cvd, data = usrds)
summary(fit)
results.tab(fit)



fit.grp <- coxph(Surv(tdeath, death)~cut(albumin.0, quantile(albumin.0, c(0.,.33,.66,.1))) + esrdtime + bmi + factor(smokegrp)
             + cholest + diabetes + pst.sbp + female + age + hist.cvd, data = usrds)
summary(fit.grp)
# Check for proportional hazards
par(mfrow=c(1,1)) ;plot(survfit(fit.grp), lty=1:3, fun="cloglog"); legend( 0.25, .5, lty=1:3, legend=c("low","med","high"), bty="n" )
# Doesn't look proportional

###########################################################################
##### Interaction model
fit.inter <- coxph(Surv(tdeath, death)~scale(albumin.0, center=T, scale=F)*hist.cvd + esrdtime + bmi + factor(smokegrp)
             + cholest + diabetes + pst.sbp + female + age, data = usrds)
summary(fit.inter)
results.tab(fit.inter)
exp(-.98*3.5 + .80 + .45 * 3.5)/exp(-.98*3.5)
exp(-.98*4.0 + .80 + .45 * 4.0)/exp(-.98*4.0)

###########################################################################
##### Imputating under nourishment with GLM
undnour.glm <- glm(undnour~.-tdeath-death-usrds.id, data=usrds, family="binomial")
summary(undnour.glm)
exp(summary(undnour.glm)$coefficients[,1])

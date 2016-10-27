###########################################################################
###########################################################################
############### STATS 255 : FINAL EXAM ####################################
############### Lars Hertel ###############################################
###########################################################################
###########################################################################


###########################################################################
##### Diagnostics Model A
fit <- coxph(Surv(tdeath, death)~albumin.0 + esrdtime + bmi + factor(smokegrp)
             + cholest + diabetes + pst.sbp + female + age + hist.cvd, data = usrds)
summary(fit)



###########################################################################
##### Functional Form: Martingale Residuals
pdf("./plots/modela_diag.pdf")
par(mfrow=c(2,2))
#### albumin linearly
fit.alb <- coxph(Surv(tdeath, death)~esrdtime + bmi + factor(smokegrp)
             + cholest + diabetes + pst.sbp + female + age + hist.cvd, data = usrds[-c(388),])

# Functional form of Serum Albumin
mresids <- residuals(fit.alb, type="martingale")
lmfit.alb <- lm(albumin.0 ~ esrdtime + bmi + factor(smokegrp)
                + cholest + diabetes + pst.sbp + female + age + hist.cvd, data = usrds[-c(388),])
ralb <- lmfit.alb$resid
ord <- order(ralb)
mresids <- mresids[ord]
ralb <- ralb[ord]
plot(ralb, mresids, xlab="LM Residual Albumin",
     ylab="Martingale Residual", main="Model A: Functional Form of Albumin")
lines(smooth.spline(ralb, mresids, df=6), col="darkred", lwd=2)
lines(ralb, fitted(lm(mresids~ralb)), col="darkblue", lwd=2)


#### albumin linearly
fit.alb <- coxph(Surv(tdeath, death)~esrdtime + bmi + factor(smokegrp)
                 + cholest + diabetes + pst.sbp + female + age + hist.cvd, data = usrds[-c(388),])

# Functional form of Serum Albumin
mresids <- residuals(fit.alb, type="martingale")
lmfit.alb <- lm(log(albumin.0) ~ esrdtime + bmi + factor(smokegrp)
                + cholest + diabetes + pst.sbp + female + age + hist.cvd, data = usrds[-c(388),])
ralb <- lmfit.alb$resid
ord <- order(ralb)
mresids <- mresids[ord]
ralb <- ralb[ord]
plot(ralb, mresids, xlab="LM Residual log(Albumin)",
     ylab="Martingale Residual", main="Model A: Functional Form of Albumin")
lines(smooth.spline(ralb, mresids, df=6), col="darkred", lwd=2)
lines(ralb, fitted(lm(mresids~ralb)), col="darkblue", lwd=2)
#identify(ralb, mresids, labels=usrds$usrds.id[ord], plot=T)

### REMOVE ROW 388
usrds$usrds.id[388]


###########################################################################
##### Outliers: Deviance Residuals

dresids <- residuals(fit, type="deviance")
lp <- predict(fit, type="lp")
plot(lp, dresids, xlab="Linear Predictor", ylab="Deviance Residual", main="Model A: Outliers")
abline(h=c(-3, 3))
# The residuals look weird


###########################################################################
##### Proportional Hazards Assumption

fit.grp <- coxph(Surv(tdeath, death)~strata(cut(albumin.0, quantile(albumin.0, c(0.,.33,.66,.1)))) + esrdtime + bmi + factor(smokegrp)
                 + cholest + diabetes + pst.sbp + female + age + hist.cvd, data = usrds)
summary(fit.grp)
# Check for proportional hazards
plot(survfit(fit.grp), lty=1:3, fun="cloglog", main="Model A: Proportional Hazards", xlab="Time in study (days)", ylab="cumulative log Hazard"); legend( 1, -1, lty=1:3, legend=c("low Albumin","med Albumin","high Albumin"), bty="n" )
# Looks "piecewise" parallel
dev.off()


###########################################################################
##### Influence

dfbeta <- residuals(fit, type="dfbeta")
head(dfbeta)
colnames(dfbeta) <- names(fit$coef)
summary(dfbeta)
ids <- complete.cases(usrds[,c("albumin.0", "esrdtime", "bmi", "smokegrp",
"cholest", "diabetes", "pst.sbp", "female", "age", "hist.cvd")])
# No single very influential observation
pdf("./plots/modela_dfbeta.pdf")
par(mfrow=c(4,3))
par(mar=c(3, 2, 2, 0) + 0.1);
for( i in 1:ncol(dfbeta) ){
  dfbeta <- residuals(fit, type="dfbeta")
  colnames(dfbeta) <- names(fit$coef)
  plot( usrds$usrds.id[ids], dfbeta[,i], xlab="Patient ID",
          ylab="Delta Beta", main=dimnames(dfbeta)[[2]][i] )
  large.idxs <- match(sort(-abs(dfbeta[,i]))[1:3], -abs(dfbeta[,i]))
  text( usrds$usrds.id[ids][large.idxs]-2, dfbeta[large.idxs,i], usrds$usrds.id[ids][large.idxs], col="darkred" )
}
dev.off()
par(mar=c(5, 4, 4, 2) + 0.1);
par(mfrow=c(1,1))


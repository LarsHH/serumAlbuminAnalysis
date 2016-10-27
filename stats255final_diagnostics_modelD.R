###########################################################################
###########################################################################
############### STATS 255 : FINAL EXAM ####################################
############### Lars Hertel ###############################################
###########################################################################
###########################################################################


###########################################################################
##### Diagnostics Model D
fit.tv <- coxph(Surv(start, stop, death)~albumin + esrdtime + bmi + factor(smokegrp)
                + cholest + diabetes + pst.sbp + female + age + hist.cvd, data = usrds.texpand)
summary(fit.tv)

###########################################################################
##### Functional Form: Martingale Residuals

#### albumin linearly
fit.tv.alb <- coxph(Surv(start, stop, death)~esrdtime + bmi + factor(smokegrp)
                    + cholest + diabetes + pst.sbp + female + age + hist.cvd, data = usrds.texpand)
pdf("./plots/modelc_diag.pdf")
par(mfrow=c(2,2))
# Functional form of Serum Albumin
mresids <- residuals(fit.tv.alb, type="martingale")
lmfit.alb <- lm(albumin ~ esrdtime + bmi + factor(smokegrp)
                + cholest + diabetes + pst.sbp + female + age + hist.cvd, data = usrds.texpand)
ralb <- lmfit.alb$resid
ord <- order(ralb)
mresids <- mresids[ord]
ralb <- ralb[ord]
plot(ralb, mresids, xlab="LM Residual Albumin",
     ylab="Martingale Residual", main="Model C: Functional Form of Albumin")
lines(smooth.spline(ralb, mresids, df=6), col="darkred", lwd=2)
lines(ralb, fitted(lm(mresids~ralb)), col="darkblue", lwd=2)


#### albumin linearly
fit.tv.alb <- coxph(Surv(start, stop, death)~esrdtime + bmi + factor(smokegrp)
                    + cholest + diabetes + pst.sbp + female + age + hist.cvd, data = usrds.texpand)

# Functional form of Serum Albumin
mresids <- residuals(fit.tv.alb, type="martingale")
lmfit.alb <- lm(log(albumin) ~ esrdtime + bmi + factor(smokegrp)
                + cholest + diabetes + pst.sbp + female + age + hist.cvd, data = usrds.texpand)
ralb <- lmfit.alb$resid
ord <- order(ralb)
mresids <- mresids[ord]
ralb <- ralb[ord]
plot(ralb, mresids, xlab="LM Residual log(Albumin)",
     ylab="Martingale Residual", main="Model C: Functional Form of Albumin")
lines(smooth.spline(ralb, mresids, df=6), col="darkred", lwd=2)
lines(ralb, fitted(lm(mresids~ralb)), col="darkblue", lwd=2)


###########################################################################
##### Outliers: Deviance Residuals

dresids <- residuals(fit.tv, type="deviance")
lp <- predict(fit.tv, type="lp")
plot(lp, dresids, xlab="Linear Predictor", ylab="Deviance Residual", main="Model C: Outliers")
abline(h=c(-3, 3))
# The residuals look weird


###########################################################################
##### Proportional Hazards Assumption

fit.grp <- coxph(Surv(start, stop, death)~strata(cut(albumin, quantile(albumin, c(0.,.33,.66,.1)))) + esrdtime + bmi + factor(smokegrp)
                 + cholest + diabetes + pst.sbp + female + age + hist.cvd, data = usrds.texpand)
summary(fit.grp)
# Check for proportional hazards
plot(survfit(fit.grp), lty=1:3, fun="cloglog", main="Model C: Proportional Hazards", xlab="Time in study (days)", ylab="cumulative log Hazard"); legend( 1, -1, lty=1:3, legend=c("low Albumin","med Albumin","high Albumin"), bty="n" )
# ?????
# Does this still make sense now???
dev.off()
par(mfrow=c(1,1))
###########################################################################
##### Influence

dfbeta <- residuals(fit.tv, type="dfbeta")
head(dfbeta)
colnames(dfbeta) <- names(fit.tv$coef)
summary(dfbeta)
ids <- complete.cases(usrds.texpand[,c("albumin.0", "esrdtime", "bmi", "smokegrp",
                                       "cholest", "diabetes", "pst.sbp", "female", "age", "hist.cvd")])
# No single very influential observation
pdf("./plots/modelc_dfbeta.pdf")
par(mfrow=c(4,3))
par(mar=c(3, 2, 2, 0) + 0.1);
for( i in 1:ncol(dfbeta) ){
  dfbeta <- residuals(fit.tv, type="dfbeta")
  colnames(dfbeta) <- names(fit.tv$coef)
  plot( usrds.texpand$usrds.id[ids], dfbeta[,i], xlab="Patient ID",
        ylab="Delta Beta", main=dimnames(dfbeta)[[2]][i] )
  large.idxs <- match(sort(-abs(dfbeta[,i]))[1:3], -abs(dfbeta[,i]))
  text( usrds.texpand$usrds.id[ids][large.idxs]-2, dfbeta[large.idxs,i], usrds.texpand$usrds.id[ids][large.idxs], col="darkred" )
}
dev.off()
par(mar=c(5, 4, 4, 2) + 0.1);
par(mfrow=c(1,1))



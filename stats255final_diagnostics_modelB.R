###########################################################################
###########################################################################
############### STATS 255 : FINAL EXAM ####################################
############### Lars Hertel ###############################################
###########################################################################
###########################################################################


###########################################################################
##### Diagnostics Model B
fit.inter <- coxph(Surv(tdeath, death)~albumin.0*hist.cvd + esrdtime + bmi + factor(smokegrp)
                   + cholest + diabetes + pst.sbp + female + age, data = usrds)
summary(fit.inter)

###########################################################################
##### Functional Form: Martingale Residuals
# Same as without interaction

###########################################################################
##### Outliers: Deviance Residuals
pdf("./plots/modelb_diag.pdf")
dresids <- residuals(fit.inter, type="deviance")
lp <- predict(fit.inter, type="lp")
plot(lp, dresids, xlab="Linear Predictor", ylab="Deviance Residual")
abline(h=c(-2.5, 2.5))
dev.off()
# Same as for model A


###########################################################################
##### Influence

dfbeta <- residuals(fit.inter, type="dfbeta")
colnames(dfbeta) <- names(fit.inter$coef)
summary(dfbeta)
summary(fit.inter)
# We have some pretty influential observations on hist.cvd now

dfbeta <- residuals(fit.inter, type="dfbeta")
head(dfbeta)
colnames(dfbeta) <- names(fit.inter$coef)
summary(dfbeta)
ids <- complete.cases(usrds[,c("albumin.0", "esrdtime", "bmi", "smokegrp",
                               "cholest", "diabetes", "pst.sbp", "female", "age", "hist.cvd")])
# No single very influential observation
pdf("./plots/modelb_dfbeta.pdf")
par(mfrow=c(4,3))
par(mar=c(3, 2, 2, 0) + 0.1);
for( i in 1:ncol(dfbeta) ){
  dfbeta <- residuals(fit.inter, type="dfbeta")
  colnames(dfbeta) <- names(fit.inter$coef)
  plot( usrds$usrds.id[ids], dfbeta[,i], xlab="Patient ID",
        ylab="Delta Beta", main=dimnames(dfbeta)[[2]][i] )
  large.idxs <- match(sort(-abs(dfbeta[,i]))[1:3], -abs(dfbeta[,i]))
  text( usrds$usrds.id[ids][large.idxs]-2, dfbeta[large.idxs,i], usrds$usrds.id[ids][large.idxs], col="darkred" )
}
dev.off()
par(mar=c(5, 4, 4, 2) + 0.1);
par(mfrow=c(1,1))

###########################################################################
##### Proportional Hazards Assumption
rm.idxs <- which(!is.na(match(usrds$usrds.id, c(583571, 596385, 648199))))
fit.inter <- coxph(Surv(tdeath, death)~scale(albumin.0, center=T, scale=F)*hist.cvd + esrdtime + bmi + factor(smokegrp)
                   + cholest + diabetes + pst.sbp + female + age, data = usrds[-rm.idxs,])
summary(fit.inter)

# Cannot stratify by interaction
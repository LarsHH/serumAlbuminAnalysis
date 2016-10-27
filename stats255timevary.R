###########################################################################
###########################################################################
############### STATS 255 : FINAL EXAM ####################################
############### Lars Hertel ###############################################
###########################################################################
###########################################################################

###########################################################################
##### Time varying albumin
usrds <- read.csv("http://www.ics.uci.edu./~dgillen/STAT255/Data/usrdsData.csv")
albumin <- read.csv("http://www.ics.uci.edu./~dgillen/STAT255/Data/LongitudinalAlbumin.csv")

usrds.texpand <- usrds[,c("tdeath", "death", "undnour", "esrdtime", "bmi",
"smokegrp", "cholest", "diabetes", "pst.sbp", "female", "age", "hist.cvd")]

unique(albumin$usrds.id)[1:10]
idxs.one <- which(as.numeric(table(albumin$usrds.id))==1)
ids.one <- usrds$usrds.id[idxs.one]
usrds.texpand[which(usrds.texpand$usrds.id==ids.one[2]),]

usrds.texpand <- merge(usrds, albumin, by="usrds.id")
usrds.texpand$start <- rep(NA, length(usrds.texpand$usrds.id))
usrds.texpand$stop <- rep(NA, length(usrds.texpand$usrds.id))

# Creating start / stop times
for(i in 1:length(unique(albumin$usrds.id))){
  id <- unique(albumin$usrds.id)[i]
  no_entries <- sum(usrds.texpand$usrds.id==id)
  idxs <- which(usrds.texpand$usrds.id==id)
  usrds.texpand$start[idxs] <- usrds.texpand$measday[idxs]
  usrds.texpand$stop[idxs] <- c(usrds.texpand$measday[idxs][2:no_entries], usrds.texpand$tdeath[idxs][1])
  #if(id%in%ids.one){
  #  usrds.texpand$stop[idxs] <- usrds.texpand$start[idxs]
  #}
}
usrds.texpand[1:100,c("usrds.id","measday", "start", "stop", "albumin")]
usrds.texpand[is.na(usrds.texpand$stop),"stop"] <- usrds.texpand[is.na(usrds.texpand$stop),"tdeath"]

###########################################################################
##### Fit model
fit.tv <- coxph(Surv(start, stop, death)~albumin + esrdtime + bmi + factor(smokegrp)
             + cholest + diabetes + pst.sbp + female + age + hist.cvd, data = usrds.texpand)
summary(fit.tv)
results.tab(fit.tv)

fit.tv.inter <- coxph(Surv(start, stop, death)~scale(albumin, center=T, scale=F)*hist.cvd + esrdtime + bmi + factor(smokegrp)
                + cholest + diabetes + pst.sbp + female + age, data = usrds.texpand)
summary(fit.tv.inter)
results.tab(fit.tv.inter)

##### Under nourished
fit.tv <- coxph(Surv(start, stop, death)~albumin + undnour+ esrdtime + bmi + factor(smokegrp)
                + cholest + diabetes + pst.sbp + female + age + hist.cvd, data = usrds.texpand)
summary(fit.tv)
results.tab(fit.tv)

# SQRT ESRD time
fit.tv.sqrt <- coxph(Surv(start, stop, death)~albumin + sqrt(esrdtime) + bmi + factor(smokegrp)
                + cholest + diabetes + pst.sbp + female + age + hist.cvd, data = usrds.texpand)
summary(fit.tv.sqrt)
results.tab(fit.tv.sqrt)

fit.tv.log <- coxph(Surv(start, stop, death)~log2(albumin) + esrdtime + bmi + factor(smokegrp)
                + cholest + diabetes + pst.sbp + female + age + hist.cvd, data = usrds.texpand)
summary(fit.tv.log)
results.tab(fit.tv.log)

### Removing BMI
fit.tv <- coxph(Surv(start, stop, death)~albumin + esrdtime + factor(smokegrp)
                + cholest + diabetes + pst.sbp + female + age + hist.cvd, data = usrds.texpand)
summary(fit.tv)
results.tab(fit.tv)

### Removing cholesterol
fit.tv <- coxph(Surv(start, stop, death)~albumin + esrdtime + bmi + factor(smokegrp)
                + diabetes + pst.sbp + female + age + hist.cvd, data = usrds.texpand)
summary(fit.tv)
results.tab(fit.tv)

### 
fit.tv <- coxph(Surv(start, stop, death)~albumin + esrdtime + bmi + factor(smokegrp)
                + cholest + diabetes + pst.sbp + age + hist.cvd, data = usrds.texpand)
summary(fit.tv)
results.tab(fit.tv)

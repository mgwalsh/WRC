# 2-compartment soil water-retention curve model with random soil/site effects
# note that this is a starter setup using partial data from central Kenya, using the <nlme> package
# M. Walsh, November 2018

# Required packages
# install.packages(c("downloader","nlme","lattice")), dependencies=TRUE)
suppressPackageStartupMessages({
  require(downloader)
  require(nlme)
  require(lattice)
})

# Data downloads -----------------------------------------------------------
# set working directory
dir.create("KE_wrc", showWarnings = F)
setwd("./KE_wrc")

# download data
download("https://www.dropbox.com/s/jgas6v78z5cfkks/KE_test.csv?raw=1", "KE_test.csv", mode = "wb")
wrc <- read.table("KE_test.csv", header = T, sep = ",") ## load "tidy" (long) version of the original data

# Complete pooling model <nls> --------------------------------------------
vwc.nls <- nls(vwc~SSbiexp(pf,s1,r1,s2,r2), data=wrc) 
summary(vwc.nls)
plot(vwc.nls, resid(.) ~ fitted(.) | sid, abline = 0, grid=F)

# No pooling model <nlsList> ----------------------------------------------
gwrc <- groupedData(vwc ~ pf | sid, as.data.frame(wrc))
vwc.lis <- nlsList(vwc~SSbiexp(pf,s1,r1,s2,r2), gwrc) 
plot(augPred(vwc.lis, level=0:1), xlab ="pF (bars)", ylab = "Volumetric water content") ## plot of site/sid level fits
plot(intervals(vwc.lis))

# Random effects model <nlme> ---------------------------------------------
init <- getInitial(vwc~SSbiexp(pf,s1,r1,s2,r2), gwrc) 
vwc.nlme <- nlme(vwc~SSbiexp(pf,s1,r1,s2,r2), data = gwrc, 
                 fixed = s1+r1+s2+r2~1, 
                 random = pdDiag(s1+r1+s2+r2~1), 
                 start = c(s1 = init[1], r1 = init[2], s2 = init[3], r2 = init[4])) 
summary(vwc.nlme) 
plot(augPred(vwc.nlme, level=0:1), xlab ="pF (bars)", ylab = "Volumetric water content") ## plot of site/sid level fits

# plot of overall <nlme> fit
par(pty="s")
plot(gwrc$vwc ~ fitted(vwc.nlme), xlim = c(0,0.8), ylim = c(0,0.8), xlab = "fitted", ylab = "measured")
abline(c(0,1))

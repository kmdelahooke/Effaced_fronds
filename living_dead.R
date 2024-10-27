##TESTING PRE-BURIAL MORTALITY

source("./surface_data_processing.R") # point_patterns, pp_eff, pp_well
source("./raster_processing_functions.R")

# Libraries
library(circular)
library(CircMLE)
library(CircStats)
library(parallel)

#-------------------------------------------------------------------------------
# (1) ORIENTATION DIFFERENCES
#-------------------------------------------------------------------------------

#replicate without fractos
#well <- lapply(well, subset, !grepl("fracto", sp, ignore.case = T))

# account for St Shotts formatting differences
well[[14]]$FrondA <- well[[14]]$StemA
eff[[14]]$FrondA <- eff[[14]]$StemA

# sample size
n_w <- unlist(lapply(well, nrow))
n_e <- unlist(lapply(eff, nrow))
write.csv(cbind.data.frame(n_w, n_e), "./sample_size_nofracto.csv")

# calculate mean orientation and plot orientations
plot_angles <- function(x, colour, title){
  angles <- circular(x$FrondA[x$FrondA != 0.0000], units='degrees', rotation='clock', zero=pi/2)
  plot.circular(angles, col=colour, stack=T, shrink = 1, main = title)
  arrows.circular(mean(angles, na.rm =T))
  return(mean(angles, na.rm = T))
}

means <- list()
for(i in 1:14){
  wa <- plot_angles(well[[i]], "cadetblue3", gsub(".csv", "", filelist[[i]]))
  wb<- plot_angles(eff[[i]], "antiquewhite3", gsub(".csv", "", filelist[[i]]))
  means[[i]]<- cbind(well = as.numeric(wa), eff = as.numeric(wb))
}
means <- plyr::ldply(means)
write.csv(means, "./circular_means_no_fracto.csv")

# test homogeneity
for(i in 1:14){
  angles_w <- circular(well[[i]]$FrondA[well[[i]]$FrondA != 0.0000], units='degrees', rotation='clock', zero=pi/2)
  angles_e <- circular(eff[[i]]$FrondA[eff[[i]]$FrondA != 0.0000], units='degrees', rotation='clock', zero=pi/2)
  print(filenames[i])
  rao.homogeneity(list(angles_w, angles_e))
}

# test uniformity: rayleigh and hermans rasson tests
angles_w <- list()
angles_e <- list()
for(i in 1:14){
  angles_w[[i]] <- circular(well[[i]]$FrondA[well[[i]]$FrondA != 0.0000], units='degrees', rotation='clock', zero=pi/2)
  angles_e[[i]] <- circular(eff[[i]]$FrondA[eff[[i]]$FrondA != 0.0000], units='degrees', rotation='clock', zero=pi/2)
 
}

raleigh_w <- plyr::ldply(lapply(angles_w, function(x){unlist(rayleigh.test(x)[1:2])}))
raleigh_e <- plyr::ldply(lapply(angles_e, function(x){unlist(rayleigh.test(x)[1:2])}))

rayleigh.res <- cbind.data.frame(raleigh_w, raleigh_e, filenames)
write.csv(rayleigh.res, "./rayleigh_well_eff_nofracto.csv")

herman_w <- plyr::ldply(lapply(angles_w, function(x){HR_test(x, iter = 999)}))
herman_e <- plyr::ldply(lapply(angles_e, function(x){HR_test(x, iter = 999)}))

herman.res <- cbind.data.frame(herman_w, herman_e, filenames)
write.csv(herman_w, "./herman_well_nofracto.csv")

# sample size angles
n_w <- unlist(lapply(angles_w, length))
n_e <- unlist(lapply(angles_e, length))
write.csv(cbind.data.frame(n_w, n_e), "./n_angles_nofracto.csv")

# test concentration
rho_w <- plyr::ldply(lapply(angles_w, rho.circular))
rho_e <- plyr::ldply(lapply(angles_e, rho.circular))

rho.res <- cbind.data.frame(rho_w,rho_e, filenames)
write.csv(rho.res, "./rho_well_eff_nofracto.csv", row.names = F)

# tests of distributional differences - Watson's U2
watson_test_comp <- function(well, eff){
  angles_w <- circular(well$FrondA[well$FrondA != 0.0000], units='degrees', rotation='clock', zero=pi/2)
  angles_e <- circular(eff$FrondA[eff$FrondA != 0.0000], units='degrees', rotation='clock', zero=pi/2)
  
  group <- as.factor(c(rep("well", length(angles_w)), rep("eff", length(angles_e))))
  angles <- c(angles_w, angles_e)
  #res <- watson.williams.test(angles, group, units = 'degrees', na.rm = T)
  #res <- watson.wheeler.test(list(angles_w, angles_e), na.rm = T)
  res <- watson.two.test(angles_w, angles_e)
  #return(cbind.data.frame(stat = res$statistic, p.value = res$alpha))
  print(res)
  return(res)
}

res.wat <- plyr::ldply(mapply(watson_test_comp, well, eff, SIMPLIFY = F))
res.wat$surface <- gsub(".csv", "", filelist)
str(res.wat)
res.wat$p.value[res.wat$p.value > 0.05]

write.csv(res.wat, "./watsonu2.csv")

#Watson's large sample non-parametric test from (Pewsey et al., 2013)
YgVal <- function(cdat, ndat, g) {
  N <- length(cdat) ; ndatcsum <- cumsum(ndat)
  delhat <- 0 ; tbar <- 0
  for (k in 1:g) {
    sample <- circular(0)
    if (k==1) {low <- 0} else
      if (k > 1) {low <- ndatcsum[k-1]}
    for (j in 1:ndat[k]) { sample[j] <- cdat[j+low] }
    tm1 <- trigonometric.moment(sample, p=1)
    tm2 <- trigonometric.moment(sample, p=2)
    Rbar1 <- tm1$rho; Rbar2 <- tm2$rho ; tbar[k] <- tm1$mu
    delhat[k] <- (1-Rbar2)/(2*Rbar1*Rbar1)
  }
  
  dhatmax <- max(delhat) ; dhatmin <- min(delhat)
  if (dhatmax/dhatmin <= 4) {
    CP <- 0 ; SP <- 0 ; dhat0 <- 0
    for (k in 1:g) {
      CP <- CP+ndat[k]*cos(tbar[k]) ; SP <- SP+ndat[k]*sin(tbar[k])
      dhat0 <- dhat0+ndat[k]*delhat[k] }
    dhat0 <- dhat0/N
    RP <- sqrt(CP*CP+SP*SP) ; Yg <- 2*(N-RP)/dhat0
    return(Yg) } else
      if (dhatmax/dhatmin > 4) {
        CM <- 0 ; SM <- 0 ; Yg <- 0
        for (k in 1:g) {
          CM <- CM+(ndat[k]*cos(tbar[k])/delhat[k])
          SM <- SM+(ndat[k]*sin(tbar[k])/delhat[k])
          Yg <- Yg+(ndat[k]/delhat[k]) }
        RM <- sqrt(CM*CM+SM*SM) ; Yg <- 2*(Yg-RM)
        return(Yg) }
}

Yg_comp <- function(well, eff){
  YgObs <- YgVal(cdat = c(well, eff), ndat = c(length(well), length(eff)), g = 2)
  res <- c(YgObs, pchisq(YgObs, 1, lower.tail=F))
  names(res)<- c("stat", "p-value")
  return(res)
}
res.large_wat <- plyr::ldply(mapply(Yg_comp, angles_w, angles_e, SIMPLIFY = F))
write.csv(res.large_wat, "./large_wat_nofracto.csv")

#-------------------------------------------------------------------------------
# (2) DENSITY DEPENDENCE RLA
#-------------------------------------------------------------------------------

#density dependence difference function
gdif <- function(X, ..., i, j){
  gidot <- pcfdot(X, ..., i=i)
  gjdot <- pcfdot(X, ..., i=j)
  #g <- pcf(X, ...)
  dif <- eval.fv(gidot - gjdot)
  return(dif)
}

# standard envelope
#pd_dendep <- vector("numeric", 14)
#pd_dendep2 <- vector("numeric", 14)
#par(mfrow = c(1,1))
#for(i in 1:14){
#  E <- envelope(pp_mark[[i]], gdif, i = "eff", j = "well", stoyan = 0.2, simulate=expression(rlabel(pp_mark[[i]])), savefuns = T, nsim = 19)
#  plot(E, main = gsub(".csv", "", filenames[i]), legend = F )
#  pd_dendep2[i]  <- dclf.test(E, rinterval = c(50, (round(max(E$r)) - 1)))$p.value
#  pd_dendep[i] <- dclf.test(E)$p.value
#}

#write.csv(cbind.data.frame(filenames, pd_dendep, pd_dendep2), "./pd_dendep_n19.csv", row.names = F)

# parallelised envelope for large numbers of simulations

cl <- makeCluster(detectCores())
clusterEvalQ(cl, {library(spatstat)})

pd_dendep <- vector("numeric", 14)
pd_dendep2 <- vector("numeric", 14)
for(i in 1:14){
  
  nsim <- round(9999/detectCores()) # number of simulations to run on each core
  
  call_envelope <- function(x, nsim){
    gdif <- function(X, ..., i, j){
      gidot <- pcfdot(X, ..., i=i)
      gjdot <- pcfdot(X, ..., i=j)
      dif <- eval.fv(gidot - gjdot)
      return(dif)
    }
    E <-  spatstat.explore::envelope(x, gdif, i = "eff", j = "well", stoyan = 0.2, simulate = expression(rlabel(x)), savepatterns = T, savefuns = T, nsim = nsim)
    return(E)
  }
  
  pp <- pp_mark[[i]]
  
  ppplist <- replicate(detectCores(), pp, simplify = FALSE)
  envlist <- parLapply(cl = cl, X = ppplist, fun = call_envelope, nsim = nsim) # separate function to avoid specifying multiple 'fun' arguments
  envfinal <- do.call(pool, c(envlist, savepatterns = T, savefuns = T)) # pool envelope
  E <- envelope.envelope(envfinal, nrank = 500, nsim = 9999, savefuns = T) # recompute envelope to have desired nsim/ rank
  plot(E, main = gsub(".csv", "", filenames[i]), legend = F )
  
  print(filenames[i])   #tracking
  
  #record goodness-of-fit
  pd_dendep[i] <- dclf.test(E)$p.value
  pd_dendep2[i] <- dclf.test(E,  rinterval = c(5, (round(max(E$r)) - 1)))$p.value
}

stopCluster(cl)

#df <- cbind.data.frame(pd_dendep[pd_dendep!=0.0000], pd_dendep2[pd_dendep2!=0.0000])
write.csv(cbind.data.frame(filenames, pd_dendep, pd_dendep2), "./pd_dendep_new_n9999.csv", row.names = F)


#-------------------------------------------------------------------------------
# (2) UNIVARIATE AND BIVARIATE RLA
#-------------------------------------------------------------------------------

cl <- makeCluster(detectCores())
clusterEvalQ(cl, {library(spatstat)})
pd_rl <- vector("numeric", 14)
pd_rl2 <- vector("numeric", 14)
for(i in 1:14){
  
  nsim <- round(9999/detectCores()) # number of simulations to run on each core
  
  call_envelope <- function(x, nsim){
    E <-  spatstat.explore::envelope(x, pcfcross, i = "well",j = "well", stoyan = 0.2, simulate = expression(rlabel(x)), savepatterns = T, savefuns = T, nsim = nsim)
    return(E)
  }
  
  pp <- pp_mark[[i]]
  
  ppplist <- replicate(detectCores(), pp, simplify = FALSE)
  envlist <- parLapply(cl = cl, X = ppplist, fun = call_envelope, nsim = nsim) # separate function to avoid specifying multiple 'fun' arguments
  envfinal <- do.call(pool, c(envlist, savepatterns = T, savefuns = T)) # pool envelope
  E <- envelope.envelope(envfinal, nrank = 500, nsim = 9999, savefuns = T) # recompute envelope to have desired nsim/ rank
  plot(E, main = gsub(".csv", "", filenames[i]), legend = F )
  pd_rl[i] <- dclf.test(E)$p.value
  pd_rl2[i] <- dclf.test(E,  rinterval = c(5, (round(max(E$r)) - 1)))$p.value
}
stopCluster(cl)

write.csv(cbind.data.frame(filenames, pd_rl, pd_rl2), "./pd_rl_well2_n9999.csv", row.names = F)


#-------------------------------------------------------------------------------
# (3) EFFACEMENT VS W-STAT
#-------------------------------------------------------------------------------
w_stat <- read.csv("./w-stat.csv")
str(w_stat)

#%effacement
p_eff <- vector("numeric", 14)
for(i in 1:14){
  p_eff[i] <- length(eff[[i]][,1])/(length(well[[i]][,1])+length(eff[[i]][,1]))*100
}

#remove cpc and MUN
df <- cbind(p_eff = p_eff[-c(3,12)], w_stat)

#regression
fit1 <- lm(p_eff ~ w_stat, df)

# check residuals
par(mfrow = c(2,2))
plot(fit1)

par(mfrow = c(1,1))
plot(p_eff ~ w_stat, df)
abline(fit1)

hist(df$w_stat)

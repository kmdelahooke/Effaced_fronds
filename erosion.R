#EROSION AND EFFACEMENT

# Dependencies
source("./surface_data_processing.R") # point_patterns, pp_eff, pp_well
source("./raster_processing_functions.R")
#source("./spatial_analysis_functions.R") # called in 'surface_data_processing.R'

# Additional Libraries
library(dplyr)
library(plyr)
library(ggpubr)

#-------------------------------------------------------------------------------
# (1) kernel smoothed intensity visualization
#-------------------------------------------------------------------------------

dens.well <- lapply(pp_well, density.ppp, sigma = 1000)
dens.eff <- lapply(pp_eff, density.ppp, sigma = 1000)

par(mfrow = c(3,2))
for (i in 1:length(dens.well)){
  plot(dens.well[[i]], main = paste(gsub(".csv", "", filelist[[i]]), "well-preserved n =", pp_well[[i]]$n))
  plot(dens.eff[[i]], main = paste(gsub(".csv", "", filelist[[i]]), "effaced n =", pp_eff[[i]]$n))
}

#-------------------------------------------------------------------------------
# (2) Plot point patterns on same map
#-------------------------------------------------------------------------------
par(mfrow = c(1,1))

for (i in 1:length(dens.well)){
  plot(pp_well[[i]], main = paste(gsub(".csv", "", filelist[[i]])), cols = "cadetblue3", pch = 20)
  plot(pp_eff[[i]], main = paste(gsub(".csv", "", filelist[[i]])), cols = "antiquewhite3", pch = 20, add = T)
}

#-------------------------------------------------------------------------------
# (4) Anisotropy - point pair-orientation distribution
#-------------------------------------------------------------------------------
par(mfrow = c(3,2))
for (i in 1:length(pp_well)){
  rose(pairorient(pp_well[[i]], r1=10, r2=1000, sigma = 10), main = paste(gsub(".csv", "", filelist[[i]]), "well-preserved n =", pp_well[[i]]$n))
  rose(pairorient(pp_eff[[i]], r1=10, r2=1000, sigma = 10), main = paste(gsub(".csv", "", filelist[[i]]), "effaced n =", pp_eff[[i]]$n))
}

# Anisotropy individual taxa
fracto <- subset(surface_data[[10]], sp == "parviscopa")
fpp <- ppp(fracto$x, fracto$y, window = windows[[10]])
par(mfrow = c(1,2))
plot(fpp)

rose(pairorient(fpp, r1 = 10, r2 = 200, sigma = 10))

#-------------------------------------------------------------------------------
# (5) Spatial autocorrelation
#-------------------------------------------------------------------------------
par(mfrow = c(3,3))

# Convert im to raster
rasts.well <- lapply(dens.well, terra::rast)
rasts.eff <- lapply(dens.eff, terra::rast)

for(i in 1:length(rasts.well)){
  names(rasts.well[[i]]) <- gsub(".csv", "", filelist[[i]])
  names(rasts.eff[[i]]) <- gsub(".csv", "", filelist[[i]])
}

# Plot correlogram (adapted from pgirmess 'plot.correlog')
plot_moran_correlog <- function(x){
  name <- names(x)
  coords <- crds(terra::aggregate(x, 10))
  distmat <- as.matrix(dist(coords))
  maxdist <- 1/2*max(distmat)
  
  x <- as.df.Spatraster(x, 30)
  
  clg <- pgirmess::correlog(x[,1:2], x[,3], nbclass = 30)
  pbpv <- rep(0.05, 30)/seq(1, 30, 1) # progressive bonferroni correction
  
  plot(clg[,1:2,drop = FALSE], type="b", xlab="distance classes", ylab = attributes(clg)$Method, main = name)
  inc<-(clg[2,1] - clg[1,1])/2
  breaks <- pretty(c(clg[1,1] - inc, clg[length(clg[,1]),1] + inc), n = length(clg[,1]), min.n = 2)
  axis(1,at=breaks)
  points(clg[clg[,3]< pbpv,1:2,drop = FALSE], pch = 19, col = "red",cex = 2)
  abline(h = 0)
  
}

lapply(rasts.well, plot_moran_correlog)
lapply(rasts.eff, plot_moran_correlog)


#-------------------------------------------------------------------------------
# (6) Relationship between well-preserved and effaced fronds
#-------------------------------------------------------------------------------

#fit heterogeneous poisson model of effaced to well preserved fronds
#AIC or ANOVA (select where appropriate) of fitted model vs CSR

fit_hp <- function(pp, den, name){
  hpp <- ppm(pp ~ den)
  #plot(effectfun(hpp), main = name, xlab = "well-preserved", ylab = "effaced")
  im <- list(dens.well[[i]])
  names(im) <- "well"
  fit0 <- ppm(pp_eff[[i]] ~ 1) # CSR
  fit1 <- ppm(pp_eff[[i]] ~ well, data = im)
  
  #res <- anova(fit0, fit1, test ="LR")
  res <- c(AIC(fit0), AIC(fit1))
  return(res)
}

par(mfrow = c(3,3))
hpp <- vector("list", 13)
for(i in 1:length(dens.well)){
  hpp[[i]] <- fit_hp(pp_eff[[i]], dens.well[[i]], gsub(".csv", "", filelist[[i]]))
}

# Evaluate fitted model
par(mfrow = c(1,1))
lapply(hpp, diagnose.ppm)
lapply(hpp, qqplot.ppm)


#Export results
res <- plyr::ldply(hpp)
names(res) <- c('csr', 'well')
res$surface <- rep(gsub(".csv", "", filelist), each = 2)

write.csv(res, "./eff_well_aic_raw.csv")

#-------------------------------------------------------------------------------
# (7) Relationship of effaced and well-preserved specimens with cartesian covariates
#-------------------------------------------------------------------------------

#A. AIC or ANOVA of fitted heterogeneous Poisson models
#-------------------------------------------------------------------------------
#Plots of relationship
#select as appropriate

res <- list()
for(i in 1:length(pp_well)){
  #pp <- pp_well[[i]]
  pp2 <- pp_eff[[i]]
  den <- density(pp)
  
  #Point of max/min erosion
  dr <- terra::rast(den)
  max <- xyFromCell(dr, where.max(dr)[1,2]) 
  min <- xyFromCell(dr, where.min(dr)[1,2])
  
  a <- max[1]
  b <- max[2]
  c <- min[1]
  d <- min[2]
  
  fmax <- function(x,y) {sqrt((x-a)^2 +(y-b)^2)}
  fmin <- function(x,y) {sqrt((x-c)^2 +(y-d)^2)}
  
  #HPP (well-preserved)
  #csr <- ppm(pp ~ 1)
  #hx <- ppm(pp ~ x)
  #hy <- ppm(pp ~ y)
  #hmin <- ppm(pp ~ fmax)
  #hmax <- ppm(pp ~ fmin)
  
  #HPP2 (effaced)
  csr2 <- ppm(pp2 ~ 1)
  hx2 <- ppm(pp2 ~ x)
  hy2 <- ppm(pp2 ~ y)
  hmin2 <- ppm(pp2 ~ fmax)
  hmax2 <- ppm(pp2 ~ fmin)
  
  #AIC
  csr <- AIC(csr2)
  x <- AIC(hx2)
  y <- AIC(hy2)
  min <- AIC(hmin2)
  max <- AIC(hmax2)
  
  df <- cbind.data.frame(csr,x,y,min,max)
  
  #Goodness of Fit
  #csr <- dclf.test(csr2)
  #x <- dclf.test(hx2)
  #y <- dclf.test(hy2)
  #min <- dclf.test(hmin2)
  #max <- dclf.test(hmax2)
  
  #df <- cbind.data.frame(csr,x,y,min,max)

  #anova (well-preserved)
  #resx <- anova(csr, hx, test ="LR")[2,4]
  #resy <- anova(csr, hy, test ="LR")[2,4]
  #resmin <- anova(csr, hmin, test ="LR")[2,4]
  #resmax <- anova(csr, hmax, test ="LR")[2,4]
  
  #anova (effaced)
  #resx2 <- anova(csr2, hx2, test ="LR")[2,4]
  #resy2 <- anova(csr2, hy2, test ="LR")[2,4]
  #resmin2 <- anova(csr2, hmin2, test ="LR")[2,4]
  #resmax2 <- anova(csr2, hmax2, test ="LR")[2,4]
  
  #df <- cbind.data.frame(resx,resy,resmin,resmax,resx2,resy2,resmin2,resmax2)
  
  # write results
  res[[i]] <- df
  
  # Plot effect functions
  #par(mfrow = c(2,2))
  
  #maxy <- max(c(effectfun(hx)$lambda, effectfun(hx2)$lambda))
  #miny <- min(c(effectfun(hx)$lambda, effectfun(hx2)$lambda))
  #plot(effectfun(hx), main = paste(gsub(".csv", "", filelist[[i]]),": well", round(resx, 10), "eff",round(resx2, 10)), xlab = "x-coordinate", ylab = "intensity well preserved", ylim = c(miny, maxy), lwd = 2, col = "cadetblue3")
  #plot(effectfun(hx2), col = "antiquewhite4", add = T, lwd = 2)
  
  #maxy <- max(c(effectfun(hy)$lambda, effectfun(hy2)$lambda))
  #miny <- min(c(effectfun(hy)$lambda, effectfun(hy2)$lambda))
  #plot(effectfun(hy), main = paste(gsub(".csv", "", filelist[[i]]),": well", round(resy, 10), "eff",round(resy2, 10)), xlab = "y-coordinate", ylab = "intensity well preserved", ylim = c(miny, maxy), lwd = 2, col = "cadetblue3")
  #plot(effectfun(hy2), col = "antiquewhite4", add = T, lwd = 2)
  
  #maxy <- max(c(effectfun(hmin)$lambda, effectfun(hmin2)$lambda))
  #miny <- min(c(effectfun(hmin)$lambda, effectfun(hmin2)$lambda))
  #plot(effectfun(hmin), main = paste(gsub(".csv", "", filelist[[i]]),": well", round(resmin, 10), "eff",round(resmin2, 10)), xlab = "distance from minimum", ylab = "intensity well preserved", ylim = c(miny, maxy), lwd = 2, col = "cadetblue3")
  #plot(effectfun(hmin2), col = "antiquewhite4", add = T, lwd = 2)
  
  #maxy <- max(c(effectfun(hmax)$lambda, effectfun(hmax2)$lambda))
  #miny <- min(c(effectfun(hmax)$lambda, effectfun(hmax2)$lambda))
  #plot(effectfun(hmax), main = paste(gsub(".csv", "", filelist[[i]]),": well", round(resmax, 10), "eff",round(resmax2, 10)), xlab = "distance from maximum", ylab = "intensity well preserved", ylim = c(miny, maxy), lwd = 2, col = "cadetblue3")
  #plot(effectfun(hmax2), col = "antiquewhite4", add = T, lwd = 2)

}

# Export results
res <- plyr::ldply(res)
write.csv(res, "./aic_eff_raw_cartesian.csv")

#  With custom covariates
#-------------------------------------------------------------------------------
par(mfrow = c(1,1))

# Specify covariates

# E surface distance from line of the stream
plot(well[[5]]$y ~ well[[5]]$x ) # change as necessary
coords <- locator()

# Custom line
plot(density(pp_well[[5]])) # change as necessary
coords <- locator()

# Save coordinates

#dips
bd <- coords
cd <- coords
ed <- coords
gd <- coords
ld <- coords
md <- coords

#temp
sd <- coords

#stream
es <- coords

# Export coordinates
convert_coords <- function(X){
  df <- cbind.data.frame(x1 = X$x[1],x2 = X$x[2], y1 = X$y[1], y2 = X$y[2] )
  return(df)
}

coords <- plyr::ldply(lapply(list(bd, cd, ed, gd, ld, md), convert_coords))
coords$surface <- c("bishop", "cpc", "E", "goldmine", "lmp", "mun")
write.csv(coords, "./dip_coordinates_raw.csv")

coords <- plyr::ldply(lapply(list(sd, es), convert_coords))
coords$surface <- c("e_stream", "ss_temp")
write.csv(coords, "./other_coordinates.csv")

# Incorporate dip into heterogeneous Poisson models

coords <- md # mun dip for example

dip_direction <- function(x, y){
  #y <- dip$coefficients[2]*x + dip$coefficients[1] #line
  #x2 <- dip$coefficients[2]*y + dip$coefficients[1] # perpendicular line
  a <- -(coords$x[2] - coords$x[1])/(coords$y[2]- coords$y[1])
  b <- -1
  c <- -a*coords$x[1] + coords$y[1]
  d <- abs(a*x + b*y + c)/sqrt(a^2 + b^2)
  return(d)
  
}

hpp <- ppm(pp_well[[12]]~ dip_direction + I(dip_direction^2)) # quadratic mun ~ dip model
#hpp <- ppm(pp_well[[12]]~ dip_direction ) # linear
plot(predict(hpp)) # check that this makes sense
AIC(hpp)

#csr null models
hpps <- lapply(pp_well, ppm, ~1)
aics <- lapply(hpps, AIC)
write.csv(cbind.data.frame(filenames, unlist(aics)), "./aic_csr_nofracto.csv", row.names = F)


#Incorporate distance line to point into HP models

coords <- es #  E stream for example

a <- (coords$y[1] - coords$y[2])/(coords$x[1]- coords$x[2])
b <- -1
c <- -a*coords$x[1] + coords$y[1]

distance_line_to_point <- function(x, y){
  #a <- (coords$x[2] - coords$x[1])/(coords$y[2]- coords$y[1])
  #b <- -1
  #c <- -a*coords$x[1] + coords$y[1]
  d <- abs(a*x + b*y + c)/sqrt(a^2 + b^2)
  return(d)
}

hpp <- ppm(pp_well[[5]]~ distance_line_to_point + I(distance_line_to_point ^2)) #quadratic
#hpp <- ppm(pp_well[[5]]~ distance_line_to_point) #linear
plot(predict(hpp)) #check correct prediction
AIC(hpp)

# Incorporate distance to point into HP models

par(mfrow = c(1,1))
plot(well[[14]]$y ~well[[14]]$x ) # st shotts
coords <- locator() #select point

plot(dens.well[[14]]) #plot
points(coords, col = "red", cex = 2, pch = 4)


a = coords$x
b = coords$y
hpp <- ppm(pp_well[[14]] ~ sqrt((x-a)^2 +(y-b)^2) + I(sqrt((x-a)^2 +(y-b)^2)^2)) # quadratic
#hpp <- ppm(pp_well[[14]] ~ sqrt((x-a)^2 +(y-b)^2)) #linear
plot(predict(hpp))
AIC(hpp)


## Cropped E surface
plot(density(pp_mark[[5]]))
plot(pp_mark[[5]]$y ~pp_mark[[5]]$x, pch = 20, col = "skyblue")
win <- clickpoly(add = T)
plot(win, add = T)
wl <- lapply(well, subset, !grepl("fracto", sp, ignore.case = T))
em <- ppp(well[[5]]$x, well[[5]]$y, win)
plot(em, pch = 20, cols = "skyblue")


# B. Kolmogorov-Smirnov tests
#-------------------------------------------------------------------------------
df <- lapply(pp_well, erosion_bias) # erosion_bias() from 'spatial_analysis_functions.R'

# isolate aspects of results
pv <- ldply(lapply(df, function(x){x$p_values}))
ts <- ldply(lapply(df, function(x){x$teststat}))

colnames(pv) <- c("x", "y", "min", "max")
colnames(ts) <- c("x", "y", "min", "max")

write.csv(pv, "./ks_pvalue_nofracto_hp.csv")
write.csv(ts, "./ks_stat_nofracto_hp.csv")


#-------------------------------------------------------------------------------
# (8) Size distribution comparison
#-------------------------------------------------------------------------------

#wl <- lapply(well, subset, !grepl("fracto", sp, ignore.case = T)) #repeat minus fractos

#SURFACE-WISE

# lengths
length_w <- lapply(well, function(x){x$FrondL + x$FrondW})
length_e <- lapply(eff, function(x){x$FrondL + x$FrondW})

#normality
norm <- lapply(length_e, function(x){c(shapiro.test(log(x))$p.value, shapiro.test(log(x))$statistic)})
norm <- plyr::ldply(norm)
names(norm) <- c("p.value", "W")

write.csv(norm, "./eff_shapiro_lengths_raw_log.csv")

# histograms
plot_hist <- function(v1, v2, name){
  #v1 <- v1$FrondL + v1$StemL
  #v2 <- v2$FrondL + v2$StemL
  maxf <- as.numeric(max(c(hist(v1, plot = F)$count, hist(v2, plot = F)$count)))
  hist(v1, col = "cadetblue3", main = name, xlim = c(0, max(c(v1,v2))), ylim = c(0, maxf), plot = T)
  hist(v2, col = "antiquewhite3", add = T)
}

mapply(plot_hist, length_w, length_e, name = filenames)

# size differnces
size_differences <- function(x, y){
  mark <- c(rep("well", length(x)),rep("eff", length(y)) )
  length <- c(unlist(x), unlist(y))
  df <- cbind.data.frame(length, mark)
  df$mark <- as.factor(df$mark)
  df$length <- as.numeric(df$length)
  
  #treatments
  df <- subset(df, length < 600) #remove outlier
  df <- subset(df, !is.infinite(length))
  #df$length <- log(df$length)
  
  #Equality of variance
  lev.res <- car::leveneTest(length ~ mark, data = df)
  bar.res <- bartlett.test(length ~ mark, data = df)
  
  #welch t-test
  w.res <- t.test(length ~mark, data = df, alternative = "two.sided", var.equal = F)
  
  # Kruskall wallis test
  #k.res <- kruskal.test(length ~mark, data = df)
  #k.res.df <- cbind.data.frame(stat = k.res$statistic, k.p.value = k.res$p.value)
  
  res.df <- cbind.data.frame(n.well = length(subset(df, mark == "well")[,1]), n.eff = length(subset(df, mark == "eff")[,1]),
                             lev.p.value = lev.res[1,3], F.value = lev.res[1,2],
                             bar.p.value = bar.res$p.value, k2 = bar.res$statistic, 
                             t.p.value = w.res$p.value, df = w.res$parameter, t = w.res$statistic,
                             mean.well = w.res$estimate[2], mean.eff = w.res$estimate[1])
  
  return(res.df)
  #return(k.res.df)
  
}
res <- plyr::ldply(mapply(size_differences, length_w, length_e, SIMPLIFY = F))
res$surface <- filenames

write.csv(res, "./surface_size_diff_retro.csv", row.names = F) #Export

#COMBINED DATASET

#A: scale to maximum of group on surface
scaled_size <- function(x){
  size <- x$FrondL + x$FrondW
  #size <- size[size < 600] #remove outlier
  return(size/max(size))
}
length_w <- lapply(well, scaled_size)
length_e <- lapply(eff, scaled_size)

#B: scale by maximum of combined on surface
sf <- mapply(function(v1, v2){return(max(v1,v2))}, length_w, length_e)
length_w <- mapply(function(x, sf){x/sf}, length_w, sf)
length_e <- mapply(function(x, sf){x/sf}, length_e, sf)

#combine
mark <- c(rep("well", length(unlist(length_w))),rep("eff", length(unlist(length_e))) )
length <- c(unlist(length_w), unlist(length_e))

df <- cbind.data.frame(length, mark)
df$mark <- as.factor(df$mark)
df$length <- as.numeric(df$length)

#treatment options
df <- subset(df, length < 600) # remove outlier
df <- subset(df, !is.infinite(length))
df$length <- log(df$length)

#test normality
shapiro.test(sample(df$length, 5000))
shapiro.test(sample(unlist(length_w), 5000))

#equality of variance
car::leveneTest(length ~ mark, data = df)
bartlett.test(length ~ mark, data = df)

#box plots
ggboxplot(df, x = "mark", y = "length", xlab = "type", ylab = "length (mm)")

# welch t-test
t.test(length ~ mark, data = df, alternative = "two.sided", var.equal = F)

# kruskal-wallis test
kruskal.test(length ~mark, data = df)

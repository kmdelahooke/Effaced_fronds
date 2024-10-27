#SYN-FOSSILISATION PROCESSES: COMPARING SCALE OF HETEROGENEITY

source("./surface_data_processing.R") # point_patterns, pp_eff, pp_well
source("./raster_processing_functions.R")

#Libraries
library(zoo)
library(spatstat)

#-------------------------------------------------------------------------------
# (1) scale of heterogeneity matgrounds
#-------------------------------------------------------------------------------

# processed texture images
filelist <- list.files("D:/postprocessed_ti_w2/", ".*.tif")
filenames <- gsub(".tif", "", filelist)
filepaths <- unlist(lapply("D:/postprocessed_ti_w2/", paste0, filelist))

maps95 <- lapply(filepaths, terra::rast)
maps95 <- lapply(maps95, terra::crop, ext(1800, 2300, 1000, 1500))
maps95 <- maps95[-c(1,6,8,15,16,17,18,20,21)] # remove those with null data


calculate_cor_length <- function(rast, res){
  
  pc_dem <- terra::aggregate(rast, res, na.rm = T) 
  
  #pc_dem <- remove_outliers(pc_dem, 0.999) 
  pc_dem <- interpolate_raster(pc_dem, 5) # interpolate cracks

  scl_res <- geodiv::scl(raster(pc_dem))[1]
  return(scl_res)
}

scl_res <- list()
for(i in c(1:8, 10,11,12,14)){
  scl_res[[i]] <- calculate_cor_length(maps95[[i]], 5)

}

write.csv(cbind.data.frame(filenames[c(1:8, 10,11,12,14)], unlist(scl_res)), "./scl.csv")

scl_res <- lapply(maps95, FUN = calculate_cor_length, res = 5)

#plot(unlist(scl_res) ~ res, ylab = "Scl", xlab = "aggregation factor")

#-------------------------------------------------------------------------------
# (2) Scale of heterogeneity fronds - min nn cross distance
#-------------------------------------------------------------------------------

# Visualisation: combine KDEs

mean_nndist <- unlist(lapply(pp_mark, function(x){mean(nndist(unmark(x)))})) # sigma as the mean interpoint distance
length(mean_nndist)
for(i in 1:14){
  dw <- density.ppp(pp_well[[2]], sigma = round(mean_nndist[i]), eps = 15) 
  de <-density.ppp(pp_eff[[2]], sigma = round(mean_nndist[i]), eps = 15) 
  
  
  cdem <- lapp(sds(rast(dw), rast(de)), function(x, y){return(x - y)})
  par(mfrow = c(1,1))
  #plot(rast(dw), main = "well_preserved")
  #plot(rast(de), main = "effaced")
  plot(cdem, main = paste(filenames[i], "well preserved - effaced"))
  plot(pp_well[[i]], pch = 20, col = "springgreen4", add = T)
  plot(pp_eff[[i]], pch = 20, col = "lightyellow1", add = T)
}


# Measuring scale of heterogeneity: closest effaced frond to each well-preserved frond

#wl <- lapply(well, subset, !grepl("fracto", sp, ignore.case = T)) # repeat without fractos
pp <-list()
for(i in 1:14){
  pp[[i]] <- ppp(wl[[i]]$x, wl[[i]]$y, windows[[i]])
}


min_cross_dist <- function(x, y){
  dists <- crossdist(x, y)
  
  dists1 <- vector('numeric', dim(dists)[1]) # distance well preserved to effaced
  for(i in 1:dim(dists)[1]){
    dists1[i] <- min(dists[i,])
  }
  
  dists2 <- vector('numeric', dim(dists)[2]) # distance effaced to well preserved
  for(i in 1:dim(dists)[2]){
    dists2[i] <- min(dists[,i])
  }

  df <- c(dists1, dists2)
  
  df <- df[!duplicated(df)] # remove diagonal duplicates
  return(cbind(min = min(df), max = max(df), median = median(df), quart_1 = summary(df)[2]))
  #return(cbind(min = min(dists2), max = max(dists2), median = median(dists2), quart_1 = summary(dists2)[2]))

}

dists <- mapply(min_cross_dist, pp, pp_eff)
dists <- as.data.frame(t(dists))
names(dists) <- c("min", "max", "median", "quart1")
dists$surface <- filenames
write.csv(dists, "./min_cross_dists_nofractos.csv") # Export

#boxplot(dists$median)
#hist(dists$median)

#-------------------------------------------------------------------------------
# (3) Scale of heterogeneity fronds - mark connection function
#-------------------------------------------------------------------------------

#difference function
par(mfrow = c(1,2))
for(i in 1:14){
  p_ij <- markconnect(pp_mark[[i]], "eff", "well", correction = "isotropic", normalise = T)
  p_ii <- markconnect(pp_mark[[i]], "eff", "eff", correction = "isotropic", normalise = T)
  
  dif <- eval.fv(p_ij - p_ii)
  plot(dif, main = filenames[i], legend = F)
  abline(h = 0, lty = 2, col = "blue")
  
  if(max(p_ii$iso) > max(p_ij$iso)){
    plot(p_ii$iso ~ p_ii$r, type = "l", ylab = "p_ii/p_ij", main = filenames[i], xlab = "distance r (mm)", ylim = c(0, max(p_ii$iso)))
    lines(p_ij$iso ~ p_ij$r, col = "red")
  }else{
    plot(p_ij$iso ~ p_ij$r, type = "l", ylab = "p_ii/p_ij", col = "red", main = filenames[i], xlab = "distance r (mm)", ylim = c(0, max(p_ij$iso)))
    lines(p_ii$iso ~ p_ii$r)
  }
  
  #find distance at which p >0
 min(p_ii[p_ii$iso > 0]$r)
  
}


scale_markconnect <- function(pp, i, j, ...){
  pij <- markconnect(pp, i, j, ...)
  pii <- markconnect(pp, i, i, ...)
  dif <- eval.fv(pij - pii)
  if(pii$iso[2] > pii$theo[1]){scale <- min(pii[pii$iso < pii$theo[1]]$r)}
  else{scale <- min(pii[pii$iso > pii$theo[1]]$r)}
  return(scale)
}

scale <-lapply(pp_mark, scale_markconnect, "eff", "well", correction = "isotropic", normalise = F)

scale <- cbind.data.frame(surface = filenames, scale = unlist(scale))
write.csv(scale, "./scale_pii_iso_nf.csv", row.names = F)



#-------------------------------------------------------------------------------
# (4) Level of heterogeneity fronds - LH*
#-------------------------------------------------------------------------------

level_heterogeneity <- function(pp, nsim){
  
  #nearest neighbour distance function

  Gw <- spatstat.explore::envelope(pp, Gest, correction = "km", nsim = 99)
  
  #remove infinite values
  Gw <- subset(Gw, is.finite(Gw$obs)) 
  
  #mean of sims
  gbar <- mapply(function(x,y){mean(c(x,y))}, Gw$lo, Gw$hi)
  
  #normalise
  wn <- Gw$r * 2*intensity(pp)^0.5
  
  # differences
  Gdiff <- abs(Gw$obs - Gw$theo)
  Gdiff_star <- abs(Gw$obs - gbar)
  
  #calculate integrals
  LH <- sum(diff(Gw$r) * rollmean(Gdiff, 2))
  LH_star <- sum(diff(Gw$r) * rollmean(Gdiff_star, 2))
  NLH_star <- sum(diff(wn) * rollmean(Gdiff_star, 2))
  
  #Maximum difference
  df <- cbind.data.frame(r = Gw$r, Gdiff)
  maxGdiff <- subset(df, Gdiff == max(Gdiff))
  df <- cbind.data.frame(r = Gw$r, Gdiff_star)
  maxGdiff_star <- subset(df, Gdiff == max(Gdiff))
  
  return(cbind(LH, LH_star, NLH_star, maxGdiff = maxGdiff$r, maxGdiff_star = maxGdiff_star$r))

}

lhs <- lapply(pp_well, level_heterogeneity, nsim = 99)
lhs <- plyr::ldply(lhs)
lhs$surface <- filenames
write.csv(lhs, "./level_het_n99_well_nofracto.csv", row.names = F)


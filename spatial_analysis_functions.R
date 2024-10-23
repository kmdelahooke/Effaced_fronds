#SPATIAL ANALYSIS FUNCTIONS

#includes:
# 1. make_windows
# 2. convert_angles
# 3. erosion_bias
# 4. erosion_aic

# Libraries
library(spatstat)
library(raster) # N.B. will need to update to terra for future functionality

#-------------------------------------------------------------------------------
# (1) Create spatstat owin windows from saved coordinates 
#-------------------------------------------------------------------------------

make_windows <- function(x){
  l <- list(x = as.numeric(x[,1]), y = as.numeric(x[,2]))
  win <- owin(poly = l)
  return(win)
}

#-------------------------------------------------------------------------------
# (2) convert angles from anticlockwise from x-axis to clockwise from y-axis
#-------------------------------------------------------------------------------

convert_angles <- function(x){
  ifelse(x < 90, {x = -x + 90}, {x = -x + 450})}


#-------------------------------------------------------------------------------
# (3) Kolmogorov-Smirnov tests of heterogeneous Poisson models with spatial 
# covariates
#-------------------------------------------------------------------------------

#select options as appropriate

erosion_bias <- function(pp) {
  pp <- unmark(pp)
  
  #Density map of everything on surface
  den <- density(pp)
  plot(den, main = "Fossil density")
  plot(pp, add = TRUE)
  plot(pp$window, add = TRUE)
  
  #Point of max/min erosion
  dr <- raster(den)
  max <- xyFromCell(dr, which.max(dr)) 
  min <- xyFromCell(dr, which.min(dr)) 

  
  maxp <- ppp(x = coords$x, y = coords$y, window = pp$window)
  plot(maxp, add = TRUE, cex = 2, pch = 8)
  minp <- ppp(x = min[1], y = min[2], window = pp$window)
  plot(minp, add = TRUE, col = 'red', cex = 2, pch = 8)
  legend("top", box.col = FALSE, legend = c("maximum density", "minimum density"), pch = 8, col = c("black", "red"))
  
  a <- coords$x
  b <- coords$y
  c <- min[1]
  d <- min[2]
  
  #null model
  csr <- ppm(pp ~ 1)
  hx <- ppm(pp ~ x)
  hy <- ppm(pp ~ y)
  hmin <- ppm(pp ~ sqrt((x-a)^2 +(y-b)^2))
  hmax <- ppm(pp ~ sqrt((x-c)^2 +(y-d)^2))
  #hline <- ppm(pp[[5]] ~ distance_line_to_point)
  
  
  #Assess model fit using spatial Kolmogorov-Smirnov tests 
  # - if p>0.05, the difference between two samples is not significant enough to say that they have different distribution
  #fmax <- function(x,y) {sqrt((x-a)^2 +(y-b)^2)}
  #fmin <- function(x,y) {sqrt((x-c)^2 +(y-d)^2)}

  
  x.test <- cdf.test(hx, test='ks')
  y.test <- cdf.test(hy, test='ks')
  max.test <- cdf.test(csr, test='ks')
  min.test <- cdf.test(hmax, test='ks')
  #line.test <- cdf.test(csr, distance_line_to_point, test = 'ks')

  lst <- list(x.test, y.test, max.test, min.test)
  
  #print results
  writeLines("Spatial Kolmogorov-Smirnov tests of model in two dimensions \n
    ------------------------------------------------------------ \n")
  print(paste("X: ", "p-value =", round(x.test$p.value, 4), ", D =", round(x.test$statistic, 4) ))
  print(paste("Y: ", "p-value =", round(y.test$p.value, 4), ", D =", round(y.test$statistic, 4) ))
  print(paste("MAX: ", "p-value =", round(max.test$p.value, 4), ", D =", round(max.test$statistic, 4) ))
  print(paste("MIN: ", "p-value =", round(min.test$p.value, 4), ", D =", round(min.test$statistic, 4) ))
  
 
  p_values <- unlist(lapply(lst, function(x){round(x$p.value, 4)}))
  teststat <- unlist(lapply(lst, function(x){round(x$statistic, 4)}))
  covariate <- c("x", "y", "max", "min")
  
  df <- cbind.data.frame(p_values, teststat, covariate)

  return(df)
}

#-------------------------------------------------------------------------------
# (4) AIC of heterogeneous Poisson models with spatial covaraites
#-------------------------------------------------------------------------------

#select options as appropriate

erosion_aic <- function(pp) {
  pp <- unmark(pp_mark)

  #Density map of everything on surface
  den <- density(pp)
  plot(den, main = "Fossil density")
  plot(pp, add = TRUE)
  plot(pp$window, add = TRUE)
  
  #Point of max/min erosion
  dr <- raster(den)
  max <- xyFromCell(dr, which.max(dr)) 
  min <- xyFromCell(dr, which.min(dr)) 
  
  maxp <- ppp(x = coords$x, y = coords$y, window = pp$window)
  plot(maxp, add = TRUE, cex = 2, pch = 8)
  minp <- ppp(x = min[1], y = min[2], window = pp$window)
  plot(minp, add = TRUE, col = 'red', cex = 2, pch = 8)
  legend("top", box.col = FALSE, legend = c("maximum density", "minimum density"), pch = 8, col = c("black", "red"))
  
  a <- max[1]
  b <- max[2]
  c <- min[1]
  d <- min[2]
  
  #null model
  csr <- ppm(pp ~ 1)
  
  #linear models
  hx <- ppm(pp ~ x)
  hy <- ppm(pp ~ y )
  hmin <- ppm(pp ~ sqrt((x-a)^2 +(y-b)^2) )
  hmax <- ppm(pp ~ sqrt((x-c)^2 +(y-d)^2) )
  #hline <- ppm(pp ~ distance_line_to_point))
  
  #quadratic models
  #hx <- ppm(pp ~ x + I(x^2))
  #hy <- ppm(pp ~ y + I(y^2))
  #hmin <- ppm(pp ~ sqrt((x-a)^2 +(y-b)^2) + I(sqrt((x-a)^2 +(y-b)^2)^2))
  #hmax <- ppm(pp ~ sqrt((x-c)^2 +(y-d)^2) + I(sqrt((x-c)^2 +(y-d)^2)^2))

  lst <- list(csr, hx, hy, hmin, hmax)
  aic <- lapply(lst, AIC)
  
  return(unlist(aic))
}

# EXAMPLE USAGE
#res <- lapply(pp_well, erosion_aic)
#res <- plyr::ldply(t(res))
#names(res) <- c("csr", "x", "y", "min", "max")
#write.csv(res, "./aic_linear_raw.csv")




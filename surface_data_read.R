## READING PRE-PARTITIONED DATA ##

#n.b. this data does not contain taxonomic information thus Fractofusus cannot be filtered out

source("./spatial_analysis_functions.R")
library(spatstat)
library(dplyr)

#-------------------------------------------------------------------------------
# (1) Read in surface data
#-------------------------------------------------------------------------------

#just points retro: eff/
#retro lengths: eff_retro/
#original: eff_raw/

filelist <- list.files(path = "./data2/", pattern = ".*.csv")

get_df <- function(filename){
  x <- (paste0("./data2/", filename))
  df <- read.table(x, sep = ",", header = T)
  return(df)
}

surface_data <- lapply(filelist, get_df) # List of dataframes

filenames <- gsub(".csv", "", filelist) # Vector of surface names

#-------------------------------------------------------------------------------
# (2) Make spatstat windows
#-------------------------------------------------------------------------------

## CREATE + EXPORT WINDOWS
# see also 'dex_outline.R'

create_windows <- function(j){
  plot(j$y ~j$x)
  win <- clickpoly(add = T)
  return(win)
}

win <- create_windows(surface_data[[5]])
coords <- data.frame(win$bdry[[1]])
write.csv(coords, "./win/E_win.csv") #creat folder to be populated

windows <- lapply(surface_data, create_windows)
coords <- lapply(windows, function(x){data.frame(x$bdry[[1]])})
windowpath <- unlist(lapply(filelist, function(x){paste("./win_raw/",x)}))
mapply(write.csv, coords, windowpath)

# READ IN WINDOWS .CSVs

window_files <- list.files(path = "./win_raw/", pattern = ".*.csv")

get_df2 <- function(filename){
  x <- (paste0("./win_raw/", filename))
  df <- read.table(x, sep = ",", header = F)
  df <- df[,2:3] # remove labels
  df <- df[-1,] # remove header
  return(df)
}

coords <- lapply(window_files, get_df2)

# MAKE WINDOWS 

windows <- lapply(coords, make_windows)

#-------------------------------------------------------------------------------
# (3) Make point patterns
#-------------------------------------------------------------------------------

#marked combined
pp_mark <- vector("list", length(surface_data))

for(i in length(well)){
  df <- surface_data
  pp_mark[[i]] <- ppp(df$x, df$y, marks = as.factor(df$group), windows[[i]])
  
}

#effaced fronds
eff <- lapply(surface_data, subset, group == "eff")

pp_eff <- vector("list", length(eff))

for(i in 1:length(eff)){
  pp_eff[[i]] <- ppp(eff[[i]]$x, eff[[i]]$y, window = windows[[i]]) 
  
}


#well-preserved fronds
well <- lapply(surface_data, subset, group == "well")

pp_well <- vector("list", length(well))

for(i in 1:length(well)){
  pp_well[[i]] <- ppp(well[[i]]$x, well[[i]]$y, window = windows[[i]])
}

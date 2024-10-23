## SURFACE POINT PATTERNS

source("./spatial_analysis_functions.R")
library(spatstat)
library(dplyr)


#-------------------------------------------------------------------------------
# (1) Read in surface data
#-------------------------------------------------------------------------------

#just points: data/
#retro: data_retro/
#original: data_raw/

filelist <- list.files(path = "./data_raw/", pattern = ".*.csv")

get_df <- function(filename){
  x <- (paste0("./data_raw/", filename))
  df <- read.table(x, sep = ",", header = T)
  return(df)
}

surface_data <- lapply(filelist, get_df) # List of dataframes

filenames <- gsub(".csv", "", filelist) # Vector of surface names

#-------------------------------------------------------------------------------
# (2) Subset effaced fronds
#-------------------------------------------------------------------------------

# Make taxon column names consistent

surface_data[[10]]$sp <- surface_data[[10]]$taxon # rename H5
surface_data[[10]] <- surface_data[[10]][,-3]

surface_data[[12]]$sp <- surface_data[[12]]$taxon # rename MUN
surface_data[[12]] <- surface_data[[12]][,-3]

# Make names consistent
surface_data[[10]]$sp <- gsub("frond", "unknown_frond", surface_data[[10]]$sp)
surface_data[[12]]$sp <- gsub("frond", "unknown_frond", surface_data[[12]]$sp)
     
#subset effaced fronds

subset_eff <- function(x){
  subset(x, sp != 'posseff' & (grepl('eff', sp, ignore.case = T)| grepl('poss', sp, ignore.case = T)| grepl('taph', sp, ignore.case = T))) # detect 'poss' or 'eff' in taxon names but remove 'posseff'
}

eff <- lapply(surface_data, subset_eff)

#table(unlist(lapply(eff, function(x){x$sp})))

# subset well-preserved fronds

source("./taxon_list.R") # get vector 'taxa'

well <- lapply(surface_data, filter, tolower(sp) %in% taxa)


#table(unlist(lapply(well, function(x){x$sp})))


#-------------------------------------------------------------------------------
# (3) Make spatstat windows
#-------------------------------------------------------------------------------

## CREATE + EXPORT WINDOWS
# see also 'dex_outline.R'

#create_windows <- function(j){
#  plot(j$y ~j$x)
#  win <- clickpoly(add = T)
#  return(win)
#}

#win <- create_windows(surface_data[[5]])
#coords <- data.frame(win$bdry[[1]])
#write.csv(coords, "./win/E_win.csv")

#windows <- lapply(surface_data, create_windows)
#coords <- lapply(windows, function(x){data.frame(x$bdry[[1]])})
#windowpath <- unlist(lapply(filelist, function(x){paste("./win_raw/",x)}))
#mapply(write.csv, coords, windowpath)

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
# (4) Make point patterns
#-------------------------------------------------------------------------------
point_patterns <- vector("list", length(surface_data))

for(i in 1:length(surface_data)){
  point_patterns[[i]] <- ppp(surface_data[[i]]$x, surface_data[[i]]$y, window = windows[[i]])
}
#marked combined
pp_mark <- vector("list", 14)

for(i in 1:13){
  df <- rbind(well[[i]], eff[[i]])
  df$group <- c(rep("well", length(well[[i]][,1])), rep("eff", length(eff[[i]][,1])))
  
  pp_mark[[i]] <- ppp(df$x, df$y, marks = as.factor(df$group), windows[[i]])
  
}

#effaced fronds
pp_eff <- vector("list", length(eff))

for(i in 1:length(eff)){
  pp_eff[[i]] <- ppp(eff[[i]]$x, eff[[i]]$y, window = windows[[i]]) 

}



#well-preserved fronds
#w/out fractos
#well <- lapply(well, subset, !grepl("fracto", sp, ignore.case = T))

pp_well <- vector("list", length(well))

for(i in 1:length(well)){
  pp_well[[i]] <- ppp(well[[i]]$x, well[[i]]$y, window = windows[[i]])
}


##plotting individual surface
df <- surface_data[[4]]
df$sp[df$sp != 'posseff' & (grepl('eff', df$sp, ignore.case = T)| grepl('poss', df$sp, ignore.case = T)| grepl('taph', df$sp, ignore.case = T))] <- "effaced_frond"
df$sp[! df$sp %in% taxa] <- "other"

mpD <- ppp(df$x, df$y, marks = df$sp, window = windows[[4]])
plot(mpD, cols = col, pch = 20 )
cols <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666", "#9E0142", "#D53E4F", "#F46D43", "#FDAE61", "#FEE08B", "#FFFFBF", "#E6F598", "#ABDDA4", "#66C2A5", "#3288BD", "#5E4FA2")
col<-c("#C6CDF7","#ECCBAE" , "#FF0000","#24281A","#E6A0C4","#7570B3","#046C9A", "#D69C4E", "#F3DF6C","#B40F20","#81A88D","#00A08A", "#5BBCD6", "#1B9E77", "#D95F02", "#35274A" ,"#E7298A", "#66A61E", "#D53E4F")
length(col)

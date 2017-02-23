setwd("F:\\R_Learning\\TestData\\CHM\\data")
getwd()

library(rgdal)
library(sp)
library(rgeos)
library(rgl)
library(lidR)
library(magrittr)
library(dtplyr)
library(dplyr)
library(raster)

LASfile <- system.file("data", "370000_5848000.las", package="lidR")
lidar15 = lidR::readLAS("370000_5848000.las")
summary(lidar15)
str(lidar15)
plot(lidar15)
lidar15@header
head(lidar15)
groundlidar15 <- lasfilter(lidar15, Classification == 2, ReturnNumber == 2)
groundlidar15
dtmgrd = grid_terrain(groundlidar15, method = "knnidw", k = 6)
plot(groundlidar15)
plot(dtmgrd)

#### -learning function ####
arth.mean <- function(x) {
  # arithmatic mean
  sum.x <- sum(x)
  n <- length(x)
  am <- sum.x/n
  return(am)
}
arth.mean(1:20)


#### -saving the output of Console ####
con <- file("test.log")
sink(con, append=TRUE)
sink(con, append=TRUE, type="message")

# This will echo all input and not truncate 150+ character lines...
source("script.R", echo=TRUE, max.deparse.length=10000)

# Restore output to console
sink() 
sink(type="message")

# And look at the log...
cat(readLines("test.log"), sep="\n")
####
## capture all the output to a file.
zz <- file("all.Rout", open = "wt")
sink(zz)
sink(zz, type = "message")
try(log("a"))
## revert output back to the console -- only then access the file!
sink(type = "message")
sink()
file.show("all.Rout")

sink(file ="F:\\R_Learning\\TestData\\CHM\\data\\save.txt")

sink()

install.packages("all.Rout")


loadhistory(file = ".Rhistory")

install.packages("jsonlite")
library(jsonlite)
plot3d(lidar15)
plot = rgl::plot3d(lidar15) #Error in as.double(x) : cannot coerce type 'S4' to vector of type 'double'

dtm15 = grid_terrain(lidar15)
plot(dtm15)
plot3d(dtm15)

canopy15 = grid_canopy(lidar15, res = 0.5) ## Using the local maximum algorithm, assigns the elevation of the highest 
                                           # return within each grid cell to the grid cell center
plot(canopy15)

lnorm = normalize(lidar15, dtm15)
plot(lnorm)

chm15 = grid_canopy(lnorm, res = 0.5, start = c(0, 0))
plot(chm15)

chm15 = as.raster(chm15)
plot(chm15)

# Canopy surface model: Local maximum algorithm with a resolution of 2 meters
lidar15 %>% grid_canopy(2) %>% plot

install.packages("arulesViz")
library(arulesViz)
install.packages("repos")
devtools::install_github("Jean-Romain/lidR", dependencies=TRUE)

install.packages("data.table")
install.packages("rgal")
install.packages("Rcpp")
install.packages("digest")

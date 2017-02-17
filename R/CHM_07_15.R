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
library(rLiDAR)
library(raster)

LASfile <- system.file("data", "370000_5848000.las", package="lidR")
lidar15 = lidR::readLAS("370000_5848000.las")
summary(lidar15)
plot(lidar15)
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




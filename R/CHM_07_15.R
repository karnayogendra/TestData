setwd("C:\\Users\\UOM74593User\\Documents\\R\\CHM\\data")
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
groundlidar15 <- lasfilter(lidar15, Classification == 2, ReturnNumber == 2)
groundlidar15
dtmgrd = grid_terrain(groundlidar15, method = "knnidw", k = 6)
plot(groundlidar15)
plot(dtmgrd)
str(dtmgrd)
summary(dtmgrd)
dtm15 <- raster("dtm_test.tiff")
plot3d(dtmgrd)

rdtmgrd = as.raster(dtmgrd)
plot(rdtmgrd)

writeRaster(rdtmgrd, filename = "dtm15.tif", format = "GTiff", overwirte = TRUE)

# SOURCE: http://neondataskills.org/R/Raster-Data-In-R/

#view coordinate reference system
dtm15@crs
# Set coordinate system
crs(r) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
## CRS arguments:
##  +proj=utm +zone=11 +datum=WGS84 +units=m +no_defs +ellps=WGS84
## +towgs84=0,0,0
# crs(dtm15) <- "+proj=MGA +zone=55 +datum=D_GDA_1994 +units=m +no_defs +Sph=GRS_1980 +toD_GDA_1994=0,0,0"
crs(dtm15) <- "+proj=utm +zone=55 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"

#view raster extent
dtm15@extent

#plot the raster
#note that this raster represents a small region of the study site
plot(dtm15, main="Digital Terrain Model, My Study Site")



#Get min and max cell values from raster
#NOTE: this code may fail if the raster is too large
cellStats(dtm15, min)

cellStats(dtm15, max)

cellStats(dtm15, range)

#we can look at the distribution of values in the raster too
hist(dtm15, main="Distribution of elevation values", col= "purple", maxpixels=21752940)

#### -Basic Raster Math ####

#multiple each pixel in the raster by 2
DTM2 <- dtm15 * 2
DTM2

#plot the new DEM
plot(DTM2, main="DTM with all values doubled")


#create a plot of our raster: Plotting Raster Data
image(dtm15)

#specify the range of values that you want to plot in the DEM
#just plot pixels between 250 and 300 m in elevation
image(dtm15, zlim=c(700,800))

#we can specify the colors too
col <- terrain.colors(5)
image(dtm15, zlim=c(700,800), main="Digital Terrain Model (DTM)", col=col)

# Breaks and Colorbars in R

#add a color map with 5 colors
col=terrain.colors(5)
#add breaks to the colormap (6 breaks = 5 segments)
brk <- c(680, 700, 720, 740, 760, 780)
plot(dtm15, col=col, breaks=brk, main="DTM with more breaks")

# Expand right side of clipping rect to make room for the legend
par(xpd = FALSE,mar=c(5.1, 4.1, 4.1, 4.5))
#DEM with a custom legend
plot(dtm15, col=col, breaks=brk, main="DTM with a Custom (buf flipped) Legend",legend = FALSE)

#turn xpd back on to force the legend to fit next to the plot.
par(xpd = TRUE)
#add a legend - but make it appear outside of the plot
legend(par()$usr[2], 5848500,legend = c("lowest", "a bit higher", "middle ground", "higher yet", "Highest"), fill = col)

#LEGEND in REVERSE oRDER

# Expand right side of clipping rect to make room for the legend
par(xpd = FALSE,mar=c(5.1, 4.1, 4.1, 4.5))
#DEM with a custom legend
plot(dtm15, col=col, breaks=brk, main="DTM with a Custom Fixed Legend",legend = FALSE)

#turn xpd back on to force the legend to fit next to the plot.
par(xpd = TRUE)
#add a legend - but make it appear outside of the plot
legend(par()$usr[2], 5848500,legend = c("Highest", "Higher yet", "Middle", "A bit higher", "Lowest"), fill = col)

#WITH FEWER BREAKS
#add a color map with 4 colors
col=terrain.colors(4)
#add breaks to the colormap (6 breaks = 5 segments)
brk <- c(700, 720, 740, 760, 780)
plot(dtm15, col=col, breaks=brk, main="DTM with fewer breaks")

## HOW TO MAKE LEGEND with unique categories?????

#### -Cropping Rasters in R ####

plot(dtm15)
#Define the extent of the crop by clicking on the plot
cropbox1 <- drawExtent()
#crop the raster, then plot the new cropped raster
dtm15crop1 <- crop(dtm15, cropbox1)

#plot the cropped extent
plot(dtm15crop1) #You can also manually assign the extent coordinates to be used to crop a raster.  We'll need the extent defined as (`xmin`, `xmax`, `ymin` , `ymax`) to do this.  This is how we'd crop using a GIS shapefile (with a rectangular shape)


#define the crop extent
cropbox2 <-c(370500,370600,5848070, 5848500)
#crop the raster
DTMcrop2 <- crop(dtm15, cropbox2)
#plot cropped DEM
plot(DTMcrop2)

#### -Write Raster #### 
rDTMcrop2 = as.raster(DTMcrop2)
plot(rDTMcrop2)

writeRaster(rdtmgrd, filename = "rDTMcrop2.tif", format = "GTiff", overwirte = TRUE)


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

lnorm = lasnormalize(lidar15, dtm15) 
###Warning message:Dataset may be invalid: 549620 points below 0 found. 
plot(lnorm)

chm15 = grid_canopy(lnorm, res = 0.5, start = c(0, 0))
plot(chm15)

rchm15 = raster::as.raster(chm15)
plot(rchm15)
writeRaster(rchm15, filename = "CHM15_M.tif", format = "GTiff", overwirte = TRUE)

writeRaster(rchm15, filename = "~/R/CHM/output/CHM15_M.tif", format = "GTiff", overwirte = TRUE)

dir.create("C:\\Users\\UOM74593User\\Documents\\R\\CHM\\data\\Raster")

#### -HOW TO WRITE new folder and file name directly into the existing directory i.e. getwd ####

#### -DSM ####
DSM15 = grid_metrics(lidar15, max (Z), 2)
plot(DSM15)
str(DSM15)
#Writing DSM as Raster
rDSM15 = as.raster(DSM15)
writeRaster(rDSM15, filename = "~/R/CHM/output/DSM15.tif", format = "GTiff", overwrite = TRUE)

#### -CHM from subtract ####
CHMNew <- subtract(DSM15, dtm15)
CHMNew
plot(CHMNew)

## CHM can be obtained from subtracting DTM from DSM directly from ArcMAP using Minus (3D analyst tools)
CHMNew <- DSM15 - dtm15
CHMNew
plot(CHMNew)
dev.off()
### DOES not work: strange plot

if(na.fill != "none"){
  ex = extent(.las)
  grid = make_grid(ex@xmin, ex@xmax, ex@ymin, ex@ymax, res)
  
  data.table::setkeyv(grid, c("X", "Y"))
  data.table::setkeyv(dsm, c("X", "Y"))
  data.table::setattr(dsm, "class", class(grid))
  
  dsm = dsm[grid]
  
  z = interpolate(dsm[!is.na(Z)], dsm[is.na(Z)], method = na.fill, ...)
  
  dsm[is.na(Z), Z := z]
  
  as.lasmetrics(dsm, res)
}

  rm(.las)
  gc()

  return(dsm) 
}

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

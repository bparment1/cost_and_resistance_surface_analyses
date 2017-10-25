############### SESYNC Research Support: Fisheries and food security ########## 
## Importing and processing data from survey for the fisheries project at SESYNC.
## 
## DATE CREATED: 06/06/2017
## DATE MODIFIED: 10/12/2017
## AUTHORS: Benoit Parmentier 
## PROJECT: Garden Wealth (urban garden)
## ISSUE: 
## TO DO:
##
## COMMIT: initial commit gDistance example
##
## Links to investigate:

###################################################
#

###### Library used

library(gtools)                              # loading some useful tools 
library(sp)                                  # Spatial pacakge with class definition by Bivand et al.
library(spdep)                               # Spatial pacakge with methods and spatial stat. by Bivand et al.
library(rgdal)                               # GDAL wrapper for R, spatial utilities
library(raster)
library(gdata)                               # various tools with xls reading, cbindX
library(rasterVis)                           # Raster plotting functions
library(parallel)                            # Parallelization of processes with multiple cores
library(maptools)                            # Tools and functions for sp and other spatial objects e.g. spCbind
library(maps)                                # Tools and data for spatial/geographic objects
library(plyr)                                # Various tools including rbind.fill
library(spgwr)                               # GWR method
library(rgeos)                               # Geometric, topologic library of functions
library(gridExtra)                           # Combining lattice plots
library(colorRamps)                          # Palette/color ramps for symbology
library(ggplot2)
library(lubridate)
library(dplyr)
library(rowr)                                # Contains cbind.fill
library(car)
library(sf)
library(gdistance)

###### Functions used in this script and sourced from other files

create_dir_fun <- function(outDir,out_suffix=NULL){
  #if out_suffix is not null then append out_suffix string
  if(!is.null(out_suffix)){
    out_name <- paste("output_",out_suffix,sep="")
    outDir <- file.path(outDir,out_name)
  }
  #create if does not exists
  if(!file.exists(outDir)){
    dir.create(outDir)
  }
  return(outDir)
}

#Used to load RData object saved within the functions produced.
load_obj <- function(f){
  env <- new.env()
  nm <- load(f, env)[1]
  env[[nm]]
}

### Other functions ####

#function_processing_data <- ".R" #PARAM 1
#script_path <- "/nfs/bparmentier-data/Data/projects/urban_garden_pursuit/scripts" #path to script #PARAM 
#source(file.path(script_path,function_processing_data)) #source all functions used in this script 1.

############################################################################
#####  Parameters and argument set up ###########

in_dir <- "/nfs/bparmentier-data/Data/projects/urban_garden_pursuit/data" #local bpy50 , param 1
out_dir <- "/nfs/bparmentier-data/Data/projects/urban_garden_pursuit/outputs" #param 2

num_cores <- 2 #param 8
create_out_dir_param=TRUE # param 9

out_suffix <-"connectivity_example_10122017" #output suffix for the files and ouptut folder #param 12

############## START SCRIPT ############################

######### PART 0: Set up the output dir ################

if(is.null(out_dir)){
  out_dir <- in_dir #output will be created in the input dir
  
}
#out_dir <- in_dir #output will be created in the input dir

out_suffix_s <- out_suffix #can modify name of output suffix
if(create_out_dir_param==TRUE){
  out_dir <- create_dir_fun(out_dir,out_suffix)
  setwd(out_dir)
}else{
  setwd(out_dir) #use previoulsy defined directory
}

### PART I READ AND PREPARE DATA #######
#set up the working directory
#Create output directory

lf_dir <- list.files(in_dir,full.names=T) #this is the list of folder with RAW data information

set.seed(123)
r <- raster(ncol=3,nrow=3)
r[] <- 1:ncell(r)
r
plot(r)

r[] <- 1
tr1 <- transition(r,transitionFunction = mean,directions=8)

tr1
plot(tr1)

r2 <- r
r2[] <- runif(9)
ncf <- function(x){max(x) - x[1] + x[2]}
tr2 <- transition(r2,ncf,4,symm=FALSE)

transitionMatrix(tr2)
tr2
image(transitionMatrix(tr2))

tr1 # with dsCMatrix object (symmetric)
tr2 # with dgCMatrix object (asymmetric)

tr3 <- tr1*tr2
tr3 <- tr1+tr2
tr3 <- tr1*3
tr3 <- sqrt(tr1)

image(transitionMatrix(tr3))

plot(raster(tr3), main="raster(tr3)", xlab="Longitude (degrees)",
     ylab="Latitude (degrees)")

#This paper describe the links between conductance, friction etc.
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2920323/

tr1C <- geoCorrection(tr1, type="c")
tr2C <- geoCorrection(tr2, type="c")

### Cost distance and least cost path
r3 <- raster(ncol=18,nrow=9)
r3 <- setValues(r3,runif(18*9)+5)
plot(r3)

tr3 <- transition(r3, mean, 4)
tr3C <- geoCorrection(tr3, type="c", multpl=FALSE, scl=TRUE)
tr3R <- geoCorrection(tr3, type="r", multpl=FALSE, scl=TRUE)

CorrMatrix <- geoCorrection(tr3, type="r", multpl=TRUE, scl=TRUE)
tr3R <- tr3 * CorrMatrix

plot(raster(tr3R))
plot(raster(tr3C))
plot(raster(tr3))

###### Let's make a mock example:

sP <- cbind(c(-100,-100,100),c(50,-50,50))

points <- SpatialPoints(sP)
  
plot(raster(tr3C))
plot(points,add=T)
text(points,c(1,2,3))

costDistance(tr3C, sP)
commuteDistance(tr3R, sP)

rSPDistance(tr3R, sP, sP, theta=1e-12, totalNet="total")

##### PASSAGE

origin <- SpatialPoints(cbind(0, 0))
rSPraster <- passage(tr3C, origin, sP[1,], theta=3)

rSPraster
plot(rSPraster)
text(points[1,],"sP1")
text(origin,"O")

### Non overlapping trajectories

r1 <- passage(tr3C, origin, sP[1,], theta=1)
r2 <- passage(tr3C, origin, sP[2,], theta=1)
rJoint <- min(r1, r2)
rDiv <- max(max(r1, r2) * (1 - min(r1, r2)) - min(r1, r2), 0)

#### Hiking example

r <- raster(system.file("external/maungawhau.grd", package="gdistance"))
plot(r)


#The Hiking Function requires the slope (m) as input, which can be calculated from the altitude
#(z) and distance between cell centres (d).
#mij = (zj − zi)/dij
#The units of altitude and distance should be identical. Here, we use meters for both. First, we
#calculate the altitudinal differences between cells. Then we use the geoCorrection function
#to divide by the distance between cells.

altDiff <- function(x){x[2] - x[1]}
hd <- transition(r, altDiff, 8, symm=FALSE)
slope <- geoCorrection(hd)

plot(raster(slope))

adj <- adjacent(r, cells=1:ncell(r), pairs=TRUE, directions=8) 
speed <- slope
speed[adj] <- 6 * exp(-3.5 * abs(slope[adj] + 0.05)) #Tobbler Hiking function
Conductance <- geoCorrection(speed)

plot(raster(Conductance))

A <- c(2667670, 6479000)
B <- c(2667800, 6479400)
AtoB <- shortestPath(Conductance, A, B, output="SpatialLines")
BtoA <- shortestPath(Conductance, B, A, output="SpatialLines")

plot(r, xlab="x coordinate (m)", ylab="y coordinate (m)",legend.lab="Altitude (masl)")
lines(AtoB, col="red", lwd=2)
lines(BtoA, col="blue")
text(A[1] - 10, A[2] - 10, "A")
text(B[1] + 10, B[2] + 10, "B")


###################### END OF SCRIPT ################
############### SESYNC Research Support: Fisheries and food security ########## 
## Importing and processing data from survey for the fisheries project at SESYNC.
## 
## DATE CREATED: 06/06/2017
## DATE MODIFIED: 04/24/2018
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
library(rgrass7)

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

out_suffix <- "connectivity_example_04242018" #output suffix for the files and ouptut folder #param 12

in_dir <- "/nfs/bparmentier-data/Data/projects/urban_garden_pursuit/data_urban_garden"
out_dir <- "/nfs/bparmentier-data/Data/projects/urban_garden_pursuit/outputs"
in_dir_grass <- "/nfs/bparmentier-data/Data/projects/urban_garden_pursuit/data_urban_garden" 
#in_dir_grass <- "/nfs/bparmentier-data/Data/projects/urban_garden_pursuit/grass_data_urban_garden"

# background reading:
# https://grass.osgeo.org/grass72/manuals/grass_database.html

gisBase <- '/usr/lib/grass72'
#gisDbase <- '/nfs/urbangi-data/grassdata'
gisDbase <- in_dir_grass #should be the same as in_dir

#location <- 'DEM_LiDAR_1ft_2010_Improved_NYC_int'
location <- 'NYC_example'
location <- 'connectivy_example'

file_format <- ".tif" #PARAM5
NA_flag_val <- -9999 #PARAM7
out_suffix <-"ny_example_04232018" #output suffix for the files and ouptu folder #PARAM 8
create_out_dir_param=TRUE #PARAM9

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

mapset <- "nyc_site_test"
# initialize a mapset for watershed estimation results
initGRASS(gisBase = gisBase, #application location
          gisDbase = gisDbase,  #database dir
          location = location, #grass location
          mapset = mapset, # grass mapset
          override = TRUE
)

### PART I READ AND PREPARE DATA #######
#set up the working directory
#Create output directory

###################### PART 2: compare with GRASS for random walk ###########

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

##### Generate Nodes:

### Add a point
A <- c(2667670, 6479000)
B <- c(2667800, 6479400)
C <- c(2667899,6478800)

net2_sp <- SpatialPoints(rbind(A,B,C))
plot(r, xlab="x coordinate (m)", ylab="y coordinate (m)",legend.lab="Altitude (masl)")
plot(net2_sp,add=T)
test <- shortestPath(Conductance, net2_sp, net2_sp, output="SpatialLines")
class(test)

plot(r)
plot(test,add=T)

dist_test <- distance(r,net2_sp)
dist_test <- distanceFromPoints(r,net2_sp)

##### Shortest path
AtoB <- shortestPath(Conductance, A, B, output="SpatialLines")
BtoA <- shortestPath(Conductance, B, A, output="SpatialLines")
#Add new path/route
BtoC <- shortestPath(Conductance, B, C, output="SpatialLines")
CtoB <- shortestPath(Conductance, B, C, output="SpatialLines")
#Add new path/route
AtoC <- shortestPath(Conductance, A, C, output="SpatialLines")
CtoA <- shortestPath(Conductance, C, A, output="SpatialLines")

##### Random walk: commute distance

altDiff <- function(x){x[2] - x[1]}
hd <- transition(r, altDiff, 8, symm=FALSE)
#Create a Transition object from the raster
tr <- transition(r,function(x) 1/mean(x),8)

test_path <- commuteDistance(tr,net2_sp) 
plot(net2_sp)
net2_sf <- as(net2_sp,"sf")
plot(net2_sf$geometry)
#View(net2_sf)
net2_sf$ID <- 1:nrow(net2_sf)

st_write(net2_sf,"network_nodes.shp",delete_layer = T)
writeRaster(r,"r_surf.tif")

#### Add GRASS code here:

execGRASS("v.in.ogr",flags = c("o","overwrite"), 
          input="network_nodes.shp", 
          output="nodes_origin")

execGRASS("r.in.gdal",flags=c("o","overwrite"), 
          input="r_surf.tif", 
          output="r_surf")
system("r.info r_surf")
r #check we have the same res, etc.

#### Set region extent and resolution first
system("g.region -p") #Exaine current region properties
#system("g.region -p") #Exaine current region properties

system("g.region rast=r_surf")
system("g.region -p")

system("v.to.rast --overwrite input=nodes_origin use=attr output=nodes_origin_surf attribute_column=ID")
system("r.info nodes_origin_surf")

system("r.mapcalc 'r_friction = 1'") #creates a raster with value 1
system("r.info r_friction")

# compute cumulative cost surfaces
system("r.walk -k elev=r_surf friction=r_friction output=walk.cost start_points=nodes_origin stop_points=nodes_origin lambda=1")


#execGRASS("r.randomwalk",flags=c("o","overwrite"), 
#          elevation="r_surf", 
#          releasemp="nodes_origin_surf",
#          output="r_surf")


#r.randomwalk help
#r.randomwalk [-abkmnpqsvx] prefix=string [cores=integer] [cellsize=float]
#[aoicoords=float,...][aoimap=name] elevation=name [releasefile=string] 
#[caserules=integer,integer,...] [releasemap=name] [depositmap=name]
#[impactmap=name] [probmap=name] [scoremap=name] [impactobjects=name] 
#[objectscores=string] models=string mparams=string [sampling=integer]
#[retain=float] [functype=integer] [backfile=string] [cdffile=string]
#[zonalfile=string] [profile=float,...] [--verbose] [--quiet]

###################### END OF SCRIPT ################
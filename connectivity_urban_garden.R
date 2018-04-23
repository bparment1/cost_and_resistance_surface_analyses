####################################   Connectivity and cost Analysis   #######################################
#########################################   Urban garden Pursuit   #######################################
# This script explores shortest path with real data from the Urban garden pursuit group.
#
#
#AUTHORS: Benoit Parmentier                                             
#DATE CREATED: 04/16/2018 
#DATE MODIFIED: 04/23/2018
#Version: 2
#PROJECT: SESYNC and AAG 2018 workshop/Short Course preparation
#TO DO:
#
#COMMIT: resizing image for test and subsetting nodes
#
#################################################################################################

###Loading R library and packages                                                      

library(gstat) #spatial interpolation and kriging methods
library(sp) # spatial/geographfic objects and functions
library(rgdal) #GDAL/OGR binding for R with functionalities
library(spdep) #spatial analyses operations, functions etc.
library(gtools) # contains mixsort and other useful functions
library(maptools) # tools to manipulate spatial data
library(parallel) # parallel computation, part of base package no
library(rasterVis) # raster visualization operations
library(raster) # raster functionalities
library(forecast) #ARIMA forecasting
library(xts) #extension for time series object and analyses
library(zoo) # time series object and analysis
library(lubridate) # dates functionality
library(colorRamps) #contains matlab.like color palette
library(rgeos) #contains topological operations
library(sphet) #contains spreg, spatial regression modeling
library(BMS) #contains hex2bin and bin2hex, Bayesian methods
library(bitops) # function for bitwise operations
library(foreign) # import datasets from SAS, spss, stata and other sources
#library(gdata) #read xls, dbf etc., not recently updated but useful
library(classInt) #methods to generate class limits
library(plyr) #data wrangling: various operations for splitting, combining data
library(readxl) #functionalities to read in excel type data
library(gdistance) #package
library(rgrass7)
library(sf)

###### Functions used in this script

#function_preprocessing_and_analyses <- "fire_alaska_analyses_preprocessing_functions_03102017.R" #PARAM 1
script_path <- "/nfs/bparmentier-data/Data/projects/urban_garden_pursuit/scripts"
#source(file.path(script_path,function_preprocessing_and_analyses)) #source all functions used in this script 1.

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

#####  Parameters and argument set up ###########

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
location <- 'nyc_site2'

origin_fname <-  "OrigNYCSurface.tif"
biosurf_fname <- "BioSurfaceFinal.tif"
origin_garden_node_fname<- "OrigGardenNodes.tif"
new_node_fname <- "NewNodes.tif"

file_format <- ".tif" #PARAM5
NA_flag_val <- -9999 #PARAM7
out_suffix <-"ny_example_04232018" #output suffix for the files and ouptu folder #PARAM 8
create_out_dir_param=TRUE #PARAM9

################# START SCRIPT ###############################

## First create an output directory

if(is.null(out_dir)){
  out_dir <- dirname(in_dir) #output will be created in the input dir
}

out_suffix_s <- out_suffix #can modify name of output suffix
if(create_out_dir_param==TRUE){
  out_dir <- create_dir_fun(out_dir,out_suffix_s)
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

#g.gisenv set="LOCATION_NAME=newlocation"


####  PART I: EXPLORE DATA READ AND DISPLAY INPUTS #######


r_origin <- raster(file.path(in_dir,origin_fname)) #"OrigNYCSurface.tif"
r_bio <- raster(file.path(in_dir,biosurf_fname)) #<- "BioSurfaceFinal.tif"
r_origin_garden <- raster(file.path(in_dir,origin_garden_node_fname)) #<- "OrigGardenNodes.tif"
r_new_node <- raster(file.path(in_dir,new_node_fname)) # <- "NewNodes.tif"

## 1) Check if we can use point for origin?

r_origin
tr1_origin <- transition(r_origin,transitionFunction = mean,directions=8)

#This generates the following error:

#> tr1_origin <- transition(r_origin,transitionFunction = mean,directions=8)
#Error in validObject(.Object) : 
#  invalid class “dsCMatrix” object: Negative value in Dim
#In addition: Warning message:
#  In .nextMethod(.Object = .Object, ... = ...) :
#  NAs introduced by coercion to integer range

freq_origin_tb<- as.data.frame(freq(r_origin_garden))
freq_new_node_tb<- as.data.frame(freq(r_new_node))

#class(freq_origin_tb)

write.table(freq_origin_tb,"freq_origin_tb.txt",sep=",")
write.table(freq_new_node_tb,"new_node_tb.txt",sep=",")

##### Convert Nodes to polygons? This is a very quick option using gdal_polygonize.py
setwd(in_dir)
#cmd_str <- "gdal_polygonize.py input.asc -f 'GeoJSON' output.json"
cmd_str <- "gdal_polygonize.py OrigGardenNodes.tif -f 'ESRI SHapefile' OrigGardenNodes.shp"

system(cmd_str)

orig_nodes_sf <- st_read("OrigGardenNodes.shp")

plot(orig_nodes_sf)
dim(orig_nodes_sf)
table(orig_nodes_sf$DN)
#View(orig_nodes_sf)

new_node_fname <- "NewNodes.tif"
cmd_str <- "gdal_polygonize.py NewNodes.tif -f 'ESRI SHapefile' NewNodes.shp"
system(cmd_str)

new_nodes_sf <- st_read("NewNodes.shp")
plot(new_nodes_sf,border="red",add=T)
dim(new_nodes_sf)
#View(new_node_fname)
dim(freq_new_node_tb)
View(new_nodes_sf)
table(new_nodes_sf$DN)

centroids_new_nodes_sf <- st_centroid(new_nodes_sf)
centroids_orig_nodes_sf <- st_centroid(orig_nodes_sf)
plot(centroids_orig_nodes_sf)
dim(centroids_orig_nodes_sf)

##Figure nodes and suitability surface
res_pix<-960
col_mfrow<-1
row_mfrow<-1
png(filename=paste("Figure_NYC_site_nodes",out_suffix,".png",sep=""),
    width=col_mfrow*res_pix,height=row_mfrow*res_pix)

plot(r_bio)
plot(centroids_orig_nodes_sf,col="blue",add=T,cex=0.5)
plot(centroids_new_nodes_sf,col="red",add=T,cex=0.5)

dev.off()

st_write(centroids_new_nodes_sf,"centroids_new_nodes.shp",delete_dsn = T) #overwrite using delete_dsn
st_write(centroids_orig_nodes_sf,"centroids_orig_nodes.shp",delete_dsn = T)


#### Let's aggregate?


#https://casoilresource.lawr.ucdavis.edu/software/grass-gis-raster-vector-and-imagery-analysis/raster-operations/simple-comparision-two-least-cost-path-approaches/

#https://gis.stackexchange.com/questions/112304/how-to-create-least-cost-path-between-two-polygons-with-grass
#https://stackoverflow.com/questions/9605827/least-cost-path-with-multiple-points
#https://grasswiki.osgeo.org/wiki/Working_with_external_data_in_GRASS_7
#http://ncsu-geoforall-lab.github.io/geospatial-modeling-course/grass/buffers_cost.html

# register (rather than import) a GeoTIFF file in GRASS GIS:
#r.external input=terra_lst1km20030314.LST_Day.tif output=modis_celsius
#r_bio <- raster(file.path(in_dir,biosurf_fname)) #<- "BioSurfaceFinal.tif"

#cmd_str <- paste("r.external")
#cmd_str <- "r.external input="
#system()
# define output directory for files resulting from subsequent calculations:
#r.external.out directory=$HOME/gisoutput/ format="GTiff"

# perform calculations (here: extract pixels > 20 deg C)
# store output directly as GeoTIFF file, hence add the .tif extension:
#r.mapcalc "warm.tif = if(modis_celsius > 20.0, modis_celsius, null() )"

# cease GDAL output connection and turn back to write standard GRASS raster files:
#r.external.out -r

# use the result elsewhere
#qgis $HOME/gisoutput/warm.tif

#v.in.ogr input=/home/user/shape_data/test_shape.shp output=grass_map 


#exec("r.in.gdal" input=E:\cdnh43e_v1.1r1.tif output=cdnh43e_v1 location=LCC

execGRASS("r.in.gdal",flags=c("o","overwrite"), input=biosurf_fname, output="biosurf")
#execGRASS("r.in.gdal",flags=c("o","overwrite"), input=origin_fname, output="biosurf")

projection(r_bio)==st_crs(centroids_new_nodes_sf)$proj4string
#note projections not defined the same way so use flag -o to ignore
execGRASS("v.in.ogr",flags=c("o","overwrite"), input="centroids_new_nodes.shp", output="centroids_new_nodes")
          #,location="nyc_site")
execGRASS("v.in.ogr",flags = c("o","overwrite"), 
          input="centroids_orig_nodes.shp", output="centroids_orig_nodes")
#,location="nyc_site")

#http://gracilis.carleton.ca/CUOSGwiki/index.php/Evaluating_Landscape_Permeability_in_Quantum)

#execGRASS("d.rast" ,map="biosurf")



r.cost_param <- list(
                     input="biosurf", 
                     output="biosurf_cost",
                     outdir="biosurf_direction",
                     start_points="centroids_orig_nodes",
                     stop_points="centroids_new_nodes"
)

execGRASS('r.cost',flags=c("k","overwrite"), parameters = r.cost_param)

execGRASS("r.cost", flags=c("k","overwrite"),input="biosurf", output="biosurf_cost",
          outdir="biosurf_direction",
          start_points="centroids_orig_nodes",
          stop_points="centroids_new_nodes" )

execGRASS("r.cost", flags=c("k","overwrite"),input="biosurf", output="biosurf_cost",
          outdir="biosurf_direction",
          start_points="centroids_orig_nodes@nyc_site_test",
          stop_points="centroids_new_nodes@nyc_site_test" )

#### describe columns of vector database
system("db.describe table=centroids_orig_nodes")
#db.describe -c table=vect_map
system("db.describe -c table=centroids_orig_nodes") # just name of columns

system("v.info centroids_orig_nodes")
system("r.info biosurf")

#### Set region extent and resolution first
system("g.region -p") #Exaine current region properties

system("g.region rast=biosurf")
system("g.region -p")
#system("v.to.rast input=centroids_orig_nodes output=centroids_orig_nodes_surf attribute_column=DN")
system("v.to.rast --overwrite input=centroids_orig_nodes use=attr output=centroids_orig_nodes_surf attribute_column=DN")
system("r.info centroids_orig_nodes_surf")

#"rast=name[,name,...]
#Set region to match this raster map")
#centroids_orig_nodes_surf

execGRASS("r.cost", flags=c("k","overwrite"),input="biosurf", output="biosurf_cost",
          outdir="biosurf_direction",
          start_raster="centroids_orig_nodes_surf")

#"centroids_new_nodes.shp"
 #execGRASS("r.drain", input="biosurf_cost", 
#          output="biosurf_drain",vector_points="centroids_new_nodes")

## 
execGRASS("r.mapcalc", "test = 'biosurf == 1'")
system("r.mapcalc 'biosurf == 1'")
#r.mapcalc "friction = 1.0"

execGRASS("r.walk", input="biosurf_cost", 
          output="biosurf_drain",vector_points="centroids_new_nodes")


# compute shortest path from start to end points
r.drain in=walk.cost out=walk.drain vector_points=end
r.drain in=cost out=cost.drain vector_points=end

# compute cumulative cost surfaces
#r.walk -k elev=elev friction=friction out=walk.cost start_points=start stop_points=end lambda=1

#r.cost -k in=slope out=cost start_points=start stop_points=end

# generate cost "friction" maps
# unity friction for r.walk
r.mapcalc "friction = 1.0"
# use slope percent as the friction for r.cost
r.slope.aspect elev=elev slope=slope format=percent

# compute cumulative cost surfaces
r.walk -k elev=elev friction=friction out=walk.cost start_points=start stop_points=end lambda=1

# Peak cost value: 14030.894105
r.cost -k in=slope out=cost start_points=start stop_points=end
# Peak cost value: 11557.934493

# compute shortest path from start to end points
r.drain in=walk.cost out=walk.drain vector_points=end
r.drain in=cost out=cost.drain vector_points=end

######################## End of Script ###########################
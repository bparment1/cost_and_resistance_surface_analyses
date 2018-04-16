####################################   Connectivity and cost Analysis   #######################################
#########################################   Urban garden Pursuit   #######################################
# This script explores shortest path with real data from the Urban garden pursuit group.
#
#Goal: Determine the ten (10) parcels of land within Clay County in the focus zone most suitable for purchase
#towards conversion to land conservation.
#
#AUTHORS: Benoit Parmentier                                             
#DATE CREATED: 04/16/2018 
#DATE MODIFIED: 04/16/2018
#Version: 2
#PROJECT: SESYNC and AAG 2018 workshop/Short Course preparation
#TO DO:
#
#COMMIT: testing code for urban garden
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
library(gdistance)

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

origin_fname <-  "OrigNYCSurface.tif"
biosurf_fname <- "BioSurfaceFinal.tif"
origin_garden_node_fname<- "OrigGardenNodes.tif"
new_node_fname <- "NewNodes.tif"

file_format <- ".tif" #PARAM5
NA_flag_val <- -9999 #PARAM7
out_suffix <-"ny_example_04162018" #output suffix for the files and ouptu folder #PARAM 8
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

####  PART I: EXPLORE DATA READ AND DISPLAY INPUTS #######


r_origin <- raster(file.path(in_dir,origin_fname)) #"OrigNYCSurface.tif"
r_bio <- raster(file.path(in_dir,biosurf_fname)) #<- "BioSurfaceFinal.tif"
r_origin_garden <- raster(file.path(in_dir,origin_garden_node_fname)) #<- "OrigGardenNodes.tif"
r_new_node <- raster(file.path(in_dir,new_node_fname)) # <- "NewNodes.tif"

stack(r_origin,r_bio,r_origin_garden,r_new_node)

plot(r_origin)



######################## End of Script ###########################
##!/bin/bash

grass /nfs/bparmentier-data/Data/projects/urban_garden_pursuit/data_urban_garden
LOCATION_NAME = connectivy_example

LOCATION = /nfs/bparmentier-data/Data/projects/urban_garden_pursuit/data_urban_garden/connectivy_example/PERMANENT
GISDBASE = /nfs/bparmentier-data/Data/projects/urban_garden_pursuit/data_urban_garden
LOCATION_NAME = connectivy_example
MAPSET = PERMANENT

LOCATION = /usr/local/share/grassdata/spearfish70/PERMANENT
GISDBASE = /usr/local/share/grassdata
LOCATION_NAME = spearfish70
MAPSET = PERMANENT

in_dir_grass <- "/nfs/bparmentier-data/Data/projects/urban_garden_pursuit/data_urban_garden" 


gisBase <- '/usr/lib/grass72'
#gisDbase <- '/nfs/urbangi-data/grassdata'
gisDbase <- in_dir_grass #should be the same as in_dir

#location <- 'DEM_LiDAR_1ft_2010_Improved_NYC_int'
location <- 'NYC_example'
location <- 'connectivy_example'

grass -text ~/grassdata/mylocation/mymapset

#### Start here:
## Type this to start grass from permanent
grass /nfs/bparmentier-data/Data/projects/urban_garden_pursuit/data_urban_garden/connectivy_example/PERMANENT
grass /nfs/bparmentier-data/Data/projects/urban_garden_pursuit/data_urban_garden/connectivy_example/nyc_site_test

#### Add GRASS code here:

v.in.ogr" -f "o" input="network_nodes.shp" output="nodes_origin"

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
#system("r.walk -k elev=r_surf friction=r_friction output=walk.cost start_points=nodes_origin stop_points=nodes_origin lambda=1")

execGRASS("r.cost", flags=c("k","overwrite"),
          input="r_surf", 
          output="r_surf_cost",
          outdir="r_surf_direction",
          start_raster="nodes_origin_surf")

# compute shortest path from start to end points
#execGRASS()
#system("r.drain in=walk.cost out=walk.drain vector_points=end")
system("r.drain input=walk.cost output=walk.drain vector_points=nodes_origin")
system("r.drain input=r_surf_cost output=cost_drain vector_points=nodes_origin")

system("r.out.gdal input=walk.drain output=path_r_walk.tif") 
system("r.out.gdal input=walk.cost output=walk_cost.tif") 
system("r.out.gdal input=cost_drain output=cost_drain.tif") 


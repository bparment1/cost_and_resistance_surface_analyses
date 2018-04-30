##!/bin/bash

grass /nfs/bparmentier-data/Data/projects/urban_garden_pursuit/data_urban_garden
LOCATION_NAME = connectivy_example

LOCATION = /nfs/bparmentier-data/Data/projects/urban_garden_pursuit/data_urban_garden/connectivy_example/PERMANENT
GISDBASE = /nfs/bparmentier-data/Data/projects/urban_garden_pursuit/data_urban_garden
LOCATION_NAME = connectivy_example
MAPSET = PERMANENT

GISDBASE = /nfs/bparmentier-data/Data/projects/urban_garden_pursuit/data_urban_garden
#LOCATION = /nfs/bparmentier-data/Data/projects/urban_garden_pursuit/data_urban_garden/connectivy_example/
LOCATION_NAME = connectivy_example
MAPSET = nyc_site_test

LOCATION_NAME=connectivy_example
GISDBASE=/nfs/bparmentier-data/Data/projects/urban_garden_pursuit/data_urban_garden
MAPSET=nyc_site_test

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

#grass75 -text ~/grassdata/mylocation/mymapset

#### Start here:
## Type this to start grass from permanent
#grass /nfs/bparmentier-data/Data/projects/urban_garden_pursuit/data_urban_garden/connectivy_example/PERMANENT
grass /nfs/bparmentier-data/Data/projects/urban_garden_pursuit/data_urban_garden/connectivy_example/nyc_site_test
grass $GISDBASE$LOCATION$MAPSET
#### Add GRASS code here:

g.gisenv #check environment
v.in.ogr -f -o input="network_nodes.shp" output="nodes_origin" #import vector files with nodes
r.in.gdal -o --overwrite input=maungawhau.tif output=r_surf #import raster image for test
          
r.info r_surf #check raster information 

#d.mon start=x0
d.mon wx0 # start a display window called wx0

d.rast r_surf # display raster
d.vect nodes_origin # add vector points on top of raster being displayed

g.region -p # examine current region

g.region rast=r_surf #make sure we are using the correct settings
g.region -p

#### Convert vector to raster
v.to.rast --overwrite input=nodes_origin use=attr output=nodes_origin_surf attribute_column=ID
r.info nodes_origin_surf

r.mapcalc 'r_friction = 1' #creates a raster with value 1
r.info r_friction
d.rast r_friction

# compute cumulative cost surfaces
r.walk -k elev=r_surf friction=r_friction output=walk.cost start_points=nodes_origin stop_points=nodes_origin lambda=1
r.walk -k --overwrite elev=r_friction friction=r_surf output=walk.cost start_points=nodes_origin stop_points=nodes_origin lambda=1

# compute shortest path from start to end points
#execGRASS()
#system("r.drain in=walk.cost out=walk.drain vector_points=end")
#system("r.drain input=walk.cost output=walk.drain vector_points=nodes_origin")
#system("r.drain input=r_surf_cost output=cost_drain vector_points=nodes_origin")

#system("r.out.gdal input=walk.drain output=path_r_walk.tif") 
#system("r.out.gdal input=walk.cost output=walk_cost.tif") 
#system("r.out.gdal input=cost_drain output=cost_drain.tif") 

############## End of script ###################
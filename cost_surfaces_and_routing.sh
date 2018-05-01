##!/bin/bash

#LOCATION = /nfs/bparmentier-data/Data/projects/urban_garden_pursuit/data_urban_garden/connectivy_example/PERMANENT
#GISDBASE = /nfs/bparmentier-data/Data/projects/urban_garden_pursuit/data_urban_garden
#LOCATION_NAME = connectivy_example
#MAPSET = PERMANENT

#LOCATION = /usr/local/share/grassdata/spearfish70/PERMANENT
#GISDBASE = /usr/local/share/grassdata
#LOCATION_NAME = spearfish70
#MAPSET = PERMANENT

#visit example
#http://ncsu-geoforall-lab.github.io/geospatial-modeling-course/grass/buffers_cost.html
#https://casoilresource.lawr.ucdavis.edu/software/grass-gis-raster-vector-and-imagery-analysis/raster-operations/simple-comparision-two-least-cost-path-approaches/

### Set path and grass environmental variable

GISDBASE=/nfs/bparmentier-data/Data/projects/urban_garden_pursuit/data_urban_garden
LOCATION_NAME=connectivy_example
MAPSET=nyc_site_test

#g.gisenv set="VARIABLE=VALUE"

# Use export to set shell environmental variable:
#export VARIABLE=value

#grass75 -text ~/grassdata/mylocation/mymapset

#### Start here:
## Type this to start grass from permanent
#grass /nfs/bparmentier-data/Data/projects/urban_garden_pursuit/data_urban_garden/connectivy_example/PERMANENT
grass /nfs/bparmentier-data/Data/projects/urban_garden_pursuit/data_urban_garden/connectivy_example/nyc_site_test
#grass $GISDBASE/$LOCATION_NAME/$MAPSET

#grass74 -c elevation.tiff -e /path/to/grassdata/test1/
#grass -c maungawhau.tif -e /nfs/bparmentier-data/Data/projects/urban_garden_pursuit/data_urban_garden/connectivy_example/nyc_site_test

g.gisenv #check environment
v.in.ogr -f -o input="network_nodes.shp" output="nodes_origin" #import vector files with nodes
r.in.gdal -o --overwrite input=maungawhau.tif output=r_surf #import raster image for test
          
r.info r_surf #check raster information 

#d.mon start=x0
d.mon start=wx0 # start a display window called wx0

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
## Does not work
#r.cost -k r_surf output=r_cost start_points=nodes_origin
#r.cost -k input=r_surf output=r_cost start_points=nodes_origin stop_points=nodes_origin

##### To make this work do this one by one:
#### Select nodes one by one
v.extract input=nodes_origin output=node_1 where="ID==1"                    
v.extract input=nodes_origin output=node_2 where="ID==2"                    
v.extract input=nodes_origin output=node_3 where="ID==3"                    

#https://casoilresource.lawr.ucdavis.edu/software/grass-gis-raster-vector-and-imagery-analysis/raster-operations/simple-comparision-two-least-cost-path-approaches/

### path 1 to 2
#r.cost -k --overwrite r_surf output=r_cost start_points=node_1 end
r.cost -k --overwrite input=r_surf output=r_cost start_points=node_1 stop_points=node_2
# getleast cost path 
r.drain --overwrite input=r_cost output=dcost_path_1_2 vector_points=node_2
r.to.vect --overwrite input=dcost_path_1_2 output=vcost_path_1_2 type=line

d.rast r_cost
d.vect vcost_path_1_2

### path 1 to 3
#r.cost -k --overwrite r_surf output=r_cost start_points=node_1 end
r.cost -k --overwrite input=r_surf output=r_cost start_points=node_1 stop_points=node_3
# getleast cost path 
r.drain --overwrite input=r_cost output=dcost_path_1_3 vector_points=node_3
r.to.vect --overwrite input=dcost_path_1_3 output=vcost_path_1_3 type=line

### path 2 to 3
#r.cost -k --overwrite r_surf output=r_cost start_points=node_1 end
r.cost -k --overwrite input=r_surf output=r_cost start_points=node_2 stop_points=node_3
# getleast cost path 
r.drain --overwrite input=r_cost output=dcost_path_2_3 vector_points=node_3
r.to.vect --overwrite input=dcost_path_2_3 output=vcost_path_2_3 type=line

#### Generate plot

d.rast r.cost
d.rast r_surf
d.vect vcost_path_1_2
d.vect vcost_path_1_3
d.vect vcost_path_2_3

#r.drain in=walk.cost out=walk.drain vector_points=end

#r.drain --overwrite input=r_surf output=cost_drain1_2 vector_points=node_2
#r.drain --overwrite input=r_surf output=cost_drain1_3 vector_points=node_3

#r.path input=r.cost start_coordinates=640206,222795 \
#    raster_path=walkpath vector_path=walkpath

#r.walk -k elev=r_surf friction=r_friction output=walk.cost start_points=nodes_origin stop_points=nodes_origin lambda=1
#r.walk -k --overwrite elev=r_surf friction=r_friction output=walk.cost start_points=nodes_origin 

# compute shortest path from start to end points
#r.drain --overwrite input=walk.cost out=walk.drain vector_points=nodes_origin

#system("r.drain in=walk.cost out=walk.drain vector_points=end")
#system("r.drain input=walk.cost output=walk.drain vector_points=nodes_origin")
#system("r.drain input=r_surf_cost output=cost_drain vector_points=nodes_origin")

#system("r.out.gdal input=walk.drain output=path_r_walk.tif") 
#system("r.out.gdal input=walk.cost output=walk_cost.tif") 
#system("r.out.gdal input=cost_drain output=cost_drain.tif") 

####### Now use randomwalk

#http://www.mergili.at/software/randomwalk_manual_20160121_main.html


#### Convert vector to raster
#v.to.rast --overwrite input=nodes_origin use=attr output=nodes_origin_surf attribute_column=ID

v.to.rast --overwrite input=node_1 use=attr output=node_1_rast attribute_column=ID                    
v.to.rast --overwrite input=node_2 use=attr output=node_2_rast attribute_column=ID                    

r.info node_2_rast
r.info node_1_rast

r.randomwalk help
r.randomwalk [-abkmnpqsvx] prefix=string [cores=integer] [cellsize=float]
[aoicoords=float,...][aoimap=name] elevation=name [releasefile=string] 
[caserules=integer,integer,...] [releasemap=name] [depositmap=name] [impactmap=name]
[probmap=name] [scoremap=name] [impactobjects=name] [objectscores=string]
models=string mparams=string [sampling=integer] [retain=float] [functype=integer]
[backfile=string] [cdffile=string] [zonalfile=string] [profile=float,...] [--verbose] [--quiet]

r.randomwalk -x prefix="rd" elevation=r_surf releasemap=node_1_rast depositmap=node_2_rast \
models=1,3,1.9,0.16,0.83,2,1,11,-9999,-9999 mparams=5,20,1000,100,10,5,2

r.randomwalk -x prefix="rd" elevation=r_surf releasemap=node_1_rast depositmap=node_2_rast models=1,3,1.9,0.16,0.83,2,1,11,-9999,-9999 mparams=5,20,1000,100,10,5,2


#releasemap

#Name of an optional input integer raster map defining the observed release area of each case. The pixel values have to correspond to the case id given in the releasefile. Zero stands for no case. With the flag x, random walks are started from all pixels of all release areas. If a release probability map (parameter probmap, flags p and x) or a release score map (parameter scoremap, flag q) is provided, the releasemap is not used for the computation, but only for visualization (flag v).

#depositmap

#Name of an optional input integer raster map defining observed deposition area of each case. Positive pixel values define areas with an observed deposition, pixel values of zero define areas without observed deposition. This map is needed for validation and visualization only (flag v). For validation purposes it is recommended to set the part of the impact area (parameter impactmap) not defined as deposition area to no data.
#############################  End of script ##################################
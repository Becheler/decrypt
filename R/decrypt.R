library("raster")
require(lattice)
require(viridis)

library("devtools")
if (!require("rspatial")) devtools::install_github('rspatial/rspatial')
library("rspatial")
library(sp)

library(dismo)
library(gstat)

plot_sampling_scheme <- function(mask, x0, r0, x, r, proj4string){
  pls <- dismo:::.generateCircles(x, d=2*r, lonlat=TRUE, crs=proj4string)
  pls0 <- dismo:::.generateCircles(x0, d=2*r0, lonlat=TRUE, crs=proj4string)
  plot(mask)
  plot(pls, add=TRUE)
  plot(pls0, add=TRUE, col="red")
}

make_movie <- function(history, ordered_times, max_N_value, working_folder){
  old_wd <- getwd()
  setwd(working_folder)
  for(i in ordered_times){
    png( paste0("N_", i, ".png"))
    print( plot(history[[i]], zlim=c(0,max_N_value), col = viridis(100), colNA="lightgrey", col) )
    dev.off() #only 129kb in size
  }
  system('n=0; ls -tr | while read i; do n=$((n+1)); mv -- "$i" "$(printf \'%04d\' "$n")"_"$i"; done')
  system("convert -delay 2 -quality 95 *.png movie.mp4")
  setwd(old_wd)
}

raw_posterior_probability <- function(data, mask, proj4string){
  points <- SpatialPoints(data[,c('lon','lat')], proj4string=proj4string)
  spdf <- SpatialPointsDataFrame(points, data)
  CA <- rasterToPolygons(mask)
  # define groups for mapping
  cuts <- c(0, 0.5, 1)
  # set up a palette of interpolated colors
  blues <- colorRampPalette(c('orange', 'blue'))
  pols <- list("sp.polygons", CA, fill = "lightgray")
  spplot(obj=spdf, zcol='p2', cuts=cuts, col.regions=blues(2), sp.layout=pols, pch=20, cex=2)
}

interpolate_posterior_probability <- function(data, mask, x0, proj4string){
  points <- SpatialPoints(data[,c('lon','lat')], proj4string=proj4string)
  spdf <- SpatialPointsDataFrame(points, data)
  CA <- rasterToPolygons(mask, fun=function(x){x >= 1 | is.na(x)})
  #CA <- rasterToPolygons(mask)
  # define groups for mapping
  cuts <- c(0, 0.5, 1)
  # set up a palette of interpolated colors
  blues <- colorRampPalette(c('orange', 'blue'))
  pols <- list("sp.polygons", CA, fill = "lightgray")
  spplot(obj=spdf, zcol='p2', cuts=cuts, col.regions=blues(2), sp.layout=pols, pch=20, cex=2)

  TA <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 ")
  dta <- spTransform(spdf, TA)
  cata <- spTransform(CA, TA)
  # Proximity variables
  v <- voronoi(dta)
  ca <- aggregate(cata)
  vca <- intersect(v, ca)
  r <- raster(cata, res=res(mask))
  vr <- rasterize(vca, r, 'p2')

  ## [inverse distance weighted interpolation]
  gs <- gstat(formula=p1~1, locations=dta)
  idw <- interpolate(r, gs)
  idwr <- mask(idw, vr)
  plot(1 - idwr)
  points(x0, pch="+", col= "red")
}

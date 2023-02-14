#' Function to create a rectangular polygon given the vertices
#' input : x.lim = x-coordinates of the polygon
#'         y.lim = y-coordinates of the polygon
#'         crs_obj = a CRS object for projection
#'         buff = numerical value of the buffer (the polygon will be enlarged by this value)
#' output : Square polygon         
Square.poly <- function(x.lim, y.lim, crs_obj, buff = 1){
  x_lim <- x.lim + c(-buff, buff)
  y_lim <- y.lim + c(-buff, buff)
  x_coords <- c(x_lim[1],x_lim[2],x_lim[2],x_lim[1],x_lim[1])
  y_coords <- c(y_lim[1],y_lim[1],y_lim[2],y_lim[2],y_lim[1])
  poly1 <- sp::Polygon(cbind(x_coords,y_coords))
  poly.bdy <- sp::Polygons(list(poly1), ID = "A")
  poly.bdy <- sp::SpatialPolygons(list(poly.bdy))
  proj4string(poly.bdy) <- crs_obj
  poly.bdy
}



#' function distance for faults
#' inputs 
#'        - x , y vector of coordinates
#'        - faults_poly SpatialPolygonDataFrame
#' output 
#'        - minimum distance between x,y and element of faults_poly.
#'                
## also takes a  (faults_poly) from which the distance is calculated.

FaultDistFn = function(x, y, faults_poly, proj_ = FALSE) {
  # turn coordinates into SpatialPoints object
  spp = SpatialPoints(data.frame(x=x, y=y)) 
  # set crs 
  proj4string(spp) = suppressWarnings(fm_crs_set_lengthunit(CRS(SRS_string='EPSG:6875'), "km"))
  # Extract values at spp coords, from our elev SpatialGridDataFrame
  FD <- gDistance(faults_poly, spp, byid = TRUE)
  ## This gives us a distance to all foult polygons - we just want the closest one
  ## Also convert to km
  min_dist <- apply(FD, 1,  min)
  
  return(min_dist)
}


FaultCov = function(x, y, faults_poly, covar) {
  # turn coordinates into SpatialPoints object
  spp = SpatialPoints(data.frame(x=x, y=y)) 
  # set crs 
  proj4string(spp) = suppressWarnings(fm_crs_set_lengthunit(CRS(SRS_string='EPSG:6875'), "km"))
  # Extract values at spp coords, from our elev SpatialGridDataFrame
  FD <- gDistance(faults_poly, spp, byid = TRUE)
  ## This gives us a distance to all foult polygons - we just want the closest one
  ## Also convert to km
  idxx <- apply(FD, 1,  which.min)
  
  covar[idxx]
}


FaultDistSingle = function(x, y, faults_poly, idx) {
  # turn coordinates into SpatialPoints object
  spp = SpatialPoints(data.frame(x=x, y=y)) 
  # set crs 
  proj4string(spp) = suppressWarnings(fm_crs_set_lengthunit(CRS(SRS_string='EPSG:6875'), "km"))
  # Extract values at spp coords, from our elev SpatialGridDataFrame
  FD <- gDistance(faults_poly, spp, byid = TRUE)
  ## This gives us a distance to all foult polygons - we just want the closest one
  ## Also convert to km
  apply(FD, 1,  function(x) x[idx])
  
  
}










library(grid)
library(gridExtra)
library(cowplot)
library(RColorBrewer)
##library(VGAM)


#' Create polygons from point locations that describe polygon edges
#'
#' Useful where faults are stored as points in a dataframe (e.g. UCERF 3 fault geometries).
#' See LineMaker for line alternative.
#' @param rownum row number for fault of interest
#' @param cols columns in dataframe to be used in creating polygons
#' @param Geom fault geometry stored in a dataframe
#' @return returns SpatialPolygonsDataFrame object
#' @keywords polygon
#' @export
#' @examples
#' To run with UCERF3 geometry, applying to all faults and skipping the first 8 columns:
#'
#' NFaults <- length(UCERFFaultGeom$Name)
#' #x <- seq(1, NFaults,1)
#' ### Apply faultgeometry function to all faults, make list of fault polygons
#' FaultList <- lapply(7:NFaults, PolygonMaker, Geom=UCERFFaultGeom, cols=9:82)
#' ### Make list into SpatialPolygons (sp faff)
#' FaultPolys <- SpatialPolygons(lapply(FaultList, function(x){x@polygons[[1]]}))
#'

PolygonMaker <- function(rownum, cols, Geom){
  ### Polygons are stored in a row with alternating lat/lon columns and NAs when no further points are necessary
  ### First column is fault name, keep this for polygon ID
  PolyName <- as.character(Geom$Name[rownum])
  ### as.numeric is probably not necessary here, because if it is you might have factors and you should be worried about if this will mess everything up...
  PolyCoords <- as.numeric(Geom[rownum, cols])
  ### Drop NAs (these are just here so num of columns is consistent)
  PolyCoordsOnly <- PolyCoords[is.na(PolyCoords) == FALSE]
  ### Extract alternating lats and lons
  PolyLats <- PolyCoordsOnly[seq(1, length(PolyCoordsOnly), 2)]
  PolyLons <- PolyCoordsOnly[seq(2, length(PolyCoordsOnly), 2)]
  ### Convert lats/lons into Polygon object
  ### NB: matrix needs c - if you don't explicitly hand it a list it will only use one dataset and just ignore the other
  Poly <- Polygon(matrix(c(PolyLons, PolyLats), nrow=length(PolyLats)))
  ### Convert Polygon object to Polygons object, attach name of Polygon/fault
  P2 <- Polygons(list(Poly), ID=PolyName)
  ### Attach CRS to polygon
  P3 <- SpatialPolygons((list(P2)),  proj4string = CRS("+proj=longlat +datum=WGS84"))
  ###  - Should get this into projected CR so we are working in km/m rather than degrees - can use sp_transform but you'll need to figure out suitable projection
  return(P3)
}


#' Create line objects from point locations that describe points of the line
#'
#' Useful where faults are stored as points in a dataframe (e.g. UCERF 3 fault geometries).
#' See PolygonMaker for polygon alternative.
#' @param rownum row number for fault of interest
#' @param Geom fault geometry stored in a dataframe
#' @param cols columns in dataframe to be used in creating polygons
#' @return returns SpatialLinesDataFrame object
#' @keywords line
#' @export
#' @examples
#' To run with UCERF3 geometry, applying to all faults and skipping the first 7 rows:
#' NFaults <- length(UCERFFaultGeom$Name)
#' FaultLineList <- lapply(1:NFaults, LineMaker, Geom=UCERFFaultGeom, cols=9:82)
#' FaultLines <- SpatialLines(lapply(FaultLineList, function(x){x@lines[[1]]}))

LineMaker <- function(rownum, Geom, cols){
  ### Polygons are stored in a row with alternating lat/lon columns and NAs when no further points are necessary
  ### First column is fault name, keep this for polygon ID
  FaultName <- as.character(Geom$Name[rownum])
  ### as.numeric is probably not necessary here, because if it is you might have factors and you should be worried about if this will mess everything up...
  FaultCoords <- as.numeric(Geom[rownum,cols])
  ### Drop NAs (these are just here so num of columns is consistent)
  FaultCoordsOnly <- FaultCoords[is.na(FaultCoords) == FALSE]
  ### Extract alternating lats and lons
  PolyLats <- FaultCoordsOnly[seq(1, length(FaultCoordsOnly), 2)]
  PolyLons <- FaultCoordsOnly[seq(2, length(FaultCoordsOnly), 2)]
  ### Convert lats/lons into lines
  ### NB: matrix needs c - if you don't explicitly hand it a list it will only use one dataset and just ignore the other
  Line <- Line(matrix(c(PolyLons, PolyLats), nrow=length(PolyLats)))
  ### Convert Polygon object to Polygons object, attach name of Polygon/fault
  L2 <- Lines(list(Line), ID=FaultName)
  ### Attach CRS to polygon
  L3 <- SpatialLines((list(L2)),  proj4string = CRS("+proj=longlat +datum=WGS84"))
  ###  - Should get this into projected CR so we are working in km/m rather than degrees - can use sp_transform but you'll need to figure out suitable projection
  return(L3)
}


#' Calculate on and off fault events and create a distance from fault map, which can be saved to specified directory
#'
#' This function identifies points that occur on and off-fault, and returns a distance from fault map.
#' @param Points spatial points dataframe of events/covering map area, should be in EPSG CRS.
#' @param FaultMap SpatialPolygonsDataFrame with fault polygons to be used. Should contain a column called `Fault` that names individual fault segments (numbers will also work - this is used to identify which fault an event is on)
#' @param Res Resolution of output in m. Default is 2500m.
#' @param saveDir Optional save directory for output, if supplied raster will be saved as rds
#' @param EPSG specify EPSG code, defaults to California 3310, see https://spatialreference.org/ref/epsg/ for other areas
#' @param xlim limits for x-axis of fault-distance raster
#' @param ylim limits for the y-axis of fault-distance raster
#' @return returns raster map of fault distances in metres
#' @keywords faultdistance
#' @export
#' @examples
#' Fd <- FaultDistMap(UCERFCat, FTDF2)
#'

FaultDistMap <- function(Points, FaultMap, Res=2500, saveDir = 0, EPSG=3310, xlim, ylim){
  PolyTest2 <- over(Points, FaultMap)
  Points$Fault <- PolyTest2$Name
  OnFault <- subset(Points, is.na(EQs$Fault) == FALSE)
  OffFault <- subset(Points, is.na(EQs$Fault) == TRUE)
  
  FLDF <- st_as_sf(FaultMap)
  FLDF <- st_transform(FLDF, EPSG)
  OF <- st_as_sf(OffFault)
  OF <- st_transform(OF, EPSG)
  
  dd2 <- gDistance(as(FLDF, "Spatial"), as(OF, "Spatial"), byid = TRUE)
  min_dist <- apply(dd2, 1, min)
  
  grdF3= expand.grid(x=seq(from=xlim[1], to= xlim[2], by=0.1), y=seq(from=ylim[1], to=ylim[2], by=0.1))
  coordinates(grdF3) = c("x", "y")
  proj4string(grdF3) = CRS(proj4string(FaultMap))
  # grid expand 'by' arguments don't seem to affect raster, so set resolution in raster layer.
  r2 = raster(grdF3, resolution= 0.05)
  #rm(grdF3)
  # ## Can potentially increase this resolution a bit but it takes a while
  r2 <- projectRaster(r2, crs=crs(FLDF), res=Res)
  dr2 <- gDistance(as(FLDF, "Spatial"), as(r2, "SpatialPoints"), byid = TRUE)
  r2[] <- apply(dr2, 1, min)
  #rm(dr2)
  #ggplot()+ gg(r2)+scale_fill_gradientn(colors=RColorBrewer::brewer.pal(9, "YlOrRd"))  +  ggtitle("Distance to fault (m) - resolution 2.5km")+ labs(x="", y="") + theme_classic() + coord_equal()
  # + gg(as(OF, "Spatial"))
  #r2km <- r2
  r2km <- setValues(r2, r2[]/1000)
  proj_rast <- projectRaster(r2km, crs=CRS(proj4string(FaultMap)))
  FaultDist <- as(proj_rast, 'SpatialPixelsDataFrame')
  if (saveDir != 0){
    saveRDS(FaultDist, file=paste0(saveDir,"FaultDist2pt5kmProj.rds"))
  }
  return(FaultDist)
}


#' Apply uniform buffer to fault polygons
#'
#' Simple function to apply uniform buffer to fault polygons using st_buffer, clean any intersections with clgeo_Clean, and crop large datasets to specifiec boundary.
#' @param Geom Fault geometry as SpatialPolygonsDataFrame
#' @param width width of buffer to be applied in metres
#' @param bdy boundary for cropping. Must be supplied as SpatialPolygon. Default is not to crop.
#' @param EPSG EPSG code for area of application. Defaults to California (3310), see https://spatialreference.org/ref/epsg/ for other areas
#' @return returns SpatialPolygonsDataFrame with specified buffer
#' @keywords buffer
#' @export
#' @examples
#' bdy <- SpatialPolygons(list(Polygons(list(Polygon(boundary[[1]]$loc)), ID="P1")), proj4string=CRS(proj4string(FTDF2)))
#' UnifBuff <- UniformBuffer(FTDF2, 5000, bdy=bdy)


UniformBuffer<- function(Geom, width, bdy = 0, EPSG = 3310){
  #FPr <- st_transform(Geom, EPSG)
  FBuff <- st_buffer(Geom, width)
  length(which(st_is_valid(FBuff) == FALSE))
  #FB <-  st_transform(FBuff, "+proj=longlat +datum=WGS84")
  FB <- as(Geom, "Spatial")
  
  #rm(QFBuff, QFB, QFPr, QFaults)
  ## have to clean up polygons to avoid one weird self-intersection which doesn't show up in sf
  FBc <- clgeo_Clean(FB)
  ### If you've skipped the cleaning step, this may crash Rstudio if you have any invalid geometry. Don't say I didn't warn you!
  ### But actually 'crop' just seems to throw an error, whereas gIntersection just killed everything.
  ### Can also use boundary polygon if you have one, but I'm not sure where this would be necessary?
  if (bdy != 0){
    FBc <- crop(FBc, bdy)
  }
  return(FBc)
}

#' Dip-dependent fault polygon buffer function
#'
#' Apply dip-dependent buffers to fault polygons. Will apply dip in metres, so check CRS!
#' Uses UCERF3 buffers: fault buffer scales with fault dip - 0km at 50 dip to 12km at 90 dip
#' @param PolygonsData a dataset of polygons stored as a spatial polygons object.
#' @param Dip a vector of dips. Must be the same length as PolygonsData
#' @return returns new polygons, buffered according to dip and stored as a SpatialPolygonsDataFrame
#' @keywords buffers
#' @export
#' @examples
#' FaultPolyBuffer <- BufferFunction(FL_T, FL_T$Dip)
#' Use st_transform to transform to EPSG (3310 for state of California, see https://spatialreference.org/ref/epsg/ and vice-versa

DipBufferFunction <- function(PolygonsData, Dip){
  DipAng <- vector(mode="numeric", length=nrow(PolygonsData))
  BufferSize <- vector(mode="numeric", length=nrow(PolygonsData))
  # Units of dist are m
  #BufferSize <- DipAng*116.7 - IF you start buffering at 0 (degrees and km) and extend to 15km at 90 degrees
  ## Here start the buffer if dip > 50 and buffer to 12km at dip=90
  ## But something really dodgy happens if buffersize=0 so set to 0.1 (that's like 10cm, should not make any difference to data!)
  for (i in 1:nrow(PolygonsData)){
    DipAng[i] <- Dip[i]
    if (DipAng[i] > 50){
      BufferSize[i] <- (DipAng[i]-50)*0.3
    }
    else {BufferSize[i] = 0.1}
  }
  #print(BufferSize)
  BufferPoly <- st_buffer(PolygonsData, dist=BufferSize*1000)
  return(BufferPoly)
}


#' Make a SpatialPolygons object from a dataframe 
#'
#' This function converts from a dataframe of x and y coordinates to a SpatialPolygons object
#' This is basically a nested series of conversions from a dataframe to a Polygon to a Polygons list to a SpatialPolygons object.
#' @param coords a dataframe of x and y coordinates describing a spatial polygon
#' @param crs desired crs for the object
#' @return returns SpatialPolygons object
#' @keywords SpatialPolygons
#' @export
#' @examples
#' RELM_poly <- read.delim("RELMTestingPolygon.txt", sep="", header = FALSE)
#' RELM <- spatial_poly_from_df(RELM_poly, CRS(SRS_string = wkt))

spatial_poly_from_df <- function(coords, crs){
  SpatialPolygons((list( Polygons(list(Polygon(coords)), ID="A1"))),  proj4string = crs)
  }

square_poly_from_bbox <- function(box, crs.obj, buff = 1){
  x.lim <- box[1,] + c(-buff, buff)
  y.lim <- box[2,] + c(-buff, buff)
  df.coords <- data.frame(x = c(x.lim[1], x.lim[2], x.lim[2], x.lim[1], x.lim[1]),
                          y = c(y.lim[1], y.lim[1], y.lim[2], y.lim[2], y.lim[1]))
  spatial_poly_from_df(df.coords, crs.obj)
}
#' Project inlabru model output to uniform grid for CSEP testing 
#'
#' This function takes in an inlabru fit object and model description and parameters required to generate a pyCSEP gridded-forecast
#' @param lgcp_fit an inlabru fit object
#' @param lgcp_model a description of an inlabru model
#' @param b_poly_latlon a boundary polygon in lat-lon coordinates for pyCSEP testing
#' @param dh grid spacing for pyCSEP forecast (in lat/lon)
#' @param mag_min minimum magnitude for the forecast
#' @param mag_max maximum magnitude in forecast 
#' @param b_est estimated b-value to be used in forecast
#' @param mesh inlabru mesh object used for the fitted model
#' @return returns a data frame consistent with those required by pyCSEP
#' @export
csep_grid_wrapper <- function(lgcp_fit, lgcp_model, b_poly_latlon, dh, mag_min, mag_max, b_est, mesh){
  ## Set up magnitude breaks
  mag.break <- seq(mag_min, mag_max, by=0.1)
  ## Find bbox of boundary polygon
  b <- as.data.frame(bbox(b_poly_latlon))
  ### Make uniform grid based on limits + 0.5*spacing to get midpoints expected in RELM polygon
  mid_grd <- expand.grid(x=seq(from=(b$min[1]-(0.5*dh)), to= (b$max[1]+(0.5*dh)), by=dh), y=seq(from=(b$min[2]-(0.5*dh)), to=(b$max[2] +  (0.5*dh)), by=dh))
  
  # Make spatial object
  coordinates(mid_grd) = ~x+y
  proj4string(mid_grd) <- CRS(proj4string(b_poly_latlon))
  # Only keep midpoints in the polygon
  pts_grd_latlon <- crop(mid_grd, b_poly_latlon)
  pts_grd_km <- spTransform(pts_grd_latlon, crs_Cal_km)
  
  ## Predict at mesh locations
  Pr_Num <- predict(lgcp_fit, mesh, lgcp_model)
  ## Project to  grid
  proj <- INLA::inla.mesh.project(mesh, pts_grd_km)
  ## Get values at grid midpoints
  N <- as.vector(proj$A %*% Pr_Num$mean)
  
  ## a value from N scaled by area of grid cells
  a2 <- log(N*(1/(dh^2)))/log(10)
  ## number of bins we want
  n_mb <- length(mag.break)
  
  ## Distribute events/bin over different magnitude bins according to GR
  bmt <- vector()
  for(i in 1:length(pts_grd_km)){
    Ns <-  a2[i] - b_est*(mag.break - 4.95)
    M_D <- c(abs(diff(10^Ns)), 10^Ns[n_mb])
    bmt <- c(bmt, M_D)
  }
  
  ## lower latitude values for grid cells
  lats_1 <- rep((as.numeric(pts_grd_latlon$y) - 0.5*dh), each=n_mb)
  ## Upper values
  lats_2 <- lats_1 + dh
  
  ## Same for lons
  lons_1 <- rep((as.numeric(pts_grd_latlon$x) - 0.5*dh), each=n_mb)
  lons_2 <- lons_1 + dh
  
  ## Set up depth bins (not currently used)
  D1 <- rep(0, length(lats_2))
  D2 <- rep(30, length(lats_2))
  
  ## Set up magnitude bins (upper bounds)
  M2 <- mag.break[2:n_mb]
  M2[n_mb] <- 10
  
  mags_1 <- rep(mag.break, length(pts_grd_latlon$x))
  mags_2 <- rep(M2, length(pts_grd_latlon$x))
  
  ## RELM flag, completely unused but expected in input.
  flag <- rep(1, length(lats_1))
  
  ### Make dataframe of results
  DFT_bm <- as.data.frame(cbind(lons_1, lons_2, lats_1, lats_2, D1, D2, mags_1, mags_2, bmt, flag))
  
  # Sort by lon (probably no longer necessary!)
  grd_forecast <- DFT_bm %>% arrange(lons_1)
  return(grd_forecast)
}

#' Generate catalogue-type forecast from inlabru fitted intensity posterior sample
#'
#' This function takes in posterior intensity sample from the inlabru `generate' function (or a mean/mode/percentile if you feel so inclined)
#' This sample is then used to create an intensity field from which we sample points with a rejection sampler
#' @param loglambda log intensity values from a generate call or a fitted model object
#' @param bdy a boundary polygon over which to generate samples
#' @param mesh inlabru mesh object used for the fitted model
#' @param crs coordinate reference system to be used - must be consistent with mesh, bdy and loglambda
#' @param n_events number of events to sample for each catalogue. Currently assumes this is a Poisson rate and randomly samples exact number
#' @param b_val a b-value or b-value distribution to be used in assigning magnitudes
#' @param m_min minimum magnitude in forecast model, used in magnitude creation by TapGRM
#' @return returns a data frame of a synthetic catalogue
#' @export
point_sampler_best <- function(loglambda, bdy, mesh,  crs, num_events, b_val, m_min, corner_mag){
  ## Number of events for single catalogue from a poisson distribution with lambda = num_events
  cat_n <- rpois(1, num_events)
  if (cat_n < 1){cat_n =1}
  ll <- as.data.frame(cbind(mesh$loc[,1], mesh$loc[,2], loglambda))
  colnames(ll) <- c("x", "y", "loglambda")
  coordinates(ll) <- c("x", "y")
  proj4string(ll) <- crs
  ## Mesh often extends outside our area of interest, crop so we are only considering sampling area bdy
  llcrop <- crop(ll, bdy)
  loglambda_max <- max(llcrop$loglambda)
  
  #print(loglambda_max)
  ## Set up a spatialpoints dataframe for our results
  samp.points <- SpatialPoints(data.frame(x = 0, y = 0))
  samp.points$mags <- 0
  samp.points <- samp.points[-1,]
  proj4string(samp.points) <- crs
  #print(proj4string(samp.points))
  num <- 0
  n1 <- 300000
  n2 <- 5000000
  ## To sample the correct number of points, keep going until the num >= cat_n
  while (num < cat_n){
    pts <- spsample(bdy, n1, "random")
    #pts <- spTransform(points, crs)
    proj <- INLA::inla.mesh.project(mesh, pts)
    lambda_ratio <- exp(as.vector(proj$A %*% loglambda) - loglambda_max)
    #print(max(lambda_ratio))
    #print(sum(proj$ok))
    keep <- proj$ok & (runif(n1) <= lambda_ratio)
    kept <- pts[keep]
    
    ## If we didn't get any events, run again with more sampled points
    while (length(kept) == 0){
      #print("No events kept - trying more locations")
      pts <- spsample(bdy, n2, "random")
      
      proj <- INLA::inla.mesh.project(mesh, pts)
      lambda_ratio <- exp(as.vector(proj$A %*% loglambda) - loglambda_max)
      #print(max(lambda_ratio))
      #print(sum(proj$ok))
      keep <- proj$ok & (runif(n2) <= lambda_ratio)
      kept <- pts[keep]
    }
    kept$mags <- rep(0, length(kept))
    
    samp.points <- rbind(samp.points, kept)
    num <- length(samp.points)
    
  }
  
  ## Keep exactly cat_n points, choose these randomly from all of the points we've kept so far
  kp <- sample(seq(1, length(samp.points), by=1), cat_n, replace=FALSE)
  samp.points <- samp.points[kp,]
  ## Get magnitudes for this catalogue
  ## If b_val is a list of possible b-values, select one at random. If it's just one b-value, use that.
  b_est <- b_val[runif(1, 1, length(b_val))]
  
  samp.points$mags <- TapGRM(cat_n, b_est, 8, m_min)
  
  return(samp.points)
}


#' Sample magnitudes from a tapered Gutenberg-Richter distribution
#'
#' This function takes a number of required magnitudes, a b-value and corner- and minimum-magnitudes and returns sampled magnitudes consistent with a tapered-GR distribution 
#' Uses method of Vere-Jones et al 2001 (https://doi.org/10.1046/j.1365-246x.2001.01348.x)
#' @param n number of magnitudes to sample
#' @param b b-value for magnitude sampling
#' @param corner_mag a corner magnitude value
#' @param m_min minimum magnitude for magnitude creation
#' @return returns n magnitudes from a TGR distribution
#' @export
TapGRM <- function(n, b, corner_mag, m_min){
  ## Convert Mc to moment and b to beta
  Mt <- 10**((3/2)*m_min + 9.1)
  Mc <- 10**((3/2)*corner_mag + 9.1)
  beta <- (2/3)*b
  
  R1 <- runif(n, 0, 1)
  R2 <- runif(n, 0, 1)
  
  M1 <- Mt*R1**(-1/beta)
  M2 <- Mt - Mc*log(R2)
  
  # Pick minimum value of two sampled options
  mom <- pmin(M1, M2)
  
  ## convert moments back to magnitudes
  ms <- (2/3)*(log10(mom) - 9.1)
  
  return(ms)
}

#' Calculate frequency-magnitude distribution, called by other Mc functions
#' @param mag vector of magnitudes
#' @param mbin binsize
#' @return dataframe with magnitude bins (mi), cumulative (cum) and non-cumulative (noncum) frequency-magnitude
#' @export
fmd <- function(mag,mbin){
  mi <- seq(min(round(mag/mbin)*mbin), max(round(mag/mbin)*mbin), mbin)
  nbm <- length(mi)
  cumnbmag <- numeric(nbm)
  nbmag <- numeric(nbm)
  for(i in 1:nbm) cumnbmag[i] <- length(which(mag > mi[i]-mbin/2))
  cumnbmagtmp <- c(cumnbmag,0)
  nbmag <- abs(diff(cumnbmagtmp))
  res <- list(m=mi, cum=cumnbmag, noncum=nbmag)
  return(res)
}

#' Calculate Mc by method of maximum curvature
#' @param mag vector of magnitudes
#' @param mbin binsize
#' @return list containing Mc estimate 
#' @export
maxc <- function(mag,mbin){
  FMD <- fmd(mag,mbin)
  Mc <- FMD$m[which(FMD$noncum == max(FMD$noncum))[1]]
  return(list(Mc=Mc))
}

#' Calculate Mc by method of modified b-value stability
#' Mc estimation by b-val Stability (MBS) [Cao & Gao, 2002], modification with Shi & Bolt [1982] uncertainty [Woesner & Wiemer, 2005]
#' This version allows the user to specify the number of required bins, which is sometimes necessary for catalogues with high Mc and many small events
#' @param mag vector of magnitudes
#' @param mbin binsize
#' @param n number of bins to consider
#' @param calc_inc_size size of increments to estimate Mc over, defaults to 0.1
#' @return Mc estimate, list of trial Mc values and associated b value at each Mc, uncertainty in b from Shi and Bolt (1982) and 5-bin average b-value (bave)
#' #' @export
mbs_mod <- function(mag,mbin, n, calc_inc_size = 0.1){
  ### If n = number of bins to consider, n-5 is number of averages
  nl <- n-5
  McBound <- maxc(mag,mbin)$Mc
  
  Mco <- McBound+(seq(n)-1)/(1/calc_inc_size)
  
  bi <- numeric(n); unc <- numeric(n)
  for(i in 1:n){
    indmag <- which(mag > Mco[i]-mbin/2)
    nbev <- length(indmag)
    bi[i] <- log10(exp(1))/(mean(mag[indmag])-(Mco[i]-mbin/2))
    unc[i] <- 2.3*bi[i]^2*sqrt(sum((mag[indmag]-mean(mag[indmag]))^2)/(nbev*(nbev-1)))
  }
  bave <- numeric(nl)
  for(i in 1:nl) bave[i] <- mean(bi[i:(i+5)])
  
  dbi_old <- abs(diff(bi))
  indMBS_old <- which(dbi_old <= 0.03)
  
  dbi <- abs(bave[1:nl]-bi[1:nl])
  indMBS <- which(dbi <= unc[1:nl])
  Mc <- Mco[indMBS[1]]
  return(list(Mc=Mc, Mco=Mco, bi=bi, unc=unc, bave=bave))
}


#' Function to calculate GR parameters a and b
#' 
#' @param mag vector of magnitudes
#' @param Mc estimate of completeness magnitude
#' @param mbin binsize
#' @return list containing maximum likelihood a and b-value estimates for GR FMD 
#' @export
calc_GR_params <- function(mag, Mc, mbin) {
  indmag <- which(mag > Mc-mbin/2)
  b <- log10(exp(1))/(mean(mag[indmag])-(Mc-mbin/2))
  a <- log10(length(indmag))+b*(Mc - mbin/2)
  params <- c(a, b)
  return(params)
}

#' Build a distribution of b-values using magnitude uncertainty
#'
#' This function takes a vector of earthquake magnitudes and some uncertainty on the magnitudes and simulates a set of b-values consistent with this uncertainty.
#' Magnitude uncertainties can be supplied as a single value representing the mean of a normal distribution (uncert="normal) or Laplace distribution (uncert="Laplace")
#' Or as a vector of values which correspond to individual magnitudes.
#' For each simulation, the completeness magnitude (mc) is calculated using the modified b-value stability method (see mbs function description)
#' 
#' @param mags vector of earthquake magnitudes
#' @param mbin size of magnitude bins, required for calculating mc at each simulation
#' @param uncert Choice of uncertainty distribution. Options are "normal" and "Laplace".
#' @param vals uncertainty values as a single standard deviation (uncert="normal") or s (uncert="Laplace"), or a vector with a value for each magnitude
#' @param iter number of required iterations (default = 1000)
#' @return returns a list of (n=iter) Mc and b-values consistent with the data
#' @export
bval_distro <- function(mags, mbin=0.01, uncert="normal", sd = 0.2, iter=1000){
  nmags <- length(mags)
  b_vals <- vector(mode="numeric", length=iter)
  Mc_new <- vector(mode="numeric", length=iter)
  
  Mc_base <- mbs_mod(mags, mbin, n=40)$Mc
  Mc_new[1] <- Mc_base
  b_vals[1] <- calc_GR_params(mags, Mc_base, mbin)[2]
  for (i in 2:iter){
    if (uncert == "normal"){
      new_mags <- rnorm(nmags, mean=mags, sd = sd)
    }
    
    if (uncert == "laplace"){
      new_mags <- rlaplace(nmags, m=mags, s=sd)
    }
    Mc_new[i] <- mbs_mod(new_mags, mbin, n=40)$Mc
    
    b_vals[i] <- calc_GR_params(new_mags, Mc_new[i], mbin)[2]
    
  }
  
  return(list(b_vals = b_vals, Mc_new=Mc_new))
}

#' Convert Mw to Ml magnitudes
#' 
#' A function to convert from Mw to Ml given Mw = b*Ml + a for the linear relationship between magnitude scales.
#' Parameters a and b, describing the intercept and gradient respectively, are required.
#' Defaults to HORUS/Iside conversion parameters from Lolli et al
#' @param mags_mw magnitudes to be converted in Mw
#' @param a intercept term for conversion of parameters between Ml and Mw
#' @param b gradient term for conversion of parameters between Ml and Mw
#' @return returns new magnitudes in Ml
#' @export
Mw_to_ML <- function(mags_mw, a=-0.164, b=1.066){
  Ml = (mags_mw - a)/b
  return(Ml)
}


#' Function to generate and save forecasts
#' 
#' This function is an all-in-one solution to sampling and saving catalogue forecasts. It essentially combines the inlabru generate call with the point_samplers function to generate stochastic samples in one call.
#' Inputs are the fitted inlabru model, a model formula and the inputs for point_sampler, as well as details of magnitude type, number of required forecast simulations and a name to save the output under.
#' 
#' @param model_fit a fitted inlabru model object to use as our spatial model
#' @param model_formula linear predictor for the model_fit object, required for proper sampling
#' @param n_samples the number of samples to draw from the posterior, so the number of required forecast simulations
#' @param mesh inlabru mesh object the model was contstructed on
#' @param bdy a boundary object for the edge of our forecast
#' @param crs local crs the model is constructed in
#' @param num_events the rate of a Poisson distribution which we'll use as the number of required events per forecast (output of inlabru prediction call)
#' @param b_val a b-value or b-value distribution from which to sample when assigning magnitudes to events
#' @param m_min minimum magnitude for TGR distribution
#' @param mag_type types of magnitude to be returned. Options are Mw or Ml
#' @param forecast_name name under which to save forecast as a string
#' @return saves forecast, no returned output
#' @export

full_forecast <- function(model_fit, model_formula, n_samples, mesh, bdy, crs, num_events, b_val, m_min, corner_mag, mag_type="Mw", forecast_name){
  
  generated_samples <- generate(model_fit, mesh, model_formula, n.samples = n_samples)
  
  samp_points <- future_apply(generated_samples, 2, point_sampler_best, bdy=bdy, mesh=mesh, crs=crs, num_event=num_events, b_val=b_val, m_min=4, corner_mag=corner_mag, future.seed = 5)
  
  Cats <- do.call(rbind, lapply(1:n_samples, function(x){
    data = spTransform(samp_points[[x]], crs_wgs84)
    SpatialPointsDataFrame(coords = data@coords,
                           data = as.data.frame(data))
  })
  )
  cat_lens <- as.numeric(unlist(lapply(samp_points, length)))
  Cats$cat_id <- rep(seq(1, length(samp_points), by = 1), times=cat_lens)
  
  
  times <- rep("2020-02-02T01:01:01.020000", length(Cats$mags))
  depth <- rep(10, length(Cats$mags))
  event_id <- seq(1, length(Cats$mags))
  stoch_cat_set <- as.data.frame(cbind(Cats$x, Cats$y, Cats$mags, times, depth, Cats$cat_id, event_id))
  
  if (mag_type == "Mw"){
  write.table(format(stoch_cat_set, digits = 6), forecast_name, quote = FALSE, row.names = FALSE, col.names = FALSE, sep=",")
  }
  else if (mag_type == "Ml"){
    Ml_mags <- Mw_to_ML(as.numeric(stoch_cat_set$V3))
    stoch_cat_setMl <- stoch_cat_set
    stoch_cat_setMl$V3 <- Ml_mags
    write.table(format(stoch_cat_setMl, digits = 6), forecast_name, quote = FALSE, row.names = FALSE, col.names = FALSE, sep=",")
  }
}


#' Plot pairplot of collection of models
#' 
#' Plots comparison plot of a list of models with median posterior plot as diagonal, pairwise median differences in the top right and pairwise variance differencs in the bottom left.
#' Median plots have consistent colour scale with colour scheme specified in function call. Pairwise difference plots use blue-white-red to show negative - zero - positive changes.
#' 
#' @param pred_list list of inlabru prediction objects to be compared
#' @param pred_names names for each prediction (list of strings)
#' @param colour_scheme an RColorBrewer colour scheme for median plots (comparison plots use red-blue colour scheme for variation)
#' @param poly a polygon object to include in the figure - this should be optional but is currently not. Sorry.
#' @param med_lims limits for the median differences (currently defaults to (-2, 2))
#' @param var_lims limits for varaince differences (currently defaults to (-3, 5))
#' @return plot comparing all models in list
#' @export

comp_pairplots <- function(pred_list, pred_names, colour_scheme, poly, med_lims =c(-2, 2), var_lims = c(-3, 5)){
  n <- length(pred_list)
  
  ## Set up colour scheme for median plots
  ranges <- list()
  for (j in 1:n){
    ranges = append(ranges, range(pred_list[[j]]$median))
  }
  
  csc <- scale_fill_gradientn(colours = brewer.pal(9, colour_scheme), limits = range(ranges))
  plots <- list()
  pl_nm <- 1
  for (i in 1:n){
    for (k in 1:n){
      if (i < k){
        pred_list[[i]]$med_diff <- pred_list[[i]]$median - pred_list[[k]]$median
        plots[[pl_nm]] <- as_grob(ggplot() + gg(pred_list[[i]]['med_diff']) + geom_sf(data = st_as_sf(poly), fill=NA) +  labs(x="Easting", y="Northing") + theme_classic() + coord_sf() + scale_fill_gradient2(low="blue", mid="white", high="red", limits=med_lims) + ggtitle(paste0(pred_names[i], "-", pred_names[k]))+ theme(plot.title = element_text(size=8))) 
        pl_nm <- pl_nm + 1
      }
      else if (i > k){
        pred_list[[i]]$var_diff <- pred_list[[i]]$var - pred_list[[k]]$var
        plots[[pl_nm]] <- as_grob(ggplot() + gg(pred_list[[i]]['var_diff']) + geom_sf(data = st_as_sf(poly), fill=NA) +  labs(x="Easting", y="Northing") + theme_classic() + coord_sf() + scale_fill_gradient2(low="blue", mid="white", high="red", limits=var_lims) + ggtitle(paste0(pred_names[i], "-", pred_names[k], " var")) + theme(plot.title = element_text(size=8))) 
        pl_nm <- pl_nm + 1
      }
      else{
        plots[[pl_nm]] <- as_grob(ggplot() + gg(pred_list[[i]]['median']) + geom_sf(data = st_as_sf(poly), fill=NA) +  labs(x="Easting", y="Northing") + theme_classic() + coord_sf() + csc + ggtitle(pred_names[i]) + theme(plot.title = element_text(size=8))) 
        pl_nm <- pl_nm + 1
      }
    }
  }
  
  
  ##core = list(fg_params=list(cex = 1.0)),
  ##colhead = list(fg_params=list(cex = 0.5)),
  ##rowhead = list(fg_params=list(cex = 0.5))
  combine <- rbind(tableGrob(t(c(pred_names[1:n])), theme = ttheme_minimal(base_size= 12), rows = ""), 
                   cbind(tableGrob(pred_names[1:n], theme = ttheme_minimal(base_size = 12)), 
                         arrangeGrob(grobs = plots),  size = "last"), size = "last")
  grid.newpage()
  grid.draw(combine)
  
  #grid.arrange(grobs=lapply(plots, grobTree), ncol=n)
  
}

#' Function to plot median and variances side by side
#'
#' @param pred_list list of inlabru prediction objects to be compared
#' @param pred_names names for each prediction (list of strings)
#' @param colour_scheme an RColorBrewer colour scheme for plots 
#' @param poly a polygon object to include in the figure - this should be optional but is currently not. Sorry.
#' @return plot with median log intensity and variance columns for each model in list
#' @export
comparison_plots <- function(list_of_predictions, pred_names, colour_scheme, poly){
  n <- length(list_of_predictions)
  
  ## Set up colour scheme for median plots
  ranges_med <- list()
  ranges_var <- list()
  for (j in 1:n){
    ranges_med = append(ranges_med, range(pred_list[[j]]$median))
    ranges_var = append(ranges_var, range(pred_list[[j]]$var))
  }
  
  csc_range_med <- scale_fill_gradientn(colours = brewer.pal(9, colour_scheme), limits = range(ranges_med))
  csc_range_var <- scale_fill_gradientn(colours = brewer.pal(9, colour_scheme), limits = range(ranges_var))
  
  pl_nm <- 1
  plots <- list()
  for (i in 1:n){
    plots[[pl_nm]] <- as_grob(ggplot() + gg(pred_list[[i]]['median']) + gg(poly) +  labs(x="Easting", y="Northing") + theme_classic() + coord_equal() + csc_range_med + ggtitle(pred_names[i]) + theme(text = element_text(size = 8), plot.title = element_text(size = 8)))
    pl_nm <- pl_nm + 1
    plots[[pl_nm]] <- as_grob(ggplot() + gg(pred_list[[i]]['var']) + gg(poly) +  labs(x="Easting", y="Northing") + theme_classic() + coord_equal() + csc_range_var + ggtitle(pred_names[i]) + theme(text = element_text(size = 8), plot.title = element_text(size = 8)))
    pl_nm <- pl_nm + 1
    
  }
  
  #d <- grid.arrange(grobs=lapply(plots, grobTree), ncol=2)
  col.titles <- c("median", "variance")
  #tt <- ttheme_default()
  #grid.table(d, theme=tt)
  #return(plots)
  
  #grid.arrange(grobs=lapply(plots, grobTree), ncol=2)
  
  grid.arrange(grobs=lapply(c(1,2), function(i) {
    arrangeGrob(grobs=plots[seq(i, 2*n, by=2)], top=col.titles[i], ncol=1)
  }), ncol=2)
  
}
  
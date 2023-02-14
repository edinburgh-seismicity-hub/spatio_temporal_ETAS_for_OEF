library(sf)
library(inlabru)
library(INLA)
library(dplyr)
library(data.table)
library(metR)
library(matrixStats)
library(mvtnorm)
library(MASS)
library(raster)
library(rgeos)
library(mnormt)
library(foreach)
## for all the functions below.

## ncore - number of cores to use - default to 1 so works for Windows

## th.p - theta parameters (vector)
## th0.p - theta0 parameters (vector)
## Sigma - 2x2 covariance matrix (matrix)
## Chol.M - Cholensky decomposition of hypothetical Sigma matrix (matrix)
## ------------------------------ ##
## xx - single x-coordinate (scalar)
## yy - single y-coordinate (scalar)
## tt - single time point (scalar)
## ------------------------------ ##
## xs - multiple x-coordinates (vector)
## ys - multiple y-coordinates (vector)
## ts - multiple times (vector)
## ------------------------------ ##
## xh - observed x-coordinates (vector)
## yh - observed y-coordinates (vector)
## th - observed times (vector)
## mh - observed magnitudes (vector)
## ------------------------------ ##
## M0 - magnitude of completeness (scalar)
## T1 - lower time interval extreme (scalar)
## T2 - upper time interval extreme (scalar)
## bdy - Polygon representing the study region (SpatialPolygon)

# magnitude triggering function (vector mh)
g.m <- function(th.p, mh, M0){
  #expm1( ( expm1(th.p[3]) + 1 ) * (mh - M0) ) + 1
  exp(exp(th.p[3])*(mh - M0))
}

# time triggering function (scalar tt and vector th)
g.t <- function(th.p, tt, th){
  output <- rep(0, length(th))
  idx.0 <- th >= tt
  
  gamma.h <- (tt - th[!idx.0])/exp(th.p[4])
  output[!idx.0] <- expm1( (-1-exp(th.p[5]))*log(gamma.h + 1) ) + 1
  output
}


## FUNCTIONS TO GET SIGMA FROM THETA BASED ON BARTLETT'S DECOMPOSITION
# function to get the values of the diagonal of A
#   mandatory input : th --> value of parameter theta 
#
#   optional input : df.p --> degrees of freedom (default set to 1)
# 
#   output : scalar, value of a, should be a Chi-square with df.p degrees of freedom 
# 
#   notes : it assumes th comes from a normal 0, 1 distribution
get_a <- function(th, df.p = 1) {
  bru_forward_transformation(qchisq, th, df.p)
}

# function to create Sigma from values of theta applying the Bartlett's decomposition
# function to get the values of the diagonal of A
#   mandatory input : theta6/7/8 --> value of parameter theta 
#                   : Chol.M --> Cholensky decomposition of a covariance matrix Sigma' 
# 
#   output : 2x2 Covariance matrix 
# 
#   notes : it assumes theta6/7/8 comes from a normal 0, 1 distribution.
#         : Chol.M is the Cholensky decomposition of Sigma' where Sigma' has to be seen as a prior 
#           covariance matrix

Sigma.from.theta <- function(th.c1, th.c2, th.n21, 
                             prior.sigma = list(Xbar = diag(c(1,1)), d = 2)){
  d_ = prior.sigma$d
  if(d_<= 1){
    stop('n_ parameter too low')
  }
  V_ <- prior.sigma$Xbar/d_
  sigma.1 = sqrt(V_[1,1])
  sigma.2 = sqrt(V_[2,2])
  ro = V_[1,2]/(sigma.1*sigma.2)
  
  c1 = sqrt(get_a(th.c1, df.p = d_ - 1 + 1))
  c2 = sqrt(get_a(th.c2, df.p = d_ - 2 + 1))
  n21 = th.n21
  
  output.11 <- (c1*sigma.1)^2
  output.22 <- (c1*sigma.2*ro + n21*sigma.2*sqrt(1 - ro^2))^2 + (1-ro^2)*((c2*sigma.2)^2)
  output.21 <- sigma.1*c1*(c1*sigma.2*ro + n21*sigma.2*sqrt(1 - ro^2))
  
  matrix(c(output.11, output.21, output.21, output.22), byrow = TRUE, ncol = 2)
}

theta.from.Sigma <- function(Sigma_, prior.sigma = list(Xbar = diag(c(1,1)), d = 2)){
  d_ <- prior.sigma$d
  V_ <- prior.sigma$Xbar/d_
  sigma.1 <- sqrt(V_[1,1])
  sigma.2 <- sqrt(V_[2,2])
  ro <- V_[1,2]/(sigma.1*sigma.2)
  
  c1 <- sqrt(Sigma_[1,1])/sigma.1
  n21 <- (Sigma_[2,1]/(sigma.1*sigma.2) - ro*c1^2)/(c1*sqrt(1-ro^2)) 
  c2 <- sqrt( (Sigma_[2,2]/(sigma.2^2) - (n21*sqrt(1-ro^2) + ro*c1)^2)/(1-ro^2) )
  c(th6 = qnorm(pchisq(c1^2, d_ - 1 + 1, log.p = TRUE), log.p = TRUE),
    th7 = qnorm(pchisq(c2^2, d_ - 2 + 1, log.p = TRUE), log.p = TRUE),
    th8 = n21)
}

get_ro <- function(th.ro, df){
  inlabru:::bru_forward_transformation(qbeta, th.ro, shape1 = (df - 1)/2,
                                       shape2 = (df - 1)/2)*(2 - 2e-6) - (1 - 1e-6)  
}


Sigma.from.theta2 <- function(th.s1, th.ro, df){
  sigma.1 <- exp(th.s1) 
  sigma.2 <- sigma.1
  ro <- get_ro(th.ro, df) 
  matrix(c(sigma.1^2, ro*sigma.1*sigma.2, ro*sigma.1*sigma.2, sigma.2^2), byrow = TRUE, ncol = 2)
}


g.s <- function(th.p, xx, yy, xh, yh, Sigma){
  if(length(xh) != length(yh)){
    stop('xh and yh has to be of the same length')
  }
  
  # mean.h <- cbind(xh, yh)
  # unlist(mclapply(1:nrow(mean.h), \(idx) dmvnorm(c(xx,yy), mean = mean.h[idx,], sigma = Sigma),
  #                 mc.cores = 3))
  mean.dif <- cbind(xx -xh, yy - yh)
  S.det <- Sigma[1,1]*Sigma[2,2] - Sigma[1,2]*Sigma[2,1] 
  S.inv <- (1/S.det)*matrix(c(Sigma[2,2], -Sigma[1,2], -Sigma[1,2], Sigma[1,1]), 
                            byrow = TRUE, ncol = 2)
  print(S.det)
  log.out <- -log(2*pi) -(1/2)*log(S.det) -(1/2)*as.numeric(diag( mean.dif %*% S.inv %*% t(mean.dif) )) 
  exp(log.out)
}


# triggering function.
g.x <- function(th.p, tt, xx, yy, th, xh, yh, mh, M0, Sigma){
  log.output <- th.p[2] + log1p(g.m(th.p, mh, M0)-1) + log1p(g.t(th.p, tt, th)-1) + log1p(g.s(th.p, xx, yy, xh, yh, Sigma)-1)
  expm1(log.output) + 1
}

## conditional intensity
lambda <- function(th.p, tt, xx, yy, th, xh, yh, mh, M0, Sigma){
  expm1(th.p[1]) + 1 + sum(g.x(th.p, tt, xx, yy, th, xh, yh, mh, M0, Sigma))
}


## integrated intensity for single point
logLambda.h.vec <- function(th.p, th, xh, yh, mh, M0, T1, T2, bdy, Sigma){
  
  if(any(th > T2)){
    return('Error: th > T2')
  }
  
  # integral in space
  lower_ <- c(bdy@bbox[1,1], bdy@bbox[2,1])
  upper_ <- c(bdy@bbox[1,2], bdy@bbox[2,2])
  p.loc <- cbind(xh, yh)
  I.s <- sapply(1:nrow(p.loc), \(idx)
                abs(pmvnorm(lower = lower_, upper = upper_, 
                            mean = p.loc[idx,], sigma = Sigma,
                            keepAttr = FALSE)))
  
  # integral in time
  Tl <- sapply(th, \(x) max(x, T1))
  gamma.l <- (Tl - th)/exp(th.p[4])
  gamma.u <- (T2 - th)/exp(th.p[4])
  w.l <- (gamma.l + 1)^(-exp(th.p[5]))  
  w.u <- (gamma.u + 1)^(-exp(th.p[5]))
  
  # output
  th.p[2] + exp(th.p[3])*(mh - M0) + th.p[4] - th.p[5] + log1p(w.l - 1) + log1p(-w.u/w.l) + log1p(I.s - 1)
}


ETAS.loglik <- function(th.p, th, xh, yh, mh, M0, T1, T2, bdy, Sigma, 
                        prec.prior = rep(1,5), ncore = 1){
  
  sd.prior <- sqrt(1/prec.prior)
  logIntegral.trig <- logLambda.h.vec(th.p, th, xh, yh, mh, M0, T1, T2, bdy, Sigma)
  Integral. <- exp(th.p[1])*(T2 - T1) + exp(logSumExp(logIntegral.trig))
  log.lambdas <- unlist(mclapply(1:length(th), \(ii) 
                                 log1p(lambda(th.p, th[ii], xh[ii], yh[ii],
                                              th, xh, yh, mh, M0, Sigma) - 1), 
                                 mc.cores = ncore))
  
  - Integral. + sum(log.lambdas) +
    sum(dnorm(th.p[1:5], sd = sd.prior, log = TRUE))
}


## integrated time-triggering function
It <- function(th.p, th, T2){
  gamma.u <- (T2 - th)/exp(th.p[4])
  exp(th.p[4] - th.p[5])*(1 - (gamma.u + 1)^(-exp(th.p[5])))
}
## inverse of integrated time-triggering function
Inv.It <- function(th.p, omega, th){
  th + exp(th.p[4])*( (1 - omega*exp(-th.p[4] + th.p[5]))^(-1/exp(th.p[5])) - 1)
}

## sampling times
sample.time <- function(th.p, n.ev, th, T2){
  if(n.ev == 0){
    df <- data.frame(ts = 1, x = 1, y = 1, mags = 1, gen = 0)
    return(df[-1,])
  }
  bound.l <- 0 #It(th.p, th, T)
  bound.u <- It(th.p, th, T2)
  unif.s <- runif(n.ev, min = bound.l, max = bound.u)
  t.sample <- Inv.It(th.p, unif.s, th)
  t.sample
}

# sampling locations
sample.loc <- function(th.p, xh, yh, n.ev, bdy, Sigma){
  
  # initialize empty SpatialPoints
  samp.points <- SpatialPoints(data.frame(x = 0, y = 0))
  samp.points <- samp.points[-1,]
  # initialize number of placed events
  num <- 0
  
  # until we placed all the events
  while (num < n.ev){
    # sample from the 2D Gaussian kernel without boundaries
    pts.matrix <- rmvnorm(n.ev, mean = c(xh, yh), sigma = Sigma)
    # transform the sample in SpatialPoints and project
    pts <- data.frame(x = pts.matrix[,1], y = pts.matrix[,2])
    coordinates(pts) <- ~ x + y
    # discard the ones outside bdy
    pts <- crop(pts, bdy)
    # merge sampled points
    if(length(pts) > 0){
      samp.points <- rbind(samp.points, pts)
      num <- length(samp.points)
    }
    else{
      #print('no retained')
    }
  }
  
  ## if we retained a number of points > n.ev, we select n.ev events at random
  kp <- sample(seq(1, length(samp.points), by=1), n.ev, replace=FALSE)
  samp.points <- samp.points[kp,]
  samp.points
}


# sampling triggered events
sample.triggered <- function(th.p, beta.p, th, xh, yh, n.ev, M0, T1, T2, bdy, 
                             Sigma){
  # if the number of events to be placed is zero returns an empty data.frame
  if(n.ev == 0){
    samp.points <- data.frame(x = 1, y = 1, ts = 1, mags = 1)
    samp.points <- samp.points[-1,]
    return(samp.points)
  }
  else{
    
    # sample times
    samp.ts <- sample.time(th.p, n.ev, th, T2)
    # sample magnitudes
    samp.mags <- rexp(n.ev, rate = beta.p) + M0
    # sample locations
    samp.locs <- sample.loc(th.p = th.p, xh = xh, yh = yh, n.ev = n.ev, 
                            bdy = bdy, Sigma = Sigma)
    
    # build output dataset
    samp.points <- data.frame(ts = samp.ts, mags = samp.mags, x = samp.locs@coords[,1],
                              y = samp.locs@coords[,2])
    # return only the ones with time different from NA (the one with NA are outside the interval T1, T2)
    # even though it should not happen given how we built sample.omori
    return(samp.points[!is.na(samp.points$ts),])
  }
}


sample.generation <- function(th.p, beta.p, Ht, M0, T1, T2, bdy,
                              Sigma, ncore = 1){
  
  # number of parents
  n.parent <- nrow(Ht)
  # calculate the aftershock rate for each parent in history
  trig.rates <- exp(logLambda.h.vec(th.p, Ht$ts, Ht$x, Ht$y, Ht$mags, M0, T1, T2, bdy,
                                    Sigma))
  # extract number of aftershock for each parent
  n.ev.v <- sapply(1:n.parent, function(x) rpois(1, trig.rates[x]))
  
  # if no aftershock has to be generated returns empty data.frame
  if(sum(n.ev.v) == 0){
    app <- data.frame(x = 1, y = 1, ts = 1, mags = 1)
    app <- app[-1,]
    return(app)
  }
  
  # identify parent with number of aftershocks > 0 
  idx.p <- which(n.ev.v > 0)
  
  #print(sample.triggered(theta.v, beta.p, Sigma, Chol.M, n.ev.v[idx.p[1]], Ht[idx.p[1],], T1, T2, M0, bdy, crsobj))
  # sample (in parallel) the aftershocks for each parent 
  sample.list <- mclapply(idx.p, function(idx) 
    sample.triggered(th.p, beta.p, Ht$ts[idx], Ht$x[idx], Ht$y[idx], n.ev.v[idx], M0, T1, T2, bdy, 
                     Sigma), mc.cores = ncore)
  
  # bind the data.frame in the list and return
  sample.pts <- bind_rows(sample.list) 
  sample.pts
}


point_sampler <- function(loglambda, bdy, mesh, num_events, 
                          crs_obj = "+proj=longlat +datum=WGS84 +no_defs"){
  
  if(is.na(proj4string(bdy))){
     proj4string(bdy) <- crs_obj
     bdy <- spTransform(bdy, crs_obj)
   }
  ## Number of events for single catalogue from a poisson distribution with lambda = num_events
  lambda_max <- max(loglambda)
  ## number of points to sample at a time - might want to adjust depending on how many points you want to actually retain.
  n.points=1000
  
  ## Set up a spatialpoints dataframe for our results
  samp.points <- SpatialPoints(data.frame(x = 0, y = 0))
  samp.points <- samp.points[-1,]
  proj4string(samp.points) <- crs_obj
  num <- 0
  
  
  ## To sample the correct number of points, keep going until the num >= cat_n
  while (num < num_events){
    pts <- spsample(bdy, n.points, "random")#points <- spsample(bdy, n.points, "random")
    ## transform to wgs84 long/lat
    pts <- spTransform(pts, crs_obj)
    ## next few lines modified directly from sample.lgcp
    proj <- INLA::inla.mesh.project(mesh, pts)
    lambda_ratio <- exp(as.vector(proj$A %*% loglambda) - lambda_max)
    keep <- proj$ok & (runif(n.points) <= lambda_ratio)
    if(sum(keep) == 0){
      print('zero keeped')
      next
    }
    kept <- pts[keep]
    
    samp.points <- rbind(samp.points, kept)
    num <- length(samp.points)
    
  }
  
  ## Keep exactly cat_n points, choose these randomly from all of the points we've kept so far
  kp <- sample(seq(1, length(samp.points), by=1), num_events, replace=FALSE)
  samp.points <- samp.points[kp,]
  
  return(samp.points)
}


sample.ETAS <- function(th.p, beta.p, M0, T1, T2, bdy, Sigma, 
                        loglambda.bkg = NULL, mesh.bkg = NULL, 
                        crs_obj = NULL, Ht = NULL, ncore = 1,
                        Unit.Area = TRUE){
  # if the upper extreme greater than lower
  if(T2 < T1){
    stop('Error - right-end of time interval greater than left-end')
  }
  
  if(Unit.Area){
    n.bkg <- rpois(1, exp(th.p[1])*(T2 - T1))
  }
  else{
    bdy.sf <- st_as_sf(bdy.sf)
    st_crs(bdy.sf) <- italy.crs
    bdy.sf <- st_transform(bdy.sf, crs_obj)
    Area.bdy <- as.numeric(st_area(bdy.sf)/(1000^2))
    n.bkg <- rpois(1, exp(th.p[1])*(T2 - T1)*Area.bdy)
  }
  
  #print(exp(th.p[1])*(T2 - T1))
  # if no background events are generated initialize an empty data.frame
  if(n.bkg == 0){
    bkg.df <- data.frame(x = 1, y = 1, ts = 1, mags = 1, gen = 0)
    bkg.df <- bkg.df[-1,]
  }
  else{
    # sample bkg events
    # if no bk.field.list element is passed it assumes uniform background rate
    if(is.null(loglambda.bkg)){
      bkg.locs <- spsample(bdy, n.bkg, 'random', iter = 10)
      # if(!is.null(crsobj)){
      #   proj4string(bkg.locs) <- crs_obj
      #   bkg.locs <- spTransform(bkg.locs, crs_obj)
      # }
      bkg.df <- data.frame(x = bkg.locs@coords[,1], 
                           y = bkg.locs@coords[,2], 
                           ts = runif(n.bkg, T1, T2), 
                           mags = rexp(n.bkg, beta.p) + M0, 
                           gen = 1)
    }
    # otherwise it samples using the information provided
    else{
      
      bkg.locs <- suppressWarnings(point_sampler(loglambda = loglambda.bkg, bdy = bdy, 
                                mesh = mesh.bkg, num_events = n.bkg,
                                crs_obj = crs_obj))
      
      bkg.df <- data.frame(x = data.frame(bkg.locs)$x, 
                           y = data.frame(bkg.locs)$y, 
                           ts = runif(n.bkg, T1, T2), 
                           mags = rexp(n.bkg, beta.p) + M0, 
                           gen = 1)
    }
  }
  
  # if known events are provided
  if(!is.null(Ht)){
    #### TODO : the function has to generate all the points.
    # sample a generation from the known events
    gen.from.past <- sample.generation(th.p, beta.p, Ht, M0, T1, T2,  bdy, Sigma, ncore)
    # if at least an aftershock is produced
    if(nrow(gen.from.past) > 0){
      # set generation
      gen.from.past$gen = 0
      # Merge first generation and background events
      Gen.list <- list(rbind(gen.from.past, bkg.df))  
    }
    else{
      Gen.list <- list(bkg.df)
    }
    
  }
  else{
    Gen.list <- list(bkg.df)
  }
  # stop if we have no background events and no events generated from known observations
  if(nrow(Gen.list[[1]]) == 0){
    #print(exp(theta.v[1])*(T2 - T1)*(area(bdy)/1000000))
    #stop('No events generated - increase theta1')
    return(Gen.list)
  }
  
  # initialize flag and gen counter
  flag = TRUE
  gen = 1
  # this goes until the condition inside the loop is met
  while(flag){
    # set parents
    parents <- Gen.list[[gen]]
    #print(c(T1,T2))
    #print(range(parents$ts))
    # generate aftershocks
    triggered <- sample.generation(th.p, beta.p, parents, 
                                   M0, T1, T2, bdy, Sigma, ncore)
    #print(nrow(triggered))
    # stop the loop if there are no more aftershocks
    if(nrow(triggered) == 0){
      flag = FALSE}
    else{
      # set generations
      triggered$gen = gen + 1
      # store new generation
      Gen.list[[gen + 1]] = triggered
      # update generation counter
      gen = gen + 1
    }
  }
  Gen.list
}


### function for model fitting

### functions for model fitting
ETAS.fit.B_bkg <- function(sample.s, N.breaks.min, max.length, M0, T1, T2, bdy, 
                       Sigma = NULL,
                       prior.mean = rep(0,5), 
                       prior.prec = rep(1,5),
                       bru.opt = list(bru_verbose = 3,
                                      bru_max_iter = 50),
                       bin.strat = 'def',
                       ncore = 1){
  
  # create different poisson counts models
  # first for background
  df.0 <- data.frame(counts = 0, exposures = 1)
  # this is the expression of log(Lambda0)
  form.0 <- counts ~ th.mu + log(T2 - T1) #+ log(Area.bdy)
  
  # first likelihood
  lik.0 <- inlabru::like(formula = form.0,
                         data = df.0,
                         family = 'poisson',
                         options = list(E = df.0$exposures))
  
  
  ### create time bins 
  
  df.j <- foreach(idx = 1:nrow(sample.s), .combine = rbind) %do% {
    
    tt <- sample.s$ts[idx]
    if(bin.strat == 'def'){
      N.breaks <- max(N.breaks.min, ceiling((T2 - tt)/max.length))
      kk <- seq_len(N.breaks) - 1
      Time.bins.matrix <- tt + cbind(kk, kk + 1) / N.breaks * (T2- tt)
    }
    
    else if(bin.strat == 'alt'){
      t.max <- min(tt + 10, T2)
      max.length_ <- (t.max - tt)/20
      N.breaks <- max(N.breaks.min, ceiling((t.max - tt)/max.length_))
      kk <- seq_len(N.breaks) - 1
      Time.bins.matrix <- tt + cbind(kk, kk + 1) / N.breaks * (t.max - tt)
      if(t.max != T2){
        Time.bins.matrix <- rbind(Time.bins.matrix, c(t.max, T2))  
      }
    }
    else{
      stop('Unknown bin strategy')
    }
    
    
    data.frame(x = rep(sample.s$x[idx], each = nrow(Time.bins.matrix)),
               y = rep(sample.s$y[idx], each = nrow(Time.bins.matrix)),
               ts = rep(sample.s$ts[idx], each = nrow(Time.bins.matrix)),
               mags = rep(sample.s$mags[idx], each = nrow(Time.bins.matrix)),
               counts = 0,
               exposures = 1,
               bin.start = Time.bins.matrix[,1],
               bin.end = Time.bins.matrix[,2])
    
  }
  
  
  logLambda.h.inla <- function(th.K, th.alpha, th.c, th.pm1, ts, x, y, mags, 
                               M0, T1.v, T2.v, bdy, Sigma, ncore_ = ncore){
    th.p <- c(0, th.K[1], th.alpha[1], th.c[1], th.pm1[1])
    unlist(mclapply(1:length(ts), \(idx)
                    logLambda.h.vec(th.p, ts[idx], x[idx], y[idx], mags[idx], M0, T1.v[idx], T2.v[idx], bdy, 
                                    Sigma),
                    mc.cores = ncore_))
  }
  # creating formula for past events contributions to integrated lambda
  form.j.part <- counts ~ logLambda.h.inla(th.K = th.K, 
                                           th.alpha = th.alpha, 
                                           th.c = th.c, 
                                           th.pm1 = th.pm1, 
                                           ts = ts, x = x, y = y, mags = mags, M0 = M0, 
                                           T1.v = bin.start, T2.v = bin.end, bdy = bdy, 
                                           Sigma = Sigma)
  # second for triggered part of the integral
  lik.j.part <- inlabru::like(formula = form.j.part,
                              data = df.j,
                              family = 'poisson',
                              options = list(E = df.j$exposures)) 
  
  # third is for the sum of the log intensities
  df.s <- data.frame(x = sample.s$x, y = sample.s$y, ts = sample.s$ts, mags = sample.s$mags, 
                     bkg = sample.s$bkg, counts = 1, exposures = 0)
  
  loglambda.inla <- function(th.mu, th.K, th.alpha, th.c, th.pm1, bkg, tt, xx, yy, 
                             th, xh, yh, mh, M0, Sigma, ncore_ = ncore){
    
    th.p.df <- data.frame(th.mu = th.mu[1] + bkg, 
                          th.K = th.K[1], 
                          th.alpha = th.alpha[1], 
                          th.c = th.c[1], 
                          th.pm1 = th.pm1[1])
    
    unlist(mclapply(1:length(tt), \(idx) 
                    log(lambda(as.numeric(th.p.df[idx,]), tt[idx], xx[idx], yy[idx], 
                               th, xh, yh, mh, M0, Sigma)),
                    mc.cores = ncore_))
  }
  # creating formula for summation part
  form.s.part <- counts ~ loglambda.inla(th.mu = th.mu, 
                                         th.K = th.K, 
                                         th.alpha = th.alpha, 
                                         th.c = th.c, 
                                         th.pm1 = th.pm1,
                                         bkg = bkg,
                                         tt = sample.s$ts, xx = sample.s$x, yy = sample.s$y,
                                         th = sample.s$ts, xh = sample.s$x, yh = sample.s$y, 
                                         mh = sample.s$mags, M0 = M0, 
                                         Sigma = Sigma)
  lik.s.part <- inlabru::like(formula = form.s.part,
                              data = df.s,
                              family = 'poisson',
                              options = list(E = df.s$exposures)) 
  cmp.part <- counts ~ -1 + 
    th.mu(1, model = 'linear', mean.linear = prior.mean[1], prec.linear = prior.prec[1]) + 
    th.K(1, model = 'linear', mean.linear = prior.mean[2], prec.linear = prior.prec[2]) +
    th.alpha(1, model = 'linear', mean.linear = prior.mean[3], prec.linear = prior.prec[3]) +
    th.c(1, model = 'linear', mean.linear = prior.mean[4], prec.linear = prior.prec[4]) +
    th.pm1(1, model = 'linear', mean.linear = prior.mean[5], prec.linear = prior.prec[5]) 
  
  bru(lik.0, lik.j.part, lik.s.part, components = cmp.part, 
      options = bru.opt)
}




### functions for model fitting
ETAS.fit.B_bkg_full <- function(sample.s, N.breaks.min, max.length, M0, T1, T2, bdy,
                               prior.sigma = list(Xbar = diag(c(1,1)), d = 2),
                               prior.mean = rep(0,5), 
                               prior.prec = rep(1,5),
                               bru.opt = list(bru_verbose = 3,
                                          bru_max_iter = 50),
                               bin.strat = 'def'){
  
  # create different poisson counts models
  # first for background
  df.0 <- data.frame(counts = 0, exposures = 1)
  # this is the expression of log(Lambda0)
  form.0 <- counts ~ th.mu + log(T2 - T1) 
  
  # first likelihood
  lik.0 <- inlabru::like(formula = form.0,
                         data = df.0,
                         family = 'poisson',
                         options = list(E = df.0$exposures))
  
  
  ### create time bins 
  
  df.j <- foreach(idx = 1:nrow(sample.s), .combine = rbind) %do% {
    
    tt <- sample.s$ts[idx]
    if(bin.strat == 'def'){
      N.breaks <- max(N.breaks.min, ceiling((T2 - tt)/max.length))
      kk <- seq_len(N.breaks) - 1
      Time.bins.matrix <- tt + cbind(kk, kk + 1) / N.breaks * (T2- tt)
    }
    
    else if(bin.strat == 'alt'){
      t.max <- min(tt + 10, T2)
      max.length_ <- (t.max - tt)/20
      N.breaks <- max(N.breaks.min, ceiling((t.max - tt)/max.length_))
      kk <- seq_len(N.breaks) - 1
      Time.bins.matrix <- tt + cbind(kk, kk + 1) / N.breaks * (t.max - tt)
      if(t.max != T2){
        Time.bins.matrix <- rbind(Time.bins.matrix, c(t.max, T2))  
      }
    }
    else{
      stop('Unknown bin strategy')
    }
    
    
    data.frame(x = rep(sample.s$x[idx], each = nrow(Time.bins.matrix)),
               y = rep(sample.s$y[idx], each = nrow(Time.bins.matrix)),
               ts = rep(sample.s$ts[idx], each = nrow(Time.bins.matrix)),
               mags = rep(sample.s$mags[idx], each = nrow(Time.bins.matrix)),
               counts = 0,
               exposures = 1,
               bin.start = Time.bins.matrix[,1],
               bin.end = Time.bins.matrix[,2])
    
  }
  
  
  logLambda.h.inla <- function(th.K, th.alpha, th.c, th.pm1, th.c1, th.c2, th.n21, ts, x, y, mags, 
                               M0, T1.v, T2.v, bdy, p.sigma = prior.sigma, ncore_ = ncore){
    th.p <- c(0, th.K[1], th.alpha[1], th.c[1], th.pm1[1])
    Sigma_ <- Sigma.from.theta(th.c1[1], th.c2[1], th.n21[1], p.sigma)
    unlist(mclapply(1:length(ts), \(idx)
                    logLambda.h.vec(th.p, ts[idx], x[idx], y[idx], mags[idx], M0, T1.v[idx], T2.v[idx], bdy, 
                                    Sigma_),
                    mc.cores = ncore_))
  }
  # creating formula for past events contributions to integrated lambda
  form.j.part <- counts ~ logLambda.h.inla(th.K = th.K, 
                                           th.alpha = th.alpha, 
                                           th.c = th.c, 
                                           th.pm1 = th.pm1, 
                                           th.c1 = th.c1,
                                           th.c2 = th.c2,
                                           th.n21 = th.n21,
                                           ts = ts, x = x, y = y, mags = mags, M0 = M0, 
                                           T1.v = bin.start, T2.v = bin.end, bdy = bdy)
  # second for triggered part of the integral
  lik.j.part <- inlabru::like(formula = form.j.part,
                              data = df.j,
                              family = 'poisson',
                              options = list(E = df.j$exposures)) 
  
  # third is for the sum of the log intensities
  df.s <- data.frame(x = sample.s$x, y = sample.s$y, ts = sample.s$ts, mags = sample.s$mags, 
                     bkg = sample.s$bkg, counts = 1, exposures = 0)
  
  loglambda.inla <- function(th.mu, th.K, th.alpha, th.c, th.pm1, th.c1, th.c2, th.n21, bkg, tt, xx, yy, 
                             th, xh, yh, mh, M0, p.sigma = prior.sigma, ncore_ = ncore){
    
    th.p.df <- data.frame(th.mu = th.mu[1] + bkg, 
                          th.K = th.K[1], 
                          th.alpha = th.alpha[1], 
                          th.c = th.c[1], 
                          th.pm1 = th.pm1[1])
    Sigma_ <- Sigma.from.theta(th.c1[1], th.c2[1], th.n21[1], p.sigma)
    
    unlist(mclapply(1:length(tt), \(idx) 
                    log(lambda(as.numeric(th.p.df[idx,]), tt[idx], xx[idx], yy[idx], 
                               th, xh, yh, mh, M0, Sigma_)),
                    mc.cores = ncore_))
  }
  
  # creating formula for summation part
  form.s.part <- counts ~ loglambda.inla(th.mu = th.mu, 
                                         th.K = th.K, 
                                         th.alpha = th.alpha, 
                                         th.c = th.c, 
                                         th.pm1 = th.pm1,
                                         th.c1 = th.c1,
                                         th.c2 = th.c2,
                                         th.n21 = th.n21,
                                         bkg = bkg,
                                         tt = sample.s$ts, xx = sample.s$x, yy = sample.s$y,
                                         th = sample.s$ts, xh = sample.s$x, yh = sample.s$y, 
                                         mh = sample.s$mags, M0 = M0)
  lik.s.part <- inlabru::like(formula = form.s.part,
                              data = df.s,
                              family = 'poisson',
                              options = list(E = df.s$exposures)) 
  cmp.part <- counts ~ -1 + 
    th.mu(1, model = 'linear', mean.linear = prior.mean[1], prec.linear = prior.prec[1]) + 
    th.K(1, model = 'linear', mean.linear = prior.mean[2], prec.linear = prior.prec[2]) +
    th.alpha(1, model = 'linear', mean.linear = prior.mean[3], prec.linear = prior.prec[3]) +
    th.c(1, model = 'linear', mean.linear = prior.mean[4], prec.linear = prior.prec[4]) +
    th.pm1(1, model = 'linear', mean.linear = prior.mean[5], prec.linear = prior.prec[5]) + 
    th.c1(1, model = 'linear', mean.linear = 0, prec.linear = 1) +
    th.c2(1, model = 'linear', mean.linear = 0, prec.linear = 1) +
    th.n21(1, model = 'linear', mean.linear = 0, prec.linear = 1) 
  
  bru(lik.0, lik.j.part, lik.s.part, components = cmp.part, 
      options = bru.opt)
}

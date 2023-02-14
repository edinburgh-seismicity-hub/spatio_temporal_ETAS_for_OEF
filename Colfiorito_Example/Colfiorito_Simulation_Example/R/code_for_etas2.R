library(sf)
library(inlabru)
library(INLA)
library(dplyr)
library(matrixStats)
library(mvtnorm)
library(MASS)
library(raster)
library(rgeos)
library(mnormt)
library(foreach)
library(viridis)


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

get_ro <- function(th.ro, df){
  inlabru:::bru_forward_transformation(qbeta, th.ro, shape1 = (df - 1)/2,
                                       shape2 = (df - 1)/2)*(2 - 2e-6) - (1 - 1e-6)  
}


Sigma.from.theta <- function(th.s1, th.ro, df){
  sigma.1 <- exp(th.s1) 
  sigma.2 <- sigma.1
  ro <- get_ro(th.ro, df) 
  matrix(c(sigma.1^2, ro*sigma.1*sigma.2, ro*sigma.1*sigma.2, sigma.2^2), byrow = TRUE, ncol = 2)
}


g.s <- function(xx, yy, xh, yh, Sigma){
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
  
  log.out <- -log(2*pi) -(1/2)*log(S.det) -(1/2)*as.numeric(diag( mean.dif %*% S.inv %*% t(mean.dif) )) 
  exp(log.out)
}


# triggering function.
g.x <- function(theta, tt, xx, yy, th, xh, yh, mh, M0, Sigma){
  # initialize output
  output <- rep(0, length(th))
  t.diff <- tt - th
  upd <- t.diff > 0
  if(sum(upd) > 0){
    output[upd] <- theta[2]*exp(theta[3]*(mh[upd] - M0))*( (1 + t.diff[upd]/theta[4])^(-theta[5]) )*g.s(xx,yy,xh[upd],yh[upd],Sigma)
    output  
  }
  else{
    output
  }
}


## conditional intensity
lambda <- function(theta, tt, xx, yy, th, xh, yh, mh, M0, Sigma){
  theta[1]  + sum(g.x(theta, tt, xx, yy, th, xh, yh, mh, M0, Sigma))
}




## integrated intensity for single point
logLambda.h.vec <- function(theta, th, xh, yh, mh, M0, T1, T2, bdy, Sigma){
  
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
  Tl <- sapply(th, \(x) max(as.numeric(x), as.numeric(T1)))
  
  gamma.T2 <- ( 1 + (T2 - th)/theta[4] )^(1-theta[5])  
  gamma.Tl <- ( 1 + (Tl - th)/theta[4] )^(1-theta[5])
  
  # output
  log(theta[2]) + theta[3]*(mh - M0) + log(theta[4]) - log(theta[5] - 1) + log(gamma.Tl - gamma.T2) + log(I.s)
}



# cross-reference paper section

# length first bin
breaks_exp <- function(tt_, T2_, coef_ = 2, max_length_, N_exp_ = 10){
  
  tt_breaks <- tt_ + max_length_*((1 + coef_)^(0:N_exp_))
  tt_breaks <- tt_breaks[tt_breaks < T2]
  if(T2 - tt_breaks[length(tt_breaks)] < max_length_){
    tt_breaks[length(tt_breaks)] = T2_
  }
  if(T2_ - tt_ < max_length_){
    c(tt_, T2_)
  }
  if(tt_breaks[length(tt_breaks)] < T2_){
    tt_breaks <- c(tt_breaks, T2_)
  }
  c(tt_,tt_breaks)
} 


### functions for model fitting
ETAS.fit.B_bkg <- function(sample.s, 
                           M0, T1, T2, bdy, 
                           Sigma = NULL,
                           N.breaks.min = 1, max.length = 1, 
                           N.tail = 1, t.jump = 1,
                           coef_exp = 2, delta_exp = 0.1,
                           N_exp = 50,
                           min.tail.length = 0.05,
                           link.fun = list(mu = \(x) x,
                                           K = \(x) x,
                                           alpha = \(x) x,
                                           cc = \(x) x,
                                           pp = \(x) x), 
                           bru.opt = list(bru_verbose = 3,
                                          bru_max_iter = 50),
                           bin.strat = 'def',
                           ncore = 1){
  
  # create different poisson counts models
  # first for background
  df.0 <- data.frame(counts = 0, exposures = 1)
  # this is the expression of log(Lambda0)
  form.0 <- counts ~ link.fun$mu(th.mu) + log(T2 - T1) #+ log(Area.bdy)
  
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
      t.max <- min(tt + t.jump, T2)
      N.breaks <- max(N.breaks.min, ceiling((t.max - tt)/max.length))
      kk <- seq_len(N.breaks) - 1
      Time.bins.matrix <- tt + cbind(kk, kk + 1) / N.breaks * (t.max - tt)
      if(t.max != T2){
        N.tail.max <- max( ceiling((T2 - t.max)/min.tail.length) , 2) 
        N.tail <- min(N.tail, N.tail.max)
        kk.tail <- seq_len(N.tail) - 1
        tail.bins.matrix <- t.max + cbind(kk.tail, kk.tail + 1) / N.tail * (T2 - t.max)
        
        Time.bins.matrix <- rbind(Time.bins.matrix, tail.bins.matrix)  
      }
    }
    else if(bin.strat == 'exp'){
      t_b <- breaks_exp(tt, T2, coef_exp, delta_exp, N_exp)
      Time.bins.matrix <- cbind(t_b[-length(t_b)], t_b[-1])
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
  print(df.j)
  logLambda.h.inla <- function(th.K, th.alpha, th.c, th.p, ts, x, y, mags, 
                               M0, T1.v, T2.v, bdy, Sigma, ncore_ = ncore){
    theta_ <- c(0, 
              link.fun$K(th.K[1]), 
              link.fun$alpha(th.alpha[1]), 
              link.fun$cc(th.c[1]), 
              link.fun$pp(th.p[1]))
    #cat('theta - LogL', theta_, '\n')
    outp<- unlist(mclapply(1:length(ts), \(idx)
                    logLambda.h.vec(theta = theta_, 
                                    th = ts[idx], 
                                    xh = x[idx], 
                                    yh = y[idx], 
                                    mh = mags[idx],
                                    M0 = M0, 
                                    T1 = T1.v[idx], 
                                    T2 = T2.v[idx], 
                                    bdy = bdy,
                                    Sigma = Sigma),
                    mc.cores = ncore_))
    #cat('logL', outp, '\n')
    outp
    }
  
  
  # creating formula for past events contributions to integrated lambda
  form.j.part <- counts ~ logLambda.h.inla(th.K = th.K, 
                                           th.alpha = th.alpha, 
                                           th.c = th.c, 
                                           th.p = th.p, 
                                           ts = ts, x = x, y = y, 
                                           mags = mags, M0 = M0, 
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
  
  loglambda.inla <- function(th.mu, th.K, th.alpha, th.c, th.p, bkg, tt, xx, yy, 
                             th, xh, yh, mh, M0, Sigma, ncore_ = ncore){
    
    
    th.p.df <- data.frame(th.mu = link.fun$mu(th.mu[1])*bkg, 
                          th.K = link.fun$K(th.K[1]), 
                          th.alpha = link.fun$alpha(th.alpha[1]), 
                          th.c = link.fun$cc(th.c[1]), 
                          th.p = link.fun$pp(th.p[1]))
    #cat('th.mu', sum(th.p.df$th.mu <0), '\n')
    
    outp <- unlist(mclapply(1:length(tt), \(idx) 
                    log(lambda(theta = as.numeric(th.p.df[idx,]), 
                               tt = tt[idx], 
                               xx = xx[idx], 
                               yy = yy[idx], 
                               th = th, xh = xh, yh = yh, mh = mh, M0 = M0, Sigma = Sigma)),
                    mc.cores = ncore_))
    #cat('logl', sum(is.na(outp)), '\n')
    outp
  }
  # creating formula for summation part
  form.s.part <- counts ~ loglambda.inla(th.mu = th.mu, 
                                         th.K = th.K, 
                                         th.alpha = th.alpha, 
                                         th.c = th.c, 
                                         th.p = th.p,
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
    th.mu(1, model = 'linear', mean.linear = 0, prec.linear = 1) + 
    th.K(1, model = 'linear', mean.linear = 0, prec.linear = 1) +
    th.alpha(1, model = 'linear', mean.linear = 0, prec.linear = 1) +
    th.c(1, model = 'linear', mean.linear = 0, prec.linear = 1) +
    th.p(1, model = 'linear', mean.linear = 0, prec.linear = 1) 
  
   bru(lik.0, lik.j.part, lik.s.part, components = cmp.part, 
        options = bru.opt)
   
  
}

Plot_grid <- function(xx, yy, delta_, n.layer, bdy_,
                      min.edge = 1, displaypl = TRUE){
  # extract boundaries
  x.min <- bdy_@bbox[1,1]
  x.max <- bdy_@bbox[1,2]
  y.min <- bdy_@bbox[2,1]
  y.max <- bdy_@bbox[2,2]
  
  n.layer <- floor(min(c(n.layer, 
                         abs(x.min - xx)/delta_,
                         abs(x.max - xx)/delta_, 
                         abs(y.min - yy)/delta_,
                         abs(y.max - yy)/delta_)))
  coords.list <- list()
  poly.list <- list()
  cat('Number of layers : ', n.layer, '\n')
  poly.name <- c()
  red_coef = 1
  while(n.layer > 0 & 
        (abs(x.min - (xx - n.layer*delta_)) <= min.edge |
         abs(x.max - (n.layer*delta_ + xx)) <= min.edge |
         abs(y.min - (yy - n.layer*delta_)) <= min.edge |
         abs(y.max - (n.layer*delta_ + yy))  <= min.edge ) ){
    n.layer = n.layer - red_coef
    cat('Reducing number of layers by ', red_coef, '\n')
    red_coef <- red_coef + 1
  }
  if(n.layer > 0){
    coeff_ = 1
    for(layer_ in 1:n.layer){
      coords.1 <- cbind(x = c(xx - coeff_*delta_, xx - coeff_*delta_, xx + (coeff_ - 1)*delta_, xx + (coeff_ - 1)*delta_, xx - coeff_*delta_),
                        y = c(yy + (coeff_ - 1)*delta_, yy + coeff_*delta_, yy + coeff_*delta_, yy + (coeff_ - 1)*delta_, yy + (coeff_ - 1)*delta_))
      coords.2 <- cbind(x = c(xx + (coeff_ - 1)*delta_, xx + (coeff_ - 1)*delta_, xx + coeff_*delta_, xx + coeff_*delta_, xx + (coeff_ - 1)*delta_),
                        y = c(yy - (coeff_ - 1)*delta_, yy + coeff_*delta_, yy + coeff_*delta_, yy - (coeff_ - 1)*delta_, yy - (coeff_ - 1)*delta_))
      coords.3 <- cbind(x = c(xx - (coeff_ - 1)*delta_, xx - (coeff_ - 1)*delta_, xx + coeff_*delta_, xx + coeff_*delta_, xx - (coeff_ - 1)*delta_),
                        y = c(yy - coeff_*delta_, yy - (coeff_ - 1)*delta_, yy - (coeff_ - 1)*delta_, yy - coeff_*delta_, yy - coeff_*delta_))
      coords.4 <- cbind(x = c(xx - coeff_*delta_, xx - coeff_*delta_, xx - (coeff_ - 1)*delta_, xx - (coeff_ - 1)*delta_, xx - coeff_*delta_),
                        y = c(yy - coeff_*delta_, yy + (coeff_ - 1)*delta_, yy + (coeff_ - 1)*delta_, yy - coeff_*delta_, yy - coeff_*delta_))
      coords.list <- append(coords.list, list(coords.1, coords.2, coords.3, coords.4))
      
      poly.1 <- sp::Polygon(coords.1)
      poly.2 <- sp::Polygon(coords.2)
      poly.3 <- sp::Polygon(coords.3)
      poly.4 <- sp::Polygon(coords.4)
      poly.list <- append(poly.list, list(poly.1, poly.2, poly.3, poly.4))
      poly.name <- c(poly.name, paste0(coeff_, c('-A', '-B', '-C', '-D')))
      coeff_ = coeff_ + 1
    }
  }
  
  coords.I <- cbind(x = c(x.min, x.min, xx + n.layer*delta_, xx + n.layer*delta_, x.min),
                    y = c(yy + n.layer*delta_, y.max, y.max, yy + n.layer*delta_, yy + n.layer*delta_))
  coords.J <- cbind(x = c(xx + n.layer*delta_, xx + n.layer*delta_, x.max, x.max, xx + n.layer*delta_),
                    y = c(yy - n.layer*delta_, y.max, y.max, yy - n.layer*delta_, yy - n.layer*delta_))
  coords.K <- cbind(x = c(xx - n.layer*delta_, xx - n.layer*delta_, x.max, x.max, xx - n.layer*delta_),
                    y = c(y.min, yy - n.layer*delta_, yy - n.layer*delta_, y.min, y.min))
  coords.L <- cbind(x = c(x.min, x.min, xx - n.layer*delta_, xx - n.layer*delta_, x.min),
                    y = c(y.min, yy + n.layer*delta_, yy + n.layer*delta_, y.min, y.min))
  
  coords.list <- append(coords.list, list(coords.I, coords.J, coords.K, coords.L))
  
  poly.I <- sp::Polygon(coords.I)
  poly.J <- sp::Polygon(coords.J)
  poly.K <- sp::Polygon(coords.K)
  poly.L <- sp::Polygon(coords.L)
  poly.list <- append(poly.list, list(poly.I, poly.J, poly.K, poly.L))
  poly.name <- c(poly.name, paste0('last', c('-A', '-B', '-C', '-D')))
  
  poly.list <- lapply(seq_along(poly.list), 
                      function(i) Polygons(list(poly.list[[i]]), ID = poly.name[i]))
  
  int_poly <- SpatialPolygons(poly.list)
  
  
  #int_poly <- sp::Polygons(poly.list, ID = 'A')
  #int_poly <- sp::SpatialPolygons(list(int_poly))
  proj4string(int_poly) <- crs_Ita_km
  if(displaypl){
    ggplot()  + gg(int_poly) + gg(ama_bdy) + geom_point(aes(xx,yy)) +  
      coord_fixed()  
  }
  else{
    int_poly
  }
}


## Function to create the spatial integration grid
#' input : xx , yy  --> coordinates of the points
#'         delta_   --> minimum bin edge length (for internal layers)
#'         n.layer  --> maximum number of layer (excluding the last one, each layer contains 4 bins) 
#'         bdy_     --> rectangular region of interest
#'         min.edge --> minimum bin edge length (only for last layer)

#' output : data.frame with columns 
#'          x , y polygon coordinates (each polygon is defined by 5 points)
#'          xx, yy point coordinates (around which the bins are created)
#'          layer.id bin identifier (1-A, 1-B, 1-C,1-D,.....,last-A, last-B, last-C, last-D)
Find_sp_grid <- function(xx, yy, delta_, n.layer, bdy_, min.edge = 0.01){
  # extract boundaries
  x.min <- bdy_@bbox[1,1]
  x.max <- bdy_@bbox[1,2]
  y.min <- bdy_@bbox[2,1]
  y.max <- bdy_@bbox[2,2]
  # if TRUE considers the rectangular defined by the bdy bbox
  if(n.layer == -1){
    coord.bdy <- bdy_@polygons[[1]]@Polygons[[1]]@coords
    coord.df <- data.frame(x = coord.bdy[,1],
                           y = coord.bdy[,2],
                           xx = xx,
                           yy = yy,
                           layer.id = 'last')
    return(coord.df)
  }
  # find max number of layers
  n.layer <- floor(min(c(n.layer, 
                         abs(x.min - xx)/delta_,
                         abs(x.max - xx)/delta_, 
                         abs(y.min - yy)/delta_,
                         abs(y.max - yy)/delta_)))

  #cat('Number of layers : ', n.layer, '\n')
  red_coef = 1
  # reduce the number of layers if last bin has an edge smaller than min.edge
  while(n.layer > 0 && 
        (abs(x.min - (xx - n.layer*delta_)) <= min.edge |
         abs(x.max - (n.layer*delta_ + xx)) <= min.edge |
         abs(y.min - (yy - n.layer*delta_)) <= min.edge |
         abs(y.max - (n.layer*delta_ + yy))  <= min.edge ) ){
    n.layer = n.layer - red_coef
    #cat('Reducing number of layers by ', red_coef, '\n')
    red_coef <- red_coef + 1
  }
  if(n.layer > 0){
    coeff_ = 1
    
    for(layer_ in 1:n.layer){
      # for each layer it creates 4 data.frames with the polygons coordinates
      coords.1 <- data.frame(x = c(xx - coeff_*delta_, xx - coeff_*delta_, xx + (coeff_ - 1)*delta_, xx + (coeff_ - 1)*delta_, xx - coeff_*delta_),
                             y = c(yy + (coeff_ - 1)*delta_, yy + coeff_*delta_, yy + coeff_*delta_, yy + (coeff_ - 1)*delta_, yy + (coeff_ - 1)*delta_),
                             xx = xx,
                             yy = yy,
                             layer.id = paste0(coeff_,'-A'))
      coords.2 <- data.frame(x = c(xx + (coeff_ - 1)*delta_, xx + (coeff_ - 1)*delta_, xx + coeff_*delta_, xx + coeff_*delta_, xx + (coeff_ - 1)*delta_),
                             y = c(yy - (coeff_ - 1)*delta_, yy + coeff_*delta_, yy + coeff_*delta_, yy - (coeff_ - 1)*delta_, yy - (coeff_ - 1)*delta_),
                             xx = xx,
                             yy = yy,
                             layer.id = paste0(coeff_,'-B'))
      coords.3 <- data.frame(x = c(xx - (coeff_ - 1)*delta_, xx - (coeff_ - 1)*delta_, xx + coeff_*delta_, xx + coeff_*delta_, xx - (coeff_ - 1)*delta_),
                             y = c(yy - coeff_*delta_, yy - (coeff_ - 1)*delta_, yy - (coeff_ - 1)*delta_, yy - coeff_*delta_, yy - coeff_*delta_),
                             xx = xx,
                             yy = yy,
                             layer.id = paste0(coeff_,'-C'))
      coords.4 <- data.frame(x = c(xx - coeff_*delta_, xx - coeff_*delta_, xx - (coeff_ - 1)*delta_, xx - (coeff_ - 1)*delta_, xx - coeff_*delta_),
                             y = c(yy - coeff_*delta_, yy + (coeff_ - 1)*delta_, yy + (coeff_ - 1)*delta_, yy - coeff_*delta_, yy - coeff_*delta_),
                             xx = xx,
                             yy = yy,
                             layer.id = paste0(coeff_,'-D'))
      if(coeff_ == 1){
        coord.df <- rbind(coords.1, coords.2, coords.3, coords.4)
      } else{
        coord.df <- rbind(coord.df, coords.1, coords.2, coords.3, coords.4)
      }
      coeff_ = coeff_ + 1
    }
  }
  # last layer is irregular and created separetely
  coords.I <- data.frame(x = c(x.min, x.min, xx + n.layer*delta_, xx + n.layer*delta_, x.min),
                         y = c(yy + n.layer*delta_, y.max, y.max, yy + n.layer*delta_, yy + n.layer*delta_),
                         xx = xx,
                         yy = yy,
                         layer.id = 'last-A')
  coords.J <- data.frame(x = c(xx + n.layer*delta_, xx + n.layer*delta_, x.max, x.max, xx + n.layer*delta_),
                         y = c(yy - n.layer*delta_, y.max, y.max, yy - n.layer*delta_, yy - n.layer*delta_),
                         xx = xx,
                         yy = yy,
                         layer.id = 'last-B')
  coords.K <- data.frame(x = c(xx - n.layer*delta_, xx - n.layer*delta_, x.max, x.max, xx - n.layer*delta_),
                         y = c(y.min, yy - n.layer*delta_, yy - n.layer*delta_, y.min, y.min),
                         xx = xx,
                         yy = yy,
                         layer.id = 'last-C')
  coords.L <- data.frame(x = c(x.min, x.min, xx - n.layer*delta_, xx - n.layer*delta_, x.min),
                         y = c(y.min, yy + n.layer*delta_, yy + n.layer*delta_, y.min, y.min),
                         xx = xx,
                         yy = yy,
                         layer.id = 'last-D')
  # if n.layer = 0 then it creates 4 polygons 
  if(n.layer == 0){
    coord.df <- rbind(coords.I, coords.J, coords.K, coords.L)
  }
  else{
    coord.df <- rbind(coord.df, coords.I, coords.J, coords.K, coords.L)
  }
  coord.df
}

# Creates spatio-temporal integration grid given a data.point 
#' input  : data.point --> row of a data.frame representing an observations (x,y,ts,mags,idx.p)
#'          delta. --> minimum spatial bin edge length (for internal layers)
#'          n.layer. --> maximum number of spatial bins layers (excluding the last one, each layer contains 4 bins) 
#'          min.edge. --> minimum spatial bin edge length (for last layer)
#'          coef.t --> coefficient for the time bins
#'          max_len.t --> length of the first time bin, multiplied by (1 + coef.t)^k to find other bins
#'          T2. --> boundary of the time region
#'          bdy. --> polygon representing the spatial region
#' output : data.frame for the poisson surrogate model representing the number of triggered events in the 
#'          region (Part II of the log-lik) for a single point. It has the same columns as data.point plus 
#'          x.min, x.max, y.min, y.max representing the spatial region to be integrated
#'          t.start, t.end representing the time region to be integrated
#'          bin.name unique name for each bin (depends on idx.p)
#'          ref_layer name for calculations, each bin with the same value of ref_layer will have the same spatial integral
#'          t.ref_layer name for calculations, each bin with the same value of t.ref_layer will have the same time integral

space.time.grid <- function(data.point, delta., n.layer., min.edge., coef.t, max_len.t,
                            T2., bdy., displaygrid = FALSE){
  xx. <- data.point$x
  yy. <- data.point$y
  tt. <- data.point$ts
  idx.p <- data.point$idx.p
  # spatial bins
  if(displaygrid){
    Plot_grid(xx = xx., yy = yy., delta_ = delta., n.layer = n.layer., 
              bdy_ =  bdy., min.edge = min.edge.)
  }
  # find coordinates of polygons : each bin is represented by 5 coordinates
  sp.bins.coords <- Find_sp_grid(xx = xx., yy = yy., delta_ = delta., n.layer = n.layer., 
                                 bdy_ =  bdy., min.edge = min.edge.)
  # store bins names
  layer.name <- unique(sp.bins.coords$layer.id)
  # for each bin we keep only x.min, x.max, y.min, y.max, the reference x,y and the bin name 
  space.bins <- foreach(l.name = layer.name, .combine = rbind) %do% {
    df.sel <- sp.bins.coords[sp.bins.coords$layer.id == l.name,]
    data.frame(x.min = min(df.sel$x),
               y.min = min(df.sel$y),
               x.max = max(df.sel$x),
               y.max = max(df.sel$y),
               bin.name = paste0(l.name,'_',idx.p))
  }
  
  # time bins
  # find bins break points
  t_b <- breaks_exp(tt., T2., coef_ = coef.t, max_length_ = max_len.t)
  time.bins <- data.frame(t.start = t_b[-length(t_b)], 
                          t.end = t_b[-1]) %>%
    mutate(t.bin.name = paste0(round(t.start,3),'-',round(t.end,3)))
  time.bins$t.ref_layer = c(1:(nrow(time.bins) - 1), paste0('last-',idx.p))
  spt.grid <- suppressWarnings(cbind(merge(time.bins, space.bins), data.point))
  spt.grid$ref_layer <- spt.grid$bin.name
  spt.grid$ref_layer[!grepl('last', spt.grid$ref_layer)] <- substr(spt.grid$ref_layer[!grepl('last', spt.grid$ref_layer)], 1, 1)
  spt.grid
}




# function to calculate the integral of 2-D Gaussian distribution for each row of the data.frame
#' input  : coord.df --> spatial.time.grid output
#'          Sigma    --> covariance matrix of the Gaussian distribution (matrix or list of matrices)
#' output : vector of length nrow(coord.df) each element is the integral of a Gaussian distribution with 
#'          given Sigma, and varying mean and integration regions.

I.s_df <- function(coord.df, Sigma = diag(c(1,1))){
  # if Sigma is just a matrix creates a list with repeated Sigma
  if(is.matrix(Sigma)){
    Sigma.l <- lapply(1:nrow(coord.df), \(x) Sigma)
  } 
  # if it is a list check if its the right length
  else if(is.list(Sigma)){ 
    Sigma.l <- Sigma
    if(length(Sigma) != nrow(coord.df)){
      stop('Incorrect list length, length(Sigma) != nrow(coord.df) - they have to be equal')
    }
  }
  # if it is not a list or a matrix return error.
  else {
    stop('Unknown Sigma format - please provide it as a matrix or a list of matrices')
  }
  # extract coordinates (the mean of the Gaussian distribution to be integrated)
  p.loc <- cbind(coord.df$x, coord.df$y)
  # extract boundaries of the polygons to be integrated over
  x.min <- coord.df$x.min
  x.max <- coord.df$x.max
  y.min <- coord.df$y.min
  y.max <- coord.df$y.max
  #print(Sigma.l)
  # calculate the integral for each row of the data.frame
  sapply(1:nrow(p.loc), \(idx)
         abs(pmvnorm(lower = c(x.min[idx], y.min[idx]), 
                     upper = c(x.max[idx], y.max[idx]), 
                     mean = p.loc[idx,], sigma = Sigma.l[[idx]],
                     keepAttr = FALSE)))
}


It_df <- function(param_, time.df){
  tth <- as.numeric(time.df$ts)
  T1b <- as.numeric(time.df$t.start)
  T2b <- as.numeric(time.df$t.end)
  param_c <- param_[4]
  param_p <- param_[5]
  T.l <- sapply(1:length(tth), \(x) max(tth[x], T1b[x]))
  fun.l <- (1 + (T.l - tth)/param_c)^(1-param_p)
  fun.u <- (1 + (T2b - tth)/param_c)^(1-param_p)
  ( param_c/ (param_p - 1) )* ( fun.l - fun.u )
}

compute.grid <- function(param., gridd, Sigma_ = diag(c(1,1))){
  
  t.names <- unique(gridd$t.ref_layer)
  sp.names <- unique(gridd$ref_layer)
  time.sel <- gridd[sapply(t.names, \(bname) which(gridd$t.ref_layer == bname)[1]),]
  space.sel <- gridd[sapply(sp.names, \(bname) which(gridd$ref_layer == bname)[1]),]
  if(is.list(Sigma_)){
    Sigma.sp <- Sigma_[sapply(sp.names, \(bname) which(gridd$ref_layer == bname)[1])]
    Is.vec <- I.s_df(space.sel, Sigma.sp)
  } else {
    Is.vec <- I.s_df(space.sel, Sigma_)
  }
  
  
  It.vec <- It_df(param_ = param., time.df = time.sel)
  
  list.ret <- list(It = It.vec[match(gridd$t.ref_layer, t.names)],
                   Is = Is.vec[match(gridd$ref_layer, sp.names)],
                   space.df = space.sel %>% 
                     mutate(Is = Is.vec),
                   time.df = time.sel %>%
                     mutate(It = It.vec))
  return(list.ret)
}




### functions for model fitting
ETAS.fit.B_bkg_smart <- function(sample.s, 
                                 M0, T1, T2, bdy, 
                                 Sigma = NULL,
                                 coef.t = 1.5,
                                 delta_exp = 0.1,
                                 max_len.t = 0.1,
                                 n.layer.sp = -1, min.edge.sp = 1, 
                                 N_exp = 50,
                                 link.fun = list(mu = \(x) x,
                                                 K = \(x) x,
                                                 alpha = \(x) x,
                                                 cc = \(x) x,
                                                 pp = \(x) x), 
                           bru.opt = list(bru_verbose = 3,
                                          bru_max_iter = 50),
                           bin.strat = 'def',
                           ncore = 1){
  
  # create different poisson counts models
  # first for background
  df.0 <- data.frame(counts = 0, exposures = 1)
  # this is the expression of log(Lambda0)
  form.0 <- counts ~ link.fun$mu(th.mu) + log(T2 - T1) #+ log(Area.bdy)
  
  # first likelihood
  lik.0 <- inlabru::like(formula = form.0,
                         data = df.0,
                         family = 'poisson',
                         options = list(E = df.0$exposures))
  
  
  ### create time bins 
  cat('Start creating spatial grid...', '\n')
  time.g.st <- Sys.time()
  df.j <- foreach(idx = 1:nrow(sample.s), .combine = rbind) %do% {
    space.time.grid(data.point = sample.s[idx,], 
                    delta. = delta_exp, 
                    n.layer. = n.layer.sp, 
                    min.edge. = 1,
                    coef.t = 1.5,
                    max_len.t = 0.1,
                    T2. = T2,
                    bdy. = bdy, displaygrid = FALSE)
  } %>%
    mutate(counts = 0,
           exposures = 1)
  cat('Finished creating spatial grid, time ', Sys.time() - time.g.st, '\n')   
  
  logLambda.h.inla <- function(th.K, th.alpha, th.c, th.p, df_grid, 
                               M0, bdy, Sigma., ncore_ = ncore){
    theta_ <- c(0, 
                link.fun$K(th.K[1]), 
                link.fun$alpha(th.alpha[1]), 
                link.fun$cc(th.c[1]), 
                link.fun$pp(th.p[1]))
    #cat('theta - LogL', theta_, '\n')
    comp.list <- compute.grid(param. = theta_, gridd = df_grid, Sigma_ = Sigma.)
    
    theta_[3]*(df_grid$mags - M0) + log1p(theta_[2] - 1) + log1p(comp.list$It - 1) + 
                 log1p(comp.list$Is - 1)
    
  }
  
  #debugonce(logLambda.h.inla)
  # creating formula for past events contributions to integrated lambda
  form.j.part <- counts ~ logLambda.h.inla(th.K = th.K, 
                                           th.alpha = th.alpha, 
                                           th.c = th.c, 
                                           th.p = th.p, 
                                           df_grid = df.j,
                                           M0 = M0, bdy = bdy, 
                                           Sigma. = Sigma)
  # second for triggered part of the integral
  lik.j.part <- inlabru::like(formula = form.j.part,
                              data = df.j,
                              family = 'poisson',
                              options = list(E = df.j$exposures)) 
  
  # third is for the sum of the log intensities
  df.s <- data.frame(x = sample.s$x, y = sample.s$y, ts = sample.s$ts, mags = sample.s$mags, 
                     bkg = sample.s$bkg, counts = 1, exposures = 0)
  
  loglambda.inla <- function(th.mu, th.K, th.alpha, th.c, th.p, bkg, tt, xx, yy, 
                             th, xh, yh, mh, M0, Sigma, ncore_ = ncore){
    
    
    th.p.df <- data.frame(th.mu = link.fun$mu(th.mu[1])*bkg, 
                          th.K = link.fun$K(th.K[1]), 
                          th.alpha = link.fun$alpha(th.alpha[1]), 
                          th.c = link.fun$cc(th.c[1]), 
                          th.p = link.fun$pp(th.p[1]))
    #cat('th.mu', sum(th.p.df$th.mu <0), '\n')
    
    outp <- unlist(mclapply(1:length(tt), \(idx) 
                            log(lambda(theta = as.numeric(th.p.df[idx,]), 
                                       tt = tt[idx], 
                                       xx = xx[idx], 
                                       yy = yy[idx], 
                                       th = th, xh = xh, yh = yh, mh = mh, M0 = M0, Sigma = Sigma)),
                            mc.cores = ncore_))
    #cat('logl', sum(is.na(outp)), '\n')
    outp
  }
  # creating formula for summation part
  form.s.part <- counts ~ loglambda.inla(th.mu = th.mu, 
                                         th.K = th.K, 
                                         th.alpha = th.alpha, 
                                         th.c = th.c, 
                                         th.p = th.p,
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
    th.mu(1, model = 'linear', mean.linear = 0, prec.linear = 1) + 
    th.K(1, model = 'linear', mean.linear = 0, prec.linear = 1) +
    th.alpha(1, model = 'linear', mean.linear = 0, prec.linear = 1) +
    th.c(1, model = 'linear', mean.linear = 0, prec.linear = 1) +
    th.p(1, model = 'linear', mean.linear = 0, prec.linear = 1) 
  
  bru(lik.0, lik.j.part, lik.s.part, components = cmp.part, 
      options = bru.opt)
  
  
}


### functions for model fitting
ETAS.fit.isotropic <- function(sample.s, 
                               M0, T1, T2, bdy,
                               coef_t = 1, delta.t = 0.1, N_exp = 50,
                               delta.sp = 0.1, n.layer.sp = -1, min.edge.sp = 1, 
                              link.fun = list(mu = \(x) x,
                                                 K = \(x) x,
                                                 alpha = \(x) x,
                                                 cc = \(x) x,
                                                 pp = \(x) x,
                                                 sigma = \(x) x), 
                                 bru.opt = list(bru_verbose = 3,
                                                bru_max_iter = 50),
                                 bin.strat = 'def',
                                 ncore = 1){
  
  # create different poisson counts models
  # first for background
  df.0 <- data.frame(counts = 0, exposures = 1)
  # this is the expression of log(Lambda0)
  form.0 <- counts ~ link.fun$mu(th.mu) + log(T2 - T1) #+ log(Area.bdy)
  
  # first likelihood
  lik.0 <- inlabru::like(formula = form.0,
                         data = df.0,
                         family = 'poisson',
                         options = list(E = df.0$exposures))
  
  
  ### create time bins 
  cat('Start creating spatial grid...', '\n')
  time.g.st <- Sys.time()
  df.j <- foreach(idx = 1:nrow(sample.s), .combine = rbind) %do% {
    space.time.grid(data.point = sample.s[idx,], 
                    delta. = delta.sp, 
                    n.layer. = n.layer.sp, 
                    min.edge. = 0.1,
                    coef.t = coef_t,
                    max_len.t = delta.t,
                    T2. = T2,
                    bdy. = bdy, displaygrid = FALSE)
  } %>%
    mutate(counts = 0,
           exposures = 1)
  cat('Finished creating spatial grid, time ', Sys.time() - time.g.st, '\n')   
  
  logLambda.h.inla <- function(th.K, th.alpha, th.c, th.p, th.sigma, df_grid, 
                               M0, bdy, ncore_ = ncore){
    theta_ <- c(0, 
                link.fun$K(th.K[1]), 
                link.fun$alpha(th.alpha[1]), 
                link.fun$cc(th.c[1]), 
                link.fun$pp(th.p[1]))

    sigma.1 <- link.fun$sigma(th.sigma[1])
    Sigma. <- matrix(c(sigma.1, 0, 0, sigma.1), byrow = TRUE, ncol = 2)
    
    
    #cat('theta - LogL', theta_, '\n')
    comp.list <- compute.grid(param. = theta_, gridd = df_grid, Sigma_ = Sigma.)
    
    # cat('Is.na = ', sum(is.na(log1p(comp.list$Is - 1))), "\n")
    # cat('It.na = ', sum(is.na(log1p(comp.list$It - 1))), '\n')
    # cat('Is.inf = ', sum(is.infinite(log1p(comp.list$Is - 1))), "\n")
    # cat('It.inf = ', sum(is.infinite(log1p(comp.list$It - 1))), '\n')
    # 
    theta_[3]*(df_grid$mags - M0) + log(theta_[2]) + log(comp.list$It) + 
      log(comp.list$Is)
    
  }
  
  #debugonce(logLambda.h.inla)
  # creating formula for past events contributions to integrated lambda
  form.j.part <- counts ~ logLambda.h.inla(th.K = th.K, 
                                           th.alpha = th.alpha, 
                                           th.c = th.c, 
                                           th.p = th.p, 
                                           th.sigma = th.sigma,
                                           df_grid = df.j,
                                           M0 = M0, bdy = bdy)
  # second for triggered part of the integral
  lik.j.part <- inlabru::like(formula = form.j.part,
                              data = df.j,
                              family = 'poisson',
                              options = list(E = df.j$exposures)) 
  
  # third is for the sum of the log intensities
  df.s <- data.frame(x = sample.s$x, y = sample.s$y, ts = sample.s$ts, mags = sample.s$mags, 
                     bkg = sample.s$bkg, counts = 1, exposures = 0)
  
  loglambda.inla <- function(th.mu, th.K, th.alpha, th.c, th.p, th.sigma, bkg, tt, xx, yy, 
                             th, xh, yh, mh, M0, ncore_ = ncore){
    
    
    th.p.df <- data.frame(th.mu = link.fun$mu(th.mu[1])*bkg, 
                          th.K = link.fun$K(th.K[1]), 
                          th.alpha = link.fun$alpha(th.alpha[1]), 
                          th.c = link.fun$cc(th.c[1]), 
                          th.p = link.fun$pp(th.p[1]))
    #cat('th.mu', sum(th.p.df$th.mu <0), '\n')
    sigma.1 <- link.fun$sigma(th.sigma[1])
    
    Sigma. <- matrix(c(sigma.1, 0, 0, sigma.1), byrow = TRUE, ncol = 2)
    outp <- unlist(mclapply(1:length(tt), \(idx) 
                            log(lambda(theta = as.numeric(th.p.df[idx,]), 
                                       tt = tt[idx], 
                                       xx = xx[idx], 
                                       yy = yy[idx], 
                                       th = th, xh = xh, yh = yh, mh = mh, M0 = M0, Sigma = Sigma.)),
                            mc.cores = ncore_))
    #cat('logl', sum(is.na(outp)), '\n')
    outp
  }
  # creating formula for summation part
  form.s.part <- counts ~ loglambda.inla(th.mu = th.mu, 
                                         th.K = th.K, 
                                         th.alpha = th.alpha, 
                                         th.c = th.c, 
                                         th.p = th.p,
                                         th.sigma = th.sigma,
                                         bkg = bkg,
                                         tt = sample.s$ts, xx = sample.s$x, yy = sample.s$y,
                                         th = sample.s$ts, xh = sample.s$x, yh = sample.s$y, 
                                         mh = sample.s$mags, M0 = M0)
  lik.s.part <- inlabru::like(formula = form.s.part,
                              data = df.s,
                              family = 'poisson',
                              options = list(E = df.s$exposures)) 
  cmp.part <- counts ~ -1 + 
    th.mu(1, model = 'linear', mean.linear = 0, prec.linear = 1) + 
    th.K(1, model = 'linear', mean.linear = 0, prec.linear = 1) +
    th.alpha(1, model = 'linear', mean.linear = 0, prec.linear = 1) +
    th.c(1, model = 'linear', mean.linear = 0, prec.linear = 1) +
    th.p(1, model = 'linear', mean.linear = 0, prec.linear = 1) +
    th.sigma(1, model = 'linear', mean.linear = 0, prec.linear = 1)
  
  bru(lik.0, lik.j.part, lik.s.part, components = cmp.part, 
      options = bru.opt)
  
  
}












### functions for model fitting
ETAS.fit.iso_mags <- function(sample.s, 
                               M0, T1, T2, bdy,
                               coef_t = 1, delta.t = 0.1, N_exp = 50,
                               delta.sp = 0.1, n.layer.sp = -1, min.edge.sp = 1, 
                               link.fun = list(mu = \(x) x,
                                               K = \(x) x,
                                               alpha = \(x) x,
                                               cc = \(x) x,
                                               pp = \(x) x,
                                               sigma = \(x) x), 
                               bru.opt = list(bru_verbose = 3,
                                              bru_max_iter = 50),
                               bin.strat = 'def',
                               ncore = 1){
  
  # create different poisson counts models
  # first for background
  df.0 <- data.frame(counts = 0, exposures = 1)
  # this is the expression of log(Lambda0)
  form.0 <- counts ~ link.fun$mu(th.mu) + log(T2 - T1) #+ log(Area.bdy)
  
  # first likelihood
  lik.0 <- inlabru::like(formula = form.0,
                         data = df.0,
                         family = 'poisson',
                         options = list(E = df.0$exposures))
  
  
  ### create time bins 
  cat('Start creating spatial grid...', '\n')
  time.g.st <- Sys.time()
  df.j <- foreach(idx = 1:nrow(sample.s), .combine = rbind) %do% {
    space.time.grid(data.point = sample.s[idx,], 
                    delta. = delta.sp, 
                    n.layer. = n.layer.sp, 
                    min.edge. = 0.1,
                    coef.t = coef_t,
                    max_len.t = delta.t,
                    T2. = T2,
                    bdy. = bdy, displaygrid = FALSE)
  } %>%
    mutate(counts = 0,
           exposures = 1,
           mags.round = round(mags, 1)) 
  
  df.j$ref_layer = paste0(df.j$ref_layer, '_mag',df.j$mags.round)
  
  cat('Finished creating spatial grid, time ', Sys.time() - time.g.st, '\n')   
  
  logLambda.h.inla <- function(th.K, th.alpha, th.c, th.p, th.sigma, df_grid, 
                               M0., bdy, ncore_ = ncore){
    theta_ <- c(0, 
                link.fun$K(th.K[1]), 
                link.fun$alpha(th.alpha[1]), 
                link.fun$cc(th.c[1]), 
                link.fun$pp(th.p[1]))
    
    Sigma. <- lapply(1:nrow(df_grid), \(idx.g) 
                     link.fun$sigma(th.sigma[1], df_grid$mags.round[idx.g], M0.) )
    #print(Sigma.[[1]])
    
    #cat('theta - LogL', theta_, '\n')
    comp.list <- compute.grid(param. = theta_, gridd = df_grid, Sigma_ = Sigma.)
  
#    print(comp.list$time.df)
 #   print(comp.list$space.df)
    #  cat('Is.na = ', sum(is.na(log1p(comp.list$Is - 1))), "\n")
    #  cat('It.na = ', sum(is.na(log1p(comp.list$It - 1))), '\n')
    #  cat('Is.inf = ', sum(is.infinite(log1p(comp.list$Is - 1))), "\n")
    #  cat('It.inf = ', sum(is.infinite(log1p(comp.list$It - 1))), '\n')
    # # 
    theta_[3]*(df_grid$mags - M0.) + log(theta_[2] + 1e-100) + 
      log(comp.list$It + 1e-100) + 
      log(comp.list$Is + 1e-100)
    }
  
  #debugonce(logLambda.h.inla)
  # creating formula for past events contributions to integrated lambda
  form.j.part <- counts ~ logLambda.h.inla(th.K = th.K, 
                                           th.alpha = th.alpha, 
                                           th.c = th.c, 
                                           th.p = th.p, 
                                           th.sigma = th.sigma,
                                           df_grid = df.j,
                                           M0. = M0, bdy = bdy)
  # second for triggered part of the integral
  lik.j.part <- inlabru::like(formula = form.j.part,
                              data = df.j,
                              family = 'poisson',
                              options = list(E = df.j$exposures)) 
  
  # third is for the sum of the log intensities
  df.s <- data.frame(x = sample.s$x, y = sample.s$y, ts = sample.s$ts, mags = sample.s$mags, 
                     bkg = sample.s$bkg, counts = 1, exposures = 0)
  
  loglambda.inla <- function(th.mu, th.K, th.alpha, th.c, th.p, th.sigma, bkg, tt, xx, yy,
                             mw,
                             th, xh, yh, mh, M0., ncore_ = ncore){
    
    
    th.p.df <- data.frame(th.mu = link.fun$mu(th.mu[1])*bkg, 
                          th.K = link.fun$K(th.K[1]), 
                          th.alpha = link.fun$alpha(th.alpha[1]), 
                          th.c = link.fun$cc(th.c[1]), 
                          th.p = link.fun$pp(th.p[1]))
    #cat('th.mu', sum(th.p.df$th.mu <0), '\n')
    
    mags.round <- round(mw, 1)
    Sigma. <- lapply(1:length(mags.round), \(idx.g) 
                     link.fun$sigma(th.sigma[1], mags.round[idx.g], M0.) )
    
    outp <- unlist(mclapply(1:length(tt), \(idx) 
                            log(lambda(theta = as.numeric(th.p.df[idx,]), 
                                       tt = tt[idx], 
                                       xx = xx[idx], 
                                       yy = yy[idx], 
                                       th = th, xh = xh, yh = yh, mh = mh, M0 = M0., 
                                       Sigma = Sigma.[[idx]])),
                            mc.cores = ncore_))
    #cat('logl', sum(is.na(outp)), '\n')
    outp
  }
  # creating formula for summation part
  form.s.part <- counts ~ loglambda.inla(th.mu = th.mu, 
                                         th.K = th.K, 
                                         th.alpha = th.alpha, 
                                         th.c = th.c, 
                                         th.p = th.p,
                                         th.sigma = th.sigma,
                                         bkg = bkg,
                                         tt = sample.s$ts, xx = sample.s$x, 
                                         yy = sample.s$y,
                                         mw = sample.s$mags,
                                         th = sample.s$ts, xh = sample.s$x, yh = sample.s$y, 
                                         mh = sample.s$mags, M0. = M0)
  lik.s.part <- inlabru::like(formula = form.s.part,
                              data = df.s,
                              family = 'poisson',
                              options = list(E = df.s$exposures)) 
  cmp.part <- counts ~ -1 + 
    th.mu(1, model = 'linear', mean.linear = 0, prec.linear = 1) + 
    th.K(1, model = 'linear', mean.linear = 0, prec.linear = 1) +
    th.alpha(1, model = 'linear', mean.linear = 0, prec.linear = 1) +
    th.c(1, model = 'linear', mean.linear = 0, prec.linear = 1) +
    th.p(1, model = 'linear', mean.linear = 0, prec.linear = 1) +
    th.sigma(1, model = 'linear', mean.linear = 0, prec.linear = 1)
  
  bru(lik.0, lik.j.part, lik.s.part, components = cmp.part, 
      options = bru.opt)
  
  
}


## integrated time-triggering function

It <- function(theta, th, T2){
  gamma.u <- (T2 - th)/theta[4]
  ( theta[4]/(theta[5] - 1) )*(1 - (gamma.u + 1)^(1-theta[5]) )
}


## inverse of integrated time-triggering function
Inv.It <- function(theta, omega, th){
  th + theta[4]*( ( 1 - omega * (1/theta[4])*(theta[5] - 1) )^( -1/(theta[5] - 1) ) - 1)
}


## sampling times
sample.time <- function(theta, n.ev, th, T2){
  if(n.ev == 0){
    df <- data.frame(ts = 1, x = 1, y = 1, mags = 1, gen = 0)
    return(df[-1,])
  }
  bound.l <- 0 #It(th.p, th, T)
  bound.u <- It(theta, th, T2)
  unif.s <- runif(n.ev, min = bound.l, max = bound.u)
  t.sample <- Inv.It(theta, unif.s, th)
  t.sample
}

# sampling locations
sample.loc <- function(xh, yh, n.ev, bdy, Sigma, crs_obj. = crs_Ita_km){
  
  # initialize empty SpatialPoints
  samp.points <- SpatialPoints(data.frame(x = 0, y = 0))
  proj4string(samp.points) <- crs_obj.
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
sample.triggered <- function(theta, beta.p, th, xh, yh, n.ev, M0, T1, T2, bdy, 
                             Sigma){
  # if the number of events to be placed is zero returns an empty data.frame
  if(n.ev == 0){
    samp.points <- data.frame(x = 1, y = 1, ts = 1, mags = 1)
    samp.points <- samp.points[-1,]
    return(samp.points)
  }
  else{
    
    # sample times
    samp.ts <- sample.time(theta, n.ev, th, T2)
    # sample magnitudes
    samp.mags <- rexp(n.ev, rate = beta.p) + M0
    # sample locations
    samp.locs <- sample.loc(xh = xh, yh = yh, n.ev = n.ev, 
                            bdy = bdy, Sigma = Sigma)
    
    # build output dataset
    samp.points <- data.frame(ts = samp.ts, mags = samp.mags, x = samp.locs@coords[,1],
                              y = samp.locs@coords[,2])
    # return only the ones with time different from NA (the one with NA are outside the interval T1, T2)
    # even though it should not happen given how we built sample.omori
    return(samp.points[!is.na(samp.points$ts),])
  }
}


sample.generation <- function(theta, beta.p, Ht, M0, T1, T2, bdy,
                              Sigma, ncore = 1){
  
  # number of parents
  n.parent <- nrow(Ht)
  # calculate the aftershock rate for each parent in history
  if(is.matrix(Sigma)){
    trig.rates <- exp(logLambda.h.vec(theta, Ht$ts, Ht$x, Ht$y, Ht$mags, M0, T1, T2, bdy,
                                      Sigma))
    Sigma_l <- lapply(1:nrow(Ht), \(x) Sigma)
  }
  if(is.null(Sigma)){
    Sigma <- lapply(1:nrow(Ht), \(x) matrix(c(theta[6]*(Ht$mags[x] - M0), 0,
                                              0, theta[6]*(Ht$mags[x] - M0)), byrow = TRUE, 
                                            ncol = 2))
    trig.rates <- unlist(mclapply(1:nrow(Ht), \(idx)
                                  exp(logLambda.h.vec(theta, Ht$ts[idx], 
                                                      Ht$x[idx], Ht$y[idx], 
                                                      Ht$mags[idx], M0, T1, T2, bdy,
                                      Sigma[[idx]])) ))
    Sigma_l <- Sigma
  }
  # extract number of aftershock for each parent
  n.ev.v <- sapply(1:n.parent, function(x) rpois(1, trig.rates[x]))
  if(sum(n.ev.v) > 1000){
    print('Warning : More than 1000 triggered events')
    app <- data.frame(x = 1, y = 1, ts = 1, mags = 1)
    app <- app[-1,]
    return(app)
  }
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
    sample.triggered(theta, beta.p, Ht$ts[idx], Ht$x[idx], Ht$y[idx], n.ev.v[idx], M0, T1, T2, bdy, 
                     Sigma_l[[idx]]), mc.cores = ncore)
  
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
  loglambda_max <- max(loglambda)
  ## number of points to sample at a time - might want to adjust depending on how many points you want to actually retain.
  n.points=10000
  
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
    lambda_ratio <- exp(as.vector(proj$A %*% loglambda) - loglambda_max)
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


sample.ETAS <- function(theta, beta.p, M0, T1, T2, bdy, Sigma = NULL, 
                        loglambda.bkg = NULL, mesh.bkg = NULL, 
                        crs_obj = NULL, Ht = NULL, ncore = 1,
                        Unit.Area = TRUE){
  
  # if the upper extreme greater than lower
  if(T2 < T1){
    stop('Error - right-end of time interval greater than left-end')
  }
  
  if(Unit.Area){
    n.bkg <- rpois(1, theta[1]*(T2 - T1))
  }
  else{
    bdy.sf <- st_as_sf(bdy.sf)
    st_crs(bdy.sf) <- italy.crs
    bdy.sf <- st_transform(bdy.sf, crs_obj)
    Area.bdy <- as.numeric(st_area(bdy.sf)/(1000^2))
    n.bkg <- rpois(1, theta[1]*(T2 - T1)*Area.bdy)
  }
  
  #cat('Background : ', n.bkg, '\n')
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
    gen.from.past <- sample.generation(theta, beta.p, Ht, M0, T1, T2,  bdy, Sigma, ncore)
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
    #cat('Gen : ', nrow(parents), '\n')
    #print(c(T1,T2))
    #print(range(parents$ts))
    # generate aftershocks
    triggered <- sample.generation(theta, beta.p, parents, 
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



##################################
 
# functions to extract information from models


Lambda.pred.mags <- function(th.mu, th.K, th.alpha, th.c, th.p, th.sigma, df.obs, 
                        T1., T2., bdy., M0., link_f){
  df.obs <- df.obs[df.obs$ts < T2.,]
  param. <- c(link_f$mu(th.mu[1]), 
              link_f$K(th.K[1]), 
              link_f$alpha(th.alpha[1]), 
              link_f$cc(th.c[1]), 
              link_f$pp(th.p[1])) 
  
  mags.round <- round(df.obs$mags, 1)
  Sigma. <- lapply(1:length(mags.round), \(idx.g) 
                   link_f$sigma(th.sigma[1], mags.round[idx.g], M0.) )
  
  bkg.part <- param.[1]*(T2 - T1)
  trig.part <- sum(unlist(mclapply(1:length(Sigma.), \(idx.p) 
                                 exp(logLambda.h.vec(theta = param.,
                                       th = df.obs$ts[idx.p],
                                       xh = df.obs$x[idx.p],
                                       yh = df.obs$y[idx.p],
                                       mh = df.obs$mags[idx.p],
                                       M0 = M0.,
                                       T1 = T1.,
                                       T2 = T2.,
                                       bdy = bdy.,
                                       Sigma = Sigma.[[idx.p]])), mc.cores = 5) 
                          ))
  list(full = bkg.part + trig.part,
       background = bkg.part, 
       triggered = trig.part)
  
} 



Lambda.pred <- function(th.mu, th.K, th.alpha, th.c, th.p, th.sigma, df.obs, 
                           T1., T2., bdy., M0., link_f){
  df.obs <- df.obs[df.obs$ts < T2.,]
  param. <- c(link_f$mu(th.mu[1]), 
              link_f$K(th.K[1]), 
              link_f$alpha(th.alpha[1]), 
              link_f$cc(th.c[1]), 
              link_f$pp(th.p[1])) 
  sigma.1 <- link_f$sigma(th.sigma[1])
  Sigma. <- matrix(c(sigma.1, 0, 
                     0, sigma.1), ncol = 2, byrow = TRUE)
  
  bkg.part <- param.[1]*(T2 - T1)
  trig.part <- sum(exp(logLambda.h.vec(theta = param.,
                                       th = df.obs$ts,
                                       xh = df.obs$x,
                                       yh = df.obs$y,
                                       mh = df.obs$mags,
                                       M0 = M0.,
                                       T1 = T1.,
                                       T2 = T2.,
                                       bdy = bdy.,
                                       Sigma = Sigma.))) 
  list(full = bkg.part + trig.part,
       background = bkg.part, 
       triggered = trig.part)
  
} 

Lambda.pred.sk <- function(th.mu, th.K, th.alpha, th.c, th.p, th.sigma, df.obs, 
                             T1., T2., bdy., M0., link_f, Sigma.){
  df.obs <- df.obs[df.obs$ts < T2.,]
  param. <- c(link_f$mu(th.mu[1]), 
              link_f$K(th.K[1]), 
              link_f$alpha(th.alpha[1]), 
              link_f$cc(th.c[1]), 
              link_f$pp(th.p[1])) 
  
  bkg.part <- param.[1]*(T2 - T1)
  trig.part <- sum(exp(logLambda.h.vec(theta = param.,
                                       th = df.obs$ts,
                                       xh = df.obs$x,
                                       yh = df.obs$y,
                                       mh = df.obs$mags,
                                       M0 = M0.,
                                       T1 = T1.,
                                       T2 = T2.,
                                       bdy = bdy.,
                                       Sigma = Sigma.))) 
  list(full = bkg.part + trig.part,
       background = bkg.part, 
       triggered = trig.part)
  
} 

Lambda_fun <- function(param., df.obs, 
                       T1., T2., bdy., M0.,
                       part = 'full'){
  df.obs <- df.obs[df.obs$ts < T2.,]
  
  Sigma_m <- matrix(c(param.[6], 0,
                      0, param.[6]), byrow = TRUE, ncol = 2)
  
  bkg.part <- param.[1]*(T2 - T1)
  trig.part <- sum(exp(logLambda.h.vec(theta = param.,
                                       th = df.obs$ts,
                                       xh = df.obs$x,
                                       yh = df.obs$y,
                                       mh = df.obs$mags,
                                       M0 = M0.,
                                       T1 = T1.,
                                       T2 = T2.,
                                       bdy = bdy.,
                                       Sigma = Sigma_m))) 
  list(full = bkg.part + trig.part,
       background = bkg.part, 
       triggered = trig.part)
  
} 

Lambda_fun.ms <- function(param., df.obs, 
                       T1., T2., bdy., M0.){
  df.obs <- df.obs[df.obs$ts < T2.,]
  
  Sigma.list <- lapply(round(df.obs$mags,1), \(mh) 
                       matrix(c(param.[6]*(mh - M0), 0,
                                0, param.[6]*(mh - M0)), byrow = TRUE, ncol = 2) ) 
  
  bkg.part <- param.[1]*(T2 - T1)
  trig.part <- sum(exp(unlist(mclapply(1:nrow(df.obs), \(idx.p)
    logLambda.h.vec(theta = param.,
                    th = df.obs$ts[idx.p],
                    xh = df.obs$x[idx.p],
                    yh = df.obs$y[idx.p],
                    mh = df.obs$mags[idx.p],
                    M0 = M0.,
                    T1 = T1.,
                    T2 = T2.,
                    bdy = bdy.,
                    Sigma = Sigma.list[[idx.p]]), mc.cores = 5)) ))
  list(full = bkg.part + trig.part,
       background = bkg.part, 
       triggered = trig.part)
  
} 


fore.batch <- function(fit_, 
                       total.data,
                       imposed.data = data.frame(),
                       T1.f, T2.f, 
                       M0, 
                       sample.param = NULL,
                       n_sample = 1000,
                       bdy_ = NULL, 
                       Sigma_ = NULL,
                       log.bkg.field_ = NULL, 
                       mesh_ = mesh_col, 
                       beta.p = beta.p_,
                       crs_ = italy.crs, ncore = 1,
                       link.fun,
                       file.name,
                       time.out.sec = 60,
                       display_ = TRUE
){
  # if posterior sample of the parameter is NULL extract it
  if(is.null(sample.param)){
    
    if(display_){
      cat('No posterior samples provided - sampling posterior ', n_sample , ' times', '\n')
    }
    
    sample.param = generate(fit_, data.frame(x = 0, y = 0, ts = 0, mags = 0), 
                            ~ c(th.mu, th.K, th.alpha, th.c, th.pm1), n.samples = n_sample)
    
    if(display_){
      cat('Finished posterior sampling', '\n')
    }
  }
  # select all the data recorded before the forecasting period and merged it with imposed data
  known.data = rbind(total.data[total.data$ts < T1.f, ], imposed.data)
  #print(known.data)
  if(display_){
    cat('Start simulating ', ncol(sample.param), ' catalogues - timeout ', time.out.sec , 
        ' seconds', '\n', 'Imposed events : ', nrow(imposed.data), ' - Known events : ',
        nrow(known.data), '\n')
  }
  
  # for each posterior sample of the parameters
  for(i in 1:ncol(sample.param)){
    
    # transform parameters in the ETAS scale
    if(is.null(Sigma_)){
      th_ = c(link.fun$mu(sample.param[1,i]),
              link.fun$K(sample.param[2,i]),
              link.fun$alpha(sample.param[3,i]),
              link.fun$cc(sample.param[4,i]),
              link.fun$pp(sample.param[5,i]),
              link.fun$sigma(sample.param[6,i])
              )
    } else {
      th_ = c(link.fun$mu(sample.param[1,i]),
              link.fun$K(sample.param[2,i]),
              link.fun$alpha(sample.param[3,i]),
              link.fun$cc(sample.param[4,i]),
              link.fun$pp(sample.param[5,i])
              )
      
    }
    # if(is.null(Sigma_)){
    #   sigma.1 <- link.fun$sigma(sample.param[6,i])
    #   Sigma_ <- matrix(c(sigma.1, 0,
    #                      0, sigma.1), ncol = 2, byrow = TRUE)
    # }
    # if no known points just sample with Ht = NULL
    if(nrow(known.data) == 0){
      ss = withTimeout(expr = sample.ETAS(theta = th_,
                                          beta.p = beta.p, M0 = M0, 
                                          T1 = T1.f, T2 = T2.f,
                                          loglambda.bkg = log.bkg.field_, 
                                          mesh.bkg = mesh_, crs_obj = crs_,
                                          bdy = bdy_, Sigma = Sigma_, Ht = NULL, ncore = ncore),
                       timeout = time.out.sec,
                       onTimeout = 'silent')
    } else{ # else set Ht = known.data
      ss = withTimeout(expr = sample.ETAS(theta = th_, 
                                          beta.p = beta.p, 
                                          M0 = M0, 
                                          T1 = T1.f, T2 = T2.f,
                                          loglambda.bkg = log.bkg.field_, 
                                          mesh.bkg = mesh_, crs_obj = crs_,
                                          bdy = bdy_, Sigma = Sigma_, Ht = known.data, 
                                          ncore = ncore),
                       timeout = time.out.sec,
                       onTimeout = 'silent')
    }
    if(is.null(ss)){
      cat('Time out reached at catalogue ', i, '\n')
      next
    } else {
      # combine generations, arrange by time, add catalogue identifier
      ss = bind_rows(ss) %>%
        arrange(ts) %>%
        mutate(cat.idx = i)
      # merge with imposed data if any
      if(!is.null(imposed.data)){
        imposed.data = imposed.data %>%
          mutate(cat.idx = i,
                 gen = 0)
        ss <- bind_rows(ss, imposed.data)  
      }
    }
    # if it is the first sample initialize the file
    if(i == 1){
      write.table(ss, file = paste0('fore_full/', file.name, '.txt'), 
                  col.names = TRUE, row.names = FALSE)
    } else { # else append the new rows to the file
      write.table(ss, file = paste0('fore_full/', file.name, '.txt'), append = TRUE, 
                  col.names = FALSE, row.names = FALSE)
    }
    # print completeness 
    if( display_ & (i/ncol(sample.param)) %% 0.25 == 0){
      cat('completed : ', i/ncol(sample.param), '\n')
    }
  }
  # return the list of samples
}



 








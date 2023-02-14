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
breaks_exp <- function(tt_, T2_, coef_ = 2, delta.t, N.max = 10){
  if(T2_ - tt_ < delta.t){
    return(c(tt_, T2_))
  }
  else{
    tt_breaks <- tt_ + delta.t*((1 + coef_)^(0:N.max))
    tt_breaks <- tt_breaks[tt_breaks < T2_]
    if(T2_ - tt_breaks[length(tt_breaks)] < delta.t){
      tt_breaks[length(tt_breaks)] = T2_
    }
    if(tt_breaks[length(tt_breaks)] < T2_){
      tt_breaks <- c(tt_breaks, T2_)
    }
    c(tt_,tt_breaks)  
  }
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
  proj4string(int_poly) <- bdy_@proj4string
  if(displaypl){
    ggplot()  + gg(int_poly) + gg(bdy_) + geom_point(aes(xx,yy)) +  
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

space.time.grid <- function(data.point, delta.sp, n.layer., min.edge., 
                            coef.t, delta.t, N.max,
                            T2., bdy., displaygrid = FALSE){
  xx. <- data.point$x
  yy. <- data.point$y
  tt. <- data.point$ts
  idx.p <- data.point$idx.p
  # spatial bins
  if(displaygrid){
    Plot_grid(xx = xx., yy = yy., delta_ = delta.sp, n.layer = n.layer., 
              bdy_ =  bdy., min.edge = min.edge.)
  }
  # find coordinates of polygons : each bin is represented by 5 coordinates
  sp.bins.coords <- Find_sp_grid(xx = xx., yy = yy., delta_ = delta.sp, 
                                 n.layer = n.layer., 
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
  t_b <- breaks_exp(tt_ = tt., T2_ = T2., coef_ = coef.t, delta.t = delta.t, N.max = N.max)
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
  vapply(1:nrow(p.loc), \(idx)
         abs(pmvnorm(lower = c(x.min[idx], y.min[idx]), 
                     upper = c(x.max[idx], y.max[idx]), 
                     mean = p.loc[idx,], sigma = Sigma.l[[idx]],
                     keepAttr = FALSE)), 0)
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


compute.grid <- function(param., list.input){
  
  if(is.list(list.input$Sigma_)){
    Sigma.sp <- list.input$Sigma_[sapply(sp.names, \(bname) which(list.input$gridd$ref_layer == bname)[1])]
    Is.vec <- I.s_df(list.input$space.sel, Sigma.sp)
  } else {
    Is.vec <- I.s_df(list.input$space.sel, list.input$Sigma_)
  }
  
  
  It.vec <- It_df(param_ = param., time.df = list.input$time.sel)
  
  list.ret <- list(It = It.vec[list.input$mapping.t],
                   Is = Is.vec[list.input$mapping.s])#,
                   #space.df = list.input$space.sel %>% 
                    # mutate(Is = Is.vec),
                   #time.df = list.input$time.sel %>%
                    # mutate(It = It.vec))
  return(list.ret)
}

## link functions

# exp copula transformation
exp.t <- function(x, rate){
  bru_forward_transformation(qexp, x, rate)
}
# gamma copula transformation
gamma.t <- function(x, a, b){
  bru_forward_transformation(qgamma, x, a, b)
}
# uniform copula transformation
unif.t <- function(x, a, b){
  bru_forward_transformation(qunif, x, min = a, max = b)
}
# log-gaussian copula transformation
loggaus.t <- function(x, m, s){
  bru_forward_transformation(qlnorm, x, meanlog = m, sdlog = s)
}

# inverse link function
inv.exp.t <- function(x, rate){
  qnorm(pexp(x, rate))
}
# gamma copula transformation
inv.gamma.t <- function(x, a, b){
  qnorm(pgamma(x, a, b))
}
# uniform copula transformation
inv.unif.t <- function(x, a, b){
  qnorm(punif(x, a, b))
}
inv.loggaus.t <- function(x, m, s){
  qnorm(plnorm(x, meanlog = m, sdlog = s))
}


### functions for model fitting
ETAS.fit.spatial <- function(sample.s, 
                               M0, T1, T2, bdy,
                               coef_t = 1, delta.t = 0.1, N_max = 10,
                               delta.sp = 0.1, n.layer.sp = -1, min.edge.sp = 1, 
                               link.functions = list(mu = \(x) x,
                                                     K = \(x) x,
                                                     alpha = \(x) x,
                                                     cc = \(x) x,
                                                     pp = \(x) x,
                                                     sigma = \(x) x), 
                               bru.opt = list(bru_verbose = 3,
                                              bru_max_iter = 50),
                               ncore = 1){
  
  # create different data.frames
  # first for background
  df.0 <- data.frame(counts = 0, exposures = 1, part = 'background')
  
  ### create time bins 
  cat('Start creating spatial grid...', '\n')
  time.g.st <- Sys.time()
  
  df.j <- foreach(idx = 1:nrow(sample.s), .combine = rbind) %do% {
    space.time.grid(data.point = sample.s[idx,], 
                    delta.sp = delta.sp, 
                    n.layer. = n.layer.sp, 
                    min.edge. = min.edge.sp,
                    coef.t = coef_t,
                    delta.t = delta.t,
                    N.max = N_max,
                    T2. = T2,
                    bdy. = bdy, displaygrid = FALSE)
  } 
  df.j$counts = 0
  df.j$exposures = 1
  df.j$part = 'triggered'
  t.names <- unique(df.j$t.ref_layer)
  sp.names <- unique(df.j$ref_layer)
  time.sel <- df.j[sapply(t.names, \(bname) which(df.j$t.ref_layer == bname)[1]),]
  space.sel <- df.j[sapply(sp.names, \(bname) which(df.j$ref_layer == bname)[1]),]
  mapping.t <- match(df.j$t.ref_layer, t.names)
  mapping.s <- match(df.j$ref_layer, sp.names)
  
  
  cat('Finished creating spatial grid, time ', Sys.time() - time.g.st, '\n')   
  
  logLambda.h.inla <- function(th.K, th.alpha, th.c, th.p, th.sigma1, 
                               th.sigma2, list.input, 
                               M0, bdy, ncore_ = ncore){
    theta_ <- c(0, 
                link.functions$K(th.K[1]), 
                link.functions$alpha(th.alpha[1]), 
                link.functions$cc(th.c[1]), 
                link.functions$pp(th.p[1]))

    sigma.1 <- link.functions$sigma(th.sigma1[1])
    sigma.2 <- link.functions$sigma(th.sigma2[1])
    Sigma. <- matrix(c(sigma.1, 0, 0, sigma.2), byrow = TRUE, ncol = 2)
    list.input$Sigma_ = Sigma.
    
    #cat('theta - LogL', theta_, '\n')
    comp.list <- compute.grid(param. = theta_, list.input = list.input)
    
    # cat('Is.na = ', sum(is.na(log1p(comp.list$Is - 1))), "\n")
    # cat('It.na = ', sum(is.na(log1p(comp.list$It - 1))), '\n')
    # cat('Is.inf = ', sum(is.infinite(log1p(comp.list$Is - 1))), "\n")
    # cat('It.inf = ', sum(is.infinite(log1p(comp.list$It - 1))), '\n')
    # 
    theta_[3]*(list.input$df_grid$mags - M0) + log(theta_[2]) + log(comp.list$It) + 
      log(comp.list$Is)
    
  }
  
  
  
  # third is for the sum of the log intensities
  df.s <- data.frame(counts = nrow(sample.s), exposures = 0, part = 'SL')
  
  
  loglambda.inla <- function(th.mu, th.K, th.alpha, th.c, th.p, th.sigma1, th.sigma2,
                             bkg, tt, xx, yy, 
                             th, xh, yh, mh, M0, ncore_ = ncore){
    
    
    th.p.df <- data.frame(th.mu = link.functions$mu(th.mu[1])*bkg, 
                          th.K = link.functions$K(th.K[1]), 
                          th.alpha = link.functions$alpha(th.alpha[1]), 
                          th.c = link.functions$cc(th.c[1]), 
                          th.p = link.functions$pp(th.p[1]))
    #cat('th.mu', sum(th.p.df$th.mu <0), '\n')
    sigma.1 <- link.functions$sigma(th.sigma1[1])
    sigma.2 <- link.functions$sigma(th.sigma2[1])
    
    Sigma. <- matrix(c(sigma.1, 0, 0, sigma.2), byrow = TRUE, ncol = 2)
    outp <- unlist(mclapply(1:length(tt), \(idx) 
                            log(lambda(theta = as.numeric(th.p.df[idx,]), 
                                       tt = tt[idx], 
                                       xx = xx[idx], 
                                       yy = yy[idx], 
                                       th = th, xh = xh, yh = yh, mh = mh, M0 = M0, Sigma = Sigma.)),
                            mc.cores = ncore_))
    #cat('logl', sum(is.na(outp)), '\n')
    mean(outp)
  }
  
  
  
  df.total <- bind_rows(df.0, df.j, df.s)
  list.input. <- list(t.names = t.names, 
                      sp.names = sp.names,
                      space.sel = space.sel,
                      time.sel = time.sel,
                      mapping.t = mapping.t,
                      mapping.s = mapping.s,
                      sample.s = sample.s,
                      idx.bkg = df.total$part == 'background',
                      idx.trig = df.total$part == 'triggered',
                      idx.sl = df.total$part == 'SL',
                      n.total = nrow(df.total),
                      df_grid = df.j)
  
  pred.fun <- function(th.mu, th.K, th.alpha, th.c, th.p, th.sigma1,  th.sigma2,
                       list.input, T1, T2, bdy, M0){
    
    out <- rep(0, list.input$n.total)
    out[list.input$idx.bkg] <- log(link.functions$mu(th.mu[1])) + log(T2 - T1)
    out[list.input$idx.trig] <- logLambda.h.inla(th.K = th.K, 
                                                 th.alpha = th.alpha, 
                                                 th.c = th.c, 
                                                 th.p = th.p, 
                                                 th.sigma1 = th.sigma1,
                                                 th.sigma2 = th.sigma2,
                                                 list.input = list.input,
                                                 M0 = M0, bdy = bdy)
    out[list.input$idx.sl] <- loglambda.inla(th.mu = th.mu, 
                                             th.K = th.K, 
                                             th.alpha = th.alpha, 
                                             th.c = th.c, 
                                             th.p = th.p,
                                             th.sigma1 = th.sigma1,
                                             th.sigma2 = th.sigma2,
                                             bkg = list.input$sample.s$bkg,
                                             tt = list.input$sample.s$ts, 
                                             xx = list.input$sample.s$x, 
                                             yy = list.input$sample.s$y,
                                             th = list.input$sample.s$ts, 
                                             xh = list.input$sample.s$x, 
                                             yh = list.input$sample.s$y, 
                                             mh = list.input$sample.s$mags, 
                                             M0 = M0)
    out
    
  }
  
  
  
  # creating formula for summation part
  form.total <- counts ~ pred.fun(th.mu = th.mu,
                                  th.K = th.K, 
                                  th.alpha = th.alpha,
                                  th.c = th.c, 
                                  th.p = th.p, 
                                  th.sigma1 = th.sigma1,
                                  th.sigma2 = th.sigma2,
                                  list.input = list.input.,
                                  T1 = T1, T2 = T2, bdy = bdy, M0 = M0)
  
  
  cmp.total <- counts ~ -1 + 
    th.mu(1, model = 'linear', mean.linear = 0, prec.linear = 1) + 
    th.K(1, model = 'linear', mean.linear = 0, prec.linear = 1) +
    th.alpha(1, model = 'linear', mean.linear = 0, prec.linear = 1) +
    th.c(1, model = 'linear', mean.linear = 0, prec.linear = 1) +
    th.p(1, model = 'linear', mean.linear = 0, prec.linear = 1) +
    th.sigma1(1, model = 'linear', mean.linear = 0, prec.linear = 1) +
    th.sigma2(1, model = 'linear', mean.linear = 0, prec.linear = 1)
  
  bru(components = cmp.total,
      formula = form.total,
      data = df.total,
      family = 'Poisson',
      options = append(bru.opt, list(E = df.total$exposure)))
}


ETAS.fit.magspace <- function(sample.s, 
                         M0, T1, T2, bdy,
                         coef_t = 1, delta.t = 0.1, N_max = 10,
                         delta.sp = 0.1, n.layer.sp = -1, min.edge.sp = 1, 
                         link.functions = list(mu = \(x) x,
                                               K = \(x) x,
                                               alpha = \(x) x,
                                               cc = \(x) x,
                                               pp = \(x) x,
                                               sigma = \(x) x), 
                         bru.opt = list(bru_verbose = 3,
                                        bru_max_iter = 50),
                         ncore = 1){
  
  # create different data.frames
  # first for background
  df.0 <- data.frame(counts = 0, exposures = 1, part = 'background')
  
  ### create time bins 
  cat('Start creating spatial grid...', '\n')
  time.g.st <- Sys.time()
  
  df.j <- foreach(idx = 1:nrow(sample.s), .combine = rbind) %do% {
    space.time.grid(data.point = sample.s[idx,], 
                    delta.sp = delta.sp, 
                    n.layer. = n.layer.sp, 
                    min.edge. = min.edge.sp,
                    coef.t = coef_t,
                    delta.t = delta.t,
                    N.max = N_max,
                    T2. = T2,
                    bdy. = bdy, displaygrid = FALSE)
  } 
  df.j$counts = 0
  df.j$exposures = 1
  df.j$part = 'triggered'
  t.names <- unique(df.j$t.ref_layer)
  sp.names <- unique(df.j$ref_layer)
  time.sel <- df.j[sapply(t.names, \(bname) which(df.j$t.ref_layer == bname)[1]),]
  space.sel <- df.j[sapply(sp.names, \(bname) which(df.j$ref_layer == bname)[1]),]
  mapping.t <- match(df.j$t.ref_layer, t.names)
  mapping.s <- match(df.j$ref_layer, sp.names)
  
  
  cat('Finished creating spatial grid, time ', Sys.time() - time.g.st, '\n')   
  
  logLambda.h.inla <- function(th.K, th.alpha, th.c, th.p, th.sigma1, 
                               th.sigma2, list.input, 
                               M0, bdy, ncore_ = ncore){
    theta_ <- c(0, 
                link.functions$K(th.K[1]), 
                th.alpha[1], 
                link.functions$cc(th.c[1]), 
                link.functions$pp(th.p[1]))
    
    sigma.1 <- link.functions$sigma(th.sigma1[1])
    sigma.2 <- link.functions$sigma(th.sigma2[1])
    Sigma. <- matrix(c(sigma.1, 0, 0, sigma.2), byrow = TRUE, ncol = 2)
    list.input$Sigma_ = Sigma.
    
    #cat('theta - LogL', theta_, '\n')
    comp.list <- compute.grid(param. = theta_, list.input = list.input)
    
    # cat('Is.na = ', sum(is.na(log1p(comp.list$Is - 1))), "\n")
    # cat('It.na = ', sum(is.na(log1p(comp.list$It - 1))), '\n')
    # cat('Is.inf = ', sum(is.infinite(log1p(comp.list$Is - 1))), "\n")
    # cat('It.inf = ', sum(is.infinite(log1p(comp.list$It - 1))), '\n')
    # 
    theta_[3]*(list.input$df_grid$mags - M0) + log(theta_[2]) + log(comp.list$It) + 
      log(comp.list$Is)
    
  }
  
  
  
  # third is for the sum of the log intensities
  df.s <- data.frame(counts = nrow(sample.s), exposures = 0, part = 'SL')
  
  
  loglambda.inla <- function(th.mu, th.K, th.alpha, th.c, th.p, th.sigma1, th.sigma2,
                             bkg, tt, xx, yy, 
                             th, xh, yh, mh, M0, ncore_ = ncore){
    
    
    th.p.df <- data.frame(th.mu = link.functions$mu(th.mu[1])*bkg, 
                          th.K = link.functions$K(th.K[1]), 
                          th.alpha = th.alpha[1], 
                          th.c = link.functions$cc(th.c[1]), 
                          th.p = link.functions$pp(th.p[1]))
    #cat('th.mu', sum(th.p.df$th.mu <0), '\n')
    sigma.1 <- link.functions$sigma(th.sigma1[1])
    sigma.2 <- link.functions$sigma(th.sigma2[1])
    
    Sigma. <- matrix(c(sigma.1, 0, 0, sigma.2), byrow = TRUE, ncol = 2)
    outp <- unlist(mclapply(1:length(tt), \(idx) 
                            log(lambda(theta = as.numeric(th.p.df[idx,]), 
                                       tt = tt[idx], 
                                       xx = xx[idx], 
                                       yy = yy[idx], 
                                       th = th, xh = xh, yh = yh, mh = mh, M0 = M0, Sigma = Sigma.)),
                            mc.cores = ncore_))
    #cat('logl', sum(is.na(outp)), '\n')
    mean(outp)
  }
  
  
  
  df.total <- bind_rows(df.0, df.j, df.s)
  list.input. <- list(t.names = t.names, 
                      sp.names = sp.names,
                      space.sel = space.sel,
                      time.sel = time.sel,
                      mapping.t = mapping.t,
                      mapping.s = mapping.s,
                      sample.s = sample.s,
                      idx.bkg = df.total$part == 'background',
                      idx.trig = df.total$part == 'triggered',
                      idx.sl = df.total$part == 'SL',
                      n.total = nrow(df.total),
                      df_grid = df.j)
  
  pred.fun <- function(th.mu, th.K, th.alpha, th.c, th.p, th.sigma1,  th.sigma2,
                       list.input, T1, T2, bdy, M0){
    
    out <- rep(0, list.input$n.total)
    out[list.input$idx.bkg] <- log(link.functions$mu(th.mu[1])) + log(T2 - T1)
    out[list.input$idx.trig] <- logLambda.h.inla(th.K = th.K, 
                                                 th.alpha = th.alpha, 
                                                 th.c = th.c, 
                                                 th.p = th.p, 
                                                 th.sigma1 = th.sigma1,
                                                 th.sigma2 = th.sigma2,
                                                 list.input = list.input,
                                                 M0 = M0, bdy = bdy)
    out[list.input$idx.sl] <- loglambda.inla(th.mu = th.mu, 
                                             th.K = th.K, 
                                             th.alpha = th.alpha, 
                                             th.c = th.c, 
                                             th.p = th.p,
                                             th.sigma1 = th.sigma1,
                                             th.sigma2 = th.sigma2,
                                             bkg = list.input$sample.s$bkg,
                                             tt = list.input$sample.s$ts, 
                                             xx = list.input$sample.s$x, 
                                             yy = list.input$sample.s$y,
                                             th = list.input$sample.s$ts, 
                                             xh = list.input$sample.s$x, 
                                             yh = list.input$sample.s$y, 
                                             mh = list.input$sample.s$mags, 
                                             M0 = M0)
    out
    
  }
  
  
  
  # creating formula for summation part
  form.total <- counts ~ pred.fun(th.mu = th.mu,
                                  th.K = th.K, 
                                  th.alpha = th.alpha,
                                  th.c = th.c, 
                                  th.p = th.p, 
                                  th.sigma1 = th.sigma1,
                                  th.sigma2 = th.sigma2,
                                  list.input = list.input.,
                                  T1 = T1, T2 = T2, bdy = bdy, M0 = M0)
  
  
  cmp.total <- counts ~ -1 + 
    th.mu(1, model = 'linear', mean.linear = 0, prec.linear = 1) + 
    th.K(1, model = 'linear', mean.linear = 0, prec.linear = 1) +
    th.alpha(1, model = 'linear', mean.linear = 0, prec.linear = 1) +
    th.c(1, model = 'linear', mean.linear = 0, prec.linear = 1) +
    th.p(1, model = 'linear', mean.linear = 0, prec.linear = 1) +
    th.sigma1(1, model = 'linear', mean.linear = 0, prec.linear = 1) +
    th.sigma2(1, model = 'linear', mean.linear = 0, prec.linear = 1)
  
  bru(components = cmp.total,
      formula = form.total,
      data = df.total,
      family = 'Poisson',
      options = append(bru.opt, list(E = df.total$exposure)))
}


### functions for model fitting
ETAS.fit.isotropic <- function(sample.s, 
                               M0, T1, T2, bdy,
                               coef_t = 1, delta.t = 0.1, N_max = 10,
                               delta.sp = 0.1, n.layer.sp = -1, min.edge.sp = 1, 
                               link.functions = list(mu = \(x) x,
                                                     K = \(x) x,
                                                     alpha = \(x) x,
                                                     cc = \(x) x,
                                                     pp = \(x) x,
                                                     sigma = \(x) x), 
                               bru.opt = list(bru_verbose = 3,
                                              bru_max_iter = 50),
                               ncore = 1){
  
  # create different data.frames
  # first for background
  df.0 <- data.frame(counts = 0, exposures = 1, part = 'background')
  
  ### create time bins 
  cat('Start creating spatial grid...', '\n')
  time.g.st <- Sys.time()
  
  df.j <- foreach(idx = 1:nrow(sample.s), .combine = rbind) %do% {
    space.time.grid(data.point = sample.s[idx,], 
                    delta.sp = delta.sp, 
                    n.layer. = n.layer.sp, 
                    min.edge. = min.edge.sp,
                    coef.t = coef_t,
                    delta.t = delta.t,
                    N.max = N_max,
                    T2. = T2,
                    bdy. = bdy, displaygrid = FALSE)
  } 
  df.j$counts = 0
  df.j$exposures = 1
  df.j$part = 'triggered'
  t.names <- unique(df.j$t.ref_layer)
  sp.names <- unique(df.j$ref_layer)
  time.sel <- df.j[sapply(t.names, \(bname) which(df.j$t.ref_layer == bname)[1]),]
  space.sel <- df.j[sapply(sp.names, \(bname) which(df.j$ref_layer == bname)[1]),]
  mapping.t <- match(df.j$t.ref_layer, t.names)
  mapping.s <- match(df.j$ref_layer, sp.names)
  
  
  cat('Finished creating spatial grid, time ', Sys.time() - time.g.st, '\n')   
  
  logLambda.h.inla <- function(th.K, th.alpha, th.c, th.p, th.sigma, list.input, 
                               M0, bdy, ncore_ = ncore){
    theta_ <- c(0, 
                link.functions$K(th.K[1]), 
                link.functions$alpha(th.alpha[1]), 
                link.functions$cc(th.c[1]), 
                link.functions$pp(th.p[1]))
    
    sigma.1 <- link.functions$sigma(th.sigma[1])
    Sigma. <- matrix(c(sigma.1, 0, 0, sigma.1), byrow = TRUE, ncol = 2)
    list.input$Sigma_ = Sigma.
    
    #cat('theta - LogL', theta_, '\n')
    comp.list <- compute.grid(param. = theta_, list.input = list.input)
    
    # cat('Is.na = ', sum(is.na(log1p(comp.list$Is - 1))), "\n")
    # cat('It.na = ', sum(is.na(log1p(comp.list$It - 1))), '\n')
    # cat('Is.inf = ', sum(is.infinite(log1p(comp.list$Is - 1))), "\n")
    # cat('It.inf = ', sum(is.infinite(log1p(comp.list$It - 1))), '\n')
    # 
    theta_[3]*(list.input$df_grid$mags - M0) + log(theta_[2]) + log(comp.list$It) + 
      log(comp.list$Is)
    
  }
  
  
  
  # third is for the sum of the log intensities
  df.s <- data.frame(counts = nrow(sample.s), exposures = 0, part = 'SL')
  
  
  loglambda.inla <- function(th.mu, th.K, th.alpha, th.c, th.p, th.sigma, bkg, tt, xx, yy, 
                             th, xh, yh, mh, M0, ncore_ = ncore){
    
    
    th.p.df <- data.frame(th.mu = link.functions$mu(th.mu[1])*bkg, 
                          th.K = link.functions$K(th.K[1]), 
                          th.alpha = link.functions$alpha(th.alpha[1]), 
                          th.c = link.functions$cc(th.c[1]), 
                          th.p = link.functions$pp(th.p[1]))
    #cat('th.mu', sum(th.p.df$th.mu <0), '\n')
    sigma.1 <- link.functions$sigma(th.sigma[1])
    
    Sigma. <- matrix(c(sigma.1, 0, 0, sigma.1), byrow = TRUE, ncol = 2)
    outp <- unlist(mclapply(1:length(tt), \(idx) 
                            log(lambda(theta = as.numeric(th.p.df[idx,]), 
                                       tt = tt[idx], 
                                       xx = xx[idx], 
                                       yy = yy[idx], 
                                       th = th, xh = xh, yh = yh, mh = mh, M0 = M0, Sigma = Sigma.)),
                            mc.cores = ncore_))
    #cat('logl', sum(is.na(outp)), '\n')
    mean(outp)
  }
  
  
  
  df.total <- bind_rows(df.0, df.j, df.s)
  list.input. <- list(t.names = t.names, 
                      sp.names = sp.names,
                      space.sel = space.sel,
                      time.sel = time.sel,
                      mapping.t = mapping.t,
                      mapping.s = mapping.s,
                      sample.s = sample.s,
                      idx.bkg = df.total$part == 'background',
                      idx.trig = df.total$part == 'triggered',
                      idx.sl = df.total$part == 'SL',
                      n.total = nrow(df.total),
                      df_grid = df.j)
  
  pred.fun <- function(th.mu, th.K, th.alpha, th.c, th.p, th.sigma, 
                       list.input, T1, T2, bdy, M0){
    
    out <- rep(0, list.input$n.total)
    out[list.input$idx.bkg] <- log(link.functions$mu(th.mu[1])) + log(T2 - T1)
    out[list.input$idx.trig] <- logLambda.h.inla(th.K = th.K, 
                                                 th.alpha = th.alpha, 
                                                 th.c = th.c, 
                                                 th.p = th.p, 
                                                 th.sigma = th.sigma,
                                                 list.input = list.input,
                                                 M0 = M0, bdy = bdy)
    out[list.input$idx.sl] <- loglambda.inla(th.mu = th.mu, 
                                             th.K = th.K, 
                                             th.alpha = th.alpha, 
                                             th.c = th.c, 
                                             th.p = th.p,
                                             th.sigma = th.sigma,
                                             bkg = list.input$sample.s$bkg,
                                             tt = list.input$sample.s$ts, 
                                             xx = list.input$sample.s$x, 
                                             yy = list.input$sample.s$y,
                                             th = list.input$sample.s$ts, 
                                             xh = list.input$sample.s$x, 
                                             yh = list.input$sample.s$y, 
                                             mh = list.input$sample.s$mags, 
                                             M0 = M0)
    out
    
  }
  
  
  
  # creating formula for summation part
  form.total <- counts ~ pred.fun(th.mu = th.mu,
                                  th.K = th.K, 
                                  th.alpha = th.alpha,
                                  th.c = th.c, 
                                  th.p = th.p, 
                                  th.sigma = th.sigma, 
                                  list.input = list.input.,
                                  T1 = T1, T2 = T2, bdy = bdy, M0 = M0)
  
  
  cmp.total <- counts ~ -1 + 
    th.mu(1, model = 'linear', mean.linear = 0, prec.linear = 1) + 
    th.K(1, model = 'linear', mean.linear = 0, prec.linear = 1) +
    th.alpha(1, model = 'linear', mean.linear = 0, prec.linear = 1) +
    th.c(1, model = 'linear', mean.linear = 0, prec.linear = 1) +
    th.p(1, model = 'linear', mean.linear = 0, prec.linear = 1) +
    th.sigma(1, model = 'linear', mean.linear = 0, prec.linear = 1)
  
  bru(components = cmp.total,
      formula = form.total,
      data = df.total,
      family = 'Poisson',
      options = append(bru.opt, list(E = df.total$exposure)))
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
sample.loc <- function(xh, yh, n.ev, bdy, Sigma, crs_obj. = NA){
  
  num <- 0
  # until we placed all the events
  while (num < n.ev){
    # sample from the 2D Gaussian kernel without boundaries
    pts.matrix <- rmvnorm(n.ev, mean = c(xh, yh), sigma = Sigma)
    # transform the sample in SpatialPoints and project
    pts <- data.frame(x = pts.matrix[,1], y = pts.matrix[,2])
    coordinates(pts) <- ~ x + y
    #proj4string(pts) <- proj4string(bdy)
    # discard the ones outside bdy
    #keep.point <- pts$x >= bdy@bbox[1,1] & pts$x <= bdy@bbox[1,2] & 
    #  pts$y >= bdy@bbox[2,1] & pts$y <= bdy@bbox[2,2]
    pts <- crop(pts, bdy)#pts[keep.point,]#
    # merge sampled points
    if(num == 0){
      samp.points <- pts
      num <- length(samp.points)
    } else{
      if(length(pts) > 0){
        samp.points <- rbind(samp.points, pts)
        num <- length(samp.points)
      }
    }
    
  }
  
  ## if we retained a number of points > n.ev, we select n.ev events at random
  kp <- sample(seq(1, length(samp.points), by=1), n.ev, replace=FALSE)
  samp.points <- samp.points[kp,]
  samp.points
}


# sampling triggered events
sample.triggered <- function(theta, beta.p, th, xh, yh, n.ev, M0, T1, T2, bdy, 
                             Sigma, crs_obj., Mc = NULL, mag.distro = 'GR'){
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
    samp.mags <- sample.magnitudes(n = n.ev, 
                                   beta.p = beta.p,
                                   M0 = M0,
                                   Mc = Mc, 
                                   distro = mag.distro)
    # sample locations
    samp.locs <- sample.loc(xh = xh, yh = yh, n.ev = n.ev, 
                            bdy = bdy, Sigma = Sigma, crs_obj. = crs_obj.)
    
    # build output dataset
    # samp.points <- data.frame(ts = samp.ts, mags = samp.mags, x = samp.locs$x,
    #                           y = samp.locs$y)
    # 
    samp.points <- data.frame(ts = samp.ts, mags = samp.mags, x = samp.locs@coords[,1],
                               y = samp.locs@coords[,2])
    # return only the ones with time different from NA (the one with NA are outside the interval T1, T2)
    # even though it should not happen given how we built sample.omori
    output <- samp.points[!is.na(samp.points$ts),]
    if(nrow(output) == 0){
      samp.points <- data.frame(x = 1, y = 1, ts = 1, mags = 1)
      samp.points <- samp.points[-1,]
      return(samp.points) 
    } else {
      return(output)  
    }
    
  }
}


sample.generation <- function(theta, beta.p, Ht, M0, T1, T2, bdy,
                              Sigma, crs_obj., ncore = 1, Mc = NULL, mag.distro = 'GR'){
  
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
                     Sigma_l[[idx]], crs_obj. = crs_obj., Mc = Mc, 
                     mag.distro = mag.distro), mc.cores = ncore)
  
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
  n.points=1000
  
  ## Set up a spatialpoints dataframe for our results
  num <- 0
  
  ## To sample the correct number of points, keep going until the num >= cat_n
  while (num < num_events){
    pts <- spsample(bdy, n.points, "random")#points <- spsample(bdy, n.points, "random")
    ## transform to wgs84 long/lat
    #pts <- spTransform(pts, crs_obj)
    ## next few lines modified directly from sample.lgcp
    proj <- INLA::inla.mesh.project(mesh, pts)
    lambda_ratio <- exp(as.vector(proj$A %*% loglambda) - loglambda_max)
    keep <- proj$ok & (runif(n.points) <= lambda_ratio)
    if(sum(keep) == 0){
      print('zero kept')
      next
    }
    kept <- pts[keep]
    if(num == 0){
      samp.points <- kept
      num <- length(samp.points)
    } else {
      samp.points <- rbind(samp.points, kept)
      num <- length(samp.points)
    }
  }
  
  ## Keep exactly cat_n points, choose these randomly from all of the points we've kept so far
  kp <- sample(seq(1, length(samp.points), by=1), num_events, replace=FALSE)
  samp.points <- samp.points[kp,]
  
  return(samp.points)
}


sample.ETAS <- function(theta, beta.p, M0, T1, T2, bdy, Sigma = NULL, 
                        loglambda.bkg = NULL, mesh.bkg = NULL, 
                        crs_obj = NULL, Ht = NULL, ncore = 1,
                        Unit.Area = TRUE, Mc = NULL, mag.distro = 'GR'){
  
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
                           mags = sample.magnitudes(n = n.bkg, 
                                                    beta.p = beta.p,
                                                    M0 = M0,
                                                    Mc = Mc, 
                                                    distro = mag.distro), 
                           gen = 'background')
    }
    # otherwise it samples using the information provided
    else{
      
      bkg.locs <- suppressWarnings(point_sampler(loglambda = loglambda.bkg, bdy = bdy, 
                                                 mesh = mesh.bkg, num_events = n.bkg,
                                                 crs_obj = crs_obj))
      
      bkg.df <- data.frame(x = data.frame(bkg.locs)$x, 
                           y = data.frame(bkg.locs)$y, 
                           ts = runif(n.bkg, T1, T2), 
                           mags = sample.magnitudes(n = n.bkg, 
                                                    beta.p = beta.p,
                                                    M0 = M0,
                                                    Mc = Mc, 
                                                    distro = mag.distro), 
                           gen = 'background')
    }
  }
  
  # if known events are provided
  if(!is.null(Ht)){
    Ht$gen = 'known'
    if(nrow(bkg.df) > 0){
      Gen.list <- list(bind_rows(bkg.df, Ht))  
    }
    else{
      Gen.list <- list(Ht)
    }
    
  }
  else{
    Gen.list <- list(bkg.df)
  }
  # stop if we have no background events and no events generated from known observations
  if(nrow(Gen.list[[1]]) == 0){
    #print(exp(theta.v[1])*(T2 - T1)*(area(bdy)/1000000))
    #stop('No events generated - increase theta1')
    warning('No event generated or provided')
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
                                   M0, T1, T2, bdy, Sigma, ncore, crs_obj. = crs_obj)
    #print(nrow(triggered))
    # stop the loop if there are no more aftershocks
    if(nrow(triggered) == 0){
      flag = FALSE}
    else{
      # set generations
      triggered$gen = as.character(gen)
      # store new generation
      Gen.list[[gen + 1]] = triggered
      # update generation counter
      gen = gen + 1
    }
  }
  #cat('Building data.frame \n')
  df.ret <- bind_rows(Gen.list)
  df.ret <- df.ret %>% 
    filter(ts >= T1, ts < T2) %>%
    arrange(ts)
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

########################################
# functions for automatic model fitting
########################################


input.file.to.list <- function(input_path){
  con <- file(input_path)
  on.exit(close(con))
  par.raw <- readLines(con)
  for(i in 1:length(par.raw)){
    row.i <- par.raw[[i]]
    if(grepl('start.date', row.i)){
      eval(parse(text = row.i))
    } 
    else if(grepl('end.date', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('magnitude.completeness', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('min.longitude', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('max.longitude', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('min.latitude', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('max.latitude', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('catalog.path', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('catalog.header', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('catalog.sep', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('catalog.skip', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('catalog.colnames', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('proj.ESPG', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('a_mu', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('b_mu', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('a_K', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('b_K', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('a_alpha', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('b_alpha', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('a_c', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('b_c', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('a_p', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('b_p', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('r_sigma', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('th.mu.init', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('th.K.init', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('th.alpha.init', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('th.c.init', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('th.p.init', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('th.sigma.init', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('max_iter', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('max_step', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('num.threads', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('coef.t', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('delta.t', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('Nmax', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('delta.sp', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('n.layer.sp', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('min.edge.sp', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('mesh.max.edge.inner', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('mesh.max.edge.outer', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('mesh.cutoff', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('prior.sigma.value', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('prior.sigma.pr', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('prior.range.value', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('prior.range.pr', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('area.shape.path', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('n.periods', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('period.length', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('start.date.fore', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('magnitude.update', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('output.name', row.i)){
      eval(parse(text = row.i))
    }
  }
  # loading catalog
  catalog <- read.table(catalog.path, 
                        header = catalog.header, 
                        sep = catalog.sep, 
                        skip = catalog.skip)
  if(!catalog.header){
    colnames(catalog) <- catalog.colnames
  }
  
  if(!('Lon' %in% colnames(catalog))){
    stop('Error in the catalog column names: please set the column name of the observed longitudes equal to "time_string"')
  }
  if(!('Lat' %in% colnames(catalog))){
    stop('Error in the catalog column names: please set the column name of the observed latitude equal to "time_string"')
  }
  #
  # catalog preparation
  start.date <- as.POSIXct(start.date)
  end.date <- as.POSIXct(end.date)
  catalog <- catalog %>% 
    mutate(time_date = as.POSIXct( paste0(Year,'-', Mo, '-', Da,' ',
                                          Ho,':', Mi,':', Se),
                                   format = '%Y-%m-%d %H:%M:%OS')) %>%
    filter(time_date >= start.date,
           time_date <= end.date,
           Lon >= min.longitude,
           Lon <= max.longitude,
           Lat >= min.latitude,
           Lat <= max.latitude,
           Mw >= magnitude.completeness) %>%
    mutate(time_diff = as.numeric(difftime(time_date, start.date, units = 'days')))
  cat('Finish loading & preparing catalog', '\n')
  
  # create data.frame for inlabru
  data.bru <- data.frame(ts = catalog$time_diff,
                         magnitudes = catalog$Mw,
                         idx.p = seq_len(nrow(catalog)),
                         x = catalog$Lon,
                         y = catalog$Lat)
  # set up time interval and magnitude of completeness 
  T1 <- 0
  T2 <- ceiling(as.numeric(difftime(end.date, start.date, units = 'days')))
  M0 <- magnitude.completeness
  
  # priors
  link.f <- list(mu = \(x) gamma.t(x, a_mu, b_mu), 
                 K = \(x) loggaus.t(x, a_K, b_K), 
                 alpha = \(x) loggaus.t(x, a_alpha, b_alpha), 
                 cc = \(x) loggaus.t(x, a_c, b_c), 
                 pp = \(x) loggaus.t(x, a_p, b_p),
                 sigma = \(x) exp.t(x, r_sigma))
  
  # initial value
  th.init <- list(th.mu = th.mu.init,
                  th.K = th.K.init,
                  th.alpha = th.alpha.init,
                  th.c = th.c.init,
                  th.p = th.p.init,
                  th.sigma = th.sigma.init)
  
  # options for inlabru 
  if(is.null(max_step)){
    bru.opt.list <- list(bru_verbose = 3, # type of visual output 
                         bru_max_iter = max_iter, # maximum number of iterations
                         num.threads = num.threads,
                         #bru_method = list(max_step = 0.5),
                         inla.mode = 'experimental', # type of inla algorithm
                         bru_initial = th.init) # parameters initial values
  } else {
    bru.opt.list <- list(bru_verbose = 3, # type of visual output 
                         bru_max_iter = max_iter, # maximum number of iterations
                         bru_method = list(max_step = max_step),
                         num.threads = num.threads,
                         inla.mode = 'experimental', # type of inla algorithm
                         bru_initial = th.init) # parameters initial values
    
  }
  # output list
  list(catalog = catalog,
       catalog.bru = data.bru,
       time.int = c(start.date, end.date),
       T12 = c(T1, T2),
       lat.int = c(min.latitude, max.latitude),
       lon.int = c(min.longitude, max.longitude),
       M0 = M0,
       ESPG = proj.ESPG,
       link.functions = link.f,
       bru.opt.list = bru.opt.list,
       coef.t = coef.t,
       delta.t = delta.t,
       Nmax = Nmax,
       delta.sp = delta.sp,
       n.layer.sp = n.layer.sp,
       min.edge.sp = min.edge.sp,
       mesh.max.edge.inner = mesh.max.edge.inner,
       mesh.max.edge.outer = mesh.max.edge.outer,
       mesh.cutoff = mesh.cutoff,
       prior.sigma.value = prior.sigma.value,
       prior.sigma.pr = prior.sigma.pr,
       prior.range.value = prior.range.value,
       prior.range.pr = prior.range.pr,
       area.shape.path = area.shape.path,
       n.periods = n.periods,
       period.length = period.length,
       start.date.fore = start.date.fore,
       magnitude.update = magnitude.update,
       output.name = output.name
  )
}


explorative.plots <- function(list.input, n.bins.hist = 50, point.size = 0.2){
  pl.time.mag <- ggplot(list.input$catalog, aes(x = time_date, y = Mw)) + 
    geom_point() + 
    xlab('date') + 
    ylab('Magnitude')
  pl.time.hist <- ggplot(list.input$catalog, aes(x = time_date)) + 
    geom_histogram(bins = n.bins.hist) + 
    xlab('date') 
  cat.sp <- list.input$catalog.bru
  coordinates(cat.sp) <- c('x', 'y')
  crs.obj <- CRS(SRS_string=paste0('EPSG:', list.input$ESPG))
  proj4string(cat.sp) <- crs.obj
  pl.space <- ggplot() + gg(cat.sp, size = point.size) + coord_equal() + 
    ylab('Latitude') + 
    xlab('Longitude')
  list(time.mag = pl.time.mag,
       time.hist = pl.time.hist,
       space = pl.space)
}

background.model <- function(list.input){
  crs.obj <- CRS('+proj=longlat')#CRS(SRS_string=paste0('EPSG:', list.input$ESPG))
  crs.obj2 <- CRS(SRS_string=paste0('EPSG:', 7794))
  crs.obj_km <- fm_crs_set_lengthunit(crs.obj2, "km")
  data.SP <- list.input$catalog.bru
  coordinates(data.SP) <- c('x', 'y')
  proj4string(data.SP) <- crs.obj
  data.km <- spTransform(data.SP, crs.obj_km)
  #ggplot() + gg(data.km)
  bdy <- square_poly_from_bbox(rbind(list.input$lon.int,
                                     list.input$lat.int), crs.obj, buff = 0.01)
  bdy.km <- spTransform(bdy, crs.obj_km)
  bdy.km <- square_poly_from_bbox(bdy.km@bbox, crs.obj = crs.obj_km, buff = 0)
  # create mesh 
  mesh_ <- inla.mesh.2d(boundary = bdy.km, 
                        max.edge = c(list.input$mesh.max.edge.inner, 
                                     list.input$mesh.max.edge.outer), 
                        cutoff = list.input$mesh.cutoff,
                        crs=crs.obj_km)
  cat('Mesh created - number of nodes:', mesh_$n, '\n')
  #ggplot() + gg(mesh_)
  # spde obj
  spde.model <- inla.spde2.pcmatern(mesh_, 
                                    prior.sigma=c(list.input$prior.sigma.value, 
                                                  list.input$prior.sigma.pr), 
                                    prior.range = c(list.input$prior.range.value, 
                                                    list.input$prior.range.pr))
  
  # create model components
  model.cmp <- coordinates~ Smooth(main = coordinates, model = spde.model) +
    Intercept(1)
  cat('Start model fitting \n')
  start.t <- Sys.time()
  # fit the model
  bru_fit_bkg <- lgcp(components = model.cmp, 
                      data = data.km,  
                      domain = list(coordinates = mesh_), 
                      samplers= bdy.km,  
                      options = list(bru_verbose = 3, 
                                     control.inla = list(strategy = 'simplified_laplace'))
  )
  cat('Finished model fitting - Time :', Sys.time() - start.t, '\n')
  list.input$fit.bkg <- bru_fit_bkg
  cat('Start normalizing background field \n')
  # background rate at the mesh points
  pred.bkg.mesh <- predict(bru_fit_bkg, mesh_, ~ exp(Smooth + Intercept))
  
  # create projection to find value of bkg field at the obs points
  proj <- inla.mesh.project(mesh_, data.km@coords)
  
  # add a bkg column to the data
  data.km$bkg <- as.numeric(proj$A %*% pred.bkg.mesh$mean)
  
  ip <- ipoints(mesh_, samplers = bdy.km)
  proj.ip <- inla.mesh.project(mesh_, ip@coords)
  bkg.ip <- as.numeric(proj.ip$A %*% pred.bkg.mesh$mean)
  # approximate integral
  integ.field <- sum(ip$weight*bkg.ip) 
  cat('background model integral : ', integ.field, '\n')
  cat('observed number of events : ', nrow(data.km), '\n')
  
  list.input$catalog.bru.km <- data.km
  list.input$catalog.bru.km$bkg <- data.km$bkg/integ.field
  list.input$mesh = mesh_
  list.input$bdy.km <- bdy.km
  cat('All finished \n')
  list.input
}

isotropic.ETAS.model <- function(list.input, num.cores = 5){
  data.bru <- data.frame(x = list.input$catalog.bru.km@coords[,1],
                         y = list.input$catalog.bru.km@coords[,2],
                         mags = list.input$catalog.bru.km$magnitudes,
                         idx.p = list.input$catalog.bru.km$idx.p,
                         ts = list.input$catalog.bru.km$ts,
                         bkg = list.input$catalog.bru.km$bkg)
  ETAS.fit.isotropic(sample.s = data.bru, # data.frame representing data points 
                     coef_t = list.input$coef.t, 
                     delta.t = list.input$delta.t,  
                     N_max = list.input$Nmax,
                     delta.sp = list.input$delta.sp,
                     n.layer.sp = list.input$n.layer.sp,
                     min.edge.sp = list.input$min.edge.sp, 
                     link.functions  = list.input$link.functions, 
                     M0 = list.input$M0, 
                     T1 = list.input$T12[1], 
                     T2 = list.input$T12[2], 
                     bdy = list.input$bdy.km, # SpatialPolygon representing study region (has to be squared)
                     bru.opt = list.input$bru.opt.list, # verbose changes what inlabru prints 
                     ncore = num.cores) # set number of cores
}


spatial.ETAS.model <- function(list.input, num.cores = 5){
  data.bru <- data.frame(x = list.input$catalog.bru.km@coords[,1],
                         y = list.input$catalog.bru.km@coords[,2],
                         mags = list.input$catalog.bru.km$magnitudes,
                         idx.p = list.input$catalog.bru.km$idx.p,
                         ts = list.input$catalog.bru.km$ts,
                         bkg = list.input$catalog.bru.km$bkg)
  ETAS.fit.spatial(sample.s = data.bru, # data.frame representing data points 
                     coef_t = list.input$coef.t, 
                     delta.t = list.input$delta.t,  
                     N_max = list.input$Nmax,
                     delta.sp = list.input$delta.sp,
                     n.layer.sp = list.input$n.layer.sp,
                     min.edge.sp = list.input$min.edge.sp, 
                     link.functions  = list.input$link.functions, 
                     M0 = list.input$M0, 
                     T1 = list.input$T12[1], 
                     T2 = list.input$T12[2], 
                     bdy = list.input$bdy.km, # SpatialPolygon representing study region (has to be squared)
                     bru.opt = list.input$bru.opt.list, # verbose changes what inlabru prints 
                     ncore = num.cores) # set number of cores
}


cov.ETAS.model <- function(list.input, num.cores = 5){
  data.bru <- data.frame(x = list.input$catalog.bru.km@coords[,1],
                         y = list.input$catalog.bru.km@coords[,2],
                         mags = list.input$catalog.bru.km$magnitudes,
                         idx.p = list.input$catalog.bru.km$idx.p,
                         ts = list.input$catalog.bru.km$ts,
                         bkg = list.input$catalog.bru.km$bkg)
  ETAS.fit.cov(sample.s = data.bru, # data.frame representing data points 
                   coef_t = list.input$coef.t, 
                   delta.t = list.input$delta.t,  
                   N_max = list.input$Nmax,
                   delta.sp = list.input$delta.sp,
                   n.layer.sp = list.input$n.layer.sp,
                   min.edge.sp = list.input$min.edge.sp, 
                   link.functions  = list.input$link.functions, 
                   M0 = list.input$M0, 
                   T1 = list.input$T12[1], 
                   T2 = list.input$T12[2], 
                   bdy = list.input$bdy.km, # SpatialPolygon representing study region (has to be squared)
                   bru.opt = list.input$bru.opt.list, # verbose changes what inlabru prints 
                   ncore = num.cores) # set number of cores
}


get_posterior <- function(link.functions, model_fit, par.names){
  post.list <- lapply(1:length(model_fit$marginals.fixed), \(x) 
                      data.frame(inla.tmarginal(link.functions[[x]], 
                                                model_fit$marginals.fixed[[x]]),
                                 param = par.names[x]))
  bind_rows(post.list)
}


get_mean <- function(list.input, model.fit){
  c(list.input$link.functions$mu(model.fit$summary.fixed$mean[1]),
    list.input$link.functions$K(model.fit$summary.fixed$mean[2]),
    list.input$link.functions$alpha(model.fit$summary.fixed$mean[3]),
    list.input$link.functions$cc(model.fit$summary.fixed$mean[4]),
    list.input$link.functions$pp(model.fit$summary.fixed$mean[5]),
    list.input$link.functions$sigma(model.fit$summary.fixed$mean[6]))
  
}

get_posterior_sample <- function(list.input, model.fit, n.samp, scale = 'ETAS'){
  if(n.samp > 1000){
    s.list <- lapply(1:ceiling(n.samp/1000), function(x) 
                       
                       generate(model.fit, data.frame(), ~ c(th.mu = th.mu,
                                                    th.K = th.K,
                                                    th.alpha = th.alpha,
                                                    th.c = th.c,
                                                    th.p = th.p,
                                                    th.sigma = th.sigma),
                                n.samples = 1000) )
    s.out <- bind_cols(s.list)
  } else{
    s.out <- generate(model.fit, data.frame(), ~ c(th.mu = th.mu,
                                                   th.K = th.K,
                                                   th.alpha = th.alpha,
                                                   th.c = th.c,
                                                   th.p = th.p,
                                                   th.sigma = th.sigma),
                      n.samples = n.samp)
  }
  if(scale == 'ETAS'){
    s.out <- data.frame(mu = list.input$link.functions$mu(s.out[1,]),
                        K = list.input$link.functions$K(s.out[2,]),
                        alpha = list.input$link.functions$alpha(s.out[3,]),
                        c = list.input$link.functions$cc(s.out[4,]),
                        p = list.input$link.functions$pp(s.out[5,]),
                        sigma = list.input$link.functions$sigma(s.out[6,]))
    
    } else if(scale == 'Internal'){
    
    s.out <- data.frame(mu = s.out[1,],
                        K = s.out[2,],
                        alpha = s.out[3,],
                        c = s.out[4,],
                        p = s.out[5,],
                        sigma = s.out[6,])
    
  } else{
    stop('Unkown scale') 
  }
  return(s.out)
}

produce_single_forecast <- function(post.samples, 
                                    start.date.cat,
                                    start.date.fore,
                                    period.len,
                                    increase.start.fore.sec = 60,
                                    list.input,
                                    T.retro,
                                    crs_obj, Mc = NULL, mag.distro = 'GR'){
  
  crs.obj <- CRS(SRS_string=paste0('EPSG:', 7794))
  crs.obj_km <- fm_crs_set_lengthunit(crs.obj, "km")
  
  start.time.fore <- as.numeric(difftime(as.POSIXct(start.date.fore, format = '%Y-%m-%d %H:%M:%OS') + increase.start.fore.sec, 
                                         as.POSIXct(start.date.cat, format = '%Y-%m-%d %H:%M:%OS'), 
                                         units = 'days'))
  
  df.know <- data.frame(
    list.input$catalog.bru.km[list.input$catalog.bru.km$ts < start.time.fore,])
  print(head(df.know[order(df.know$magnitudes, decreasing = TRUE),], 10))
  cat.list <- vector('list', nrow(post.samples))
  
  for(i in seq_len(nrow(post.samples))){
    Sigma.matrix <- matrix(c(post.samples[i,6], 0,
                             0, post.samples[i,6]), byrow = TRUE, ncol = 2)
    cat('catalog : ', i, '\n')
    #cat('param : ', as.numeric(post.samples[i,]), '\n')
    syn.catalog <- sample.ETAS(theta = as.numeric(post.samples[i,]),
                               beta.p = list.input$beta.p,
                               M0 = list.input$M0,
                               T1 = max(start.time.fore - T.retro, 0),
                               T2 = start.time.fore + period.len, 
                               bdy = list.input$bdy.km,
                               Sigma = Sigma.matrix,
                               loglambda.bkg = list.input$loglambda.bkg.at.mesh,
                               mesh.bkg = list.input$mesh, 
                               crs_obj = crs.obj, 
                               Ht = df.know,
                               ncore = 5, 
                               Unit.Area = TRUE, Mc = Mc, 
                               mag.distro = mag.distro
    )
    if(nrow(syn.catalog) > 0){
      syn.catalog$idx.cat <- i
      cat.list[[i]] <- syn.catalog  
    } else {
      next
    }
  }
  return(cat.list)
}

get_productivity_mag <- function(list.input, model.fit){
  
  data.df <- as.data.frame(list.input$catalog.bru.km)
  
  Lambda.trig <- function(th.mu, th.K, th.alpha, th.c, th.p, th.sigma, list.input){
    theta.par <- c(list.input$link.functions$mu(th.mu),
                   list.input$link.functions$K(th.K),
                   list.input$link.functions$alpha(th.alpha),
                   list.input$link.functions$cc(th.c),
                   list.input$link.functions$pp(th.p),
                   list.input$link.functions$sigma(th.sigma))
    Sigma.matrix <- matrix(c(theta.par[6], 0,
                             0, theta.par[6]), byrow = TRUE, ncol = 2)
    exp(logLambda.h.vec(theta = theta.par, 
                        th = data.df$ts, 
                        xh = data.df$x,
                        yh = data.df$y,
                        mh = data.df$magnitudes, 
                        M0 = list.input$M0, 
                        T1 = list.input$T12[1], 
                        T2 = list.input$T12[2], 
                        bdy = list.input$bdy.km,
                        Sigma = Sigma.matrix) )
  }
  
  preds <- predict(model.fit, data.frame(), ~ Lambda.trig(th.mu, th.K, th.alpha, 
                                                          th.c, th.p, th.sigma, 
                                                          list.input))
  preds$magnitudes <- data.df$magnitudes
  preds
}


convert.forecast <- function(date.string, start.cat.date, M.min, period, area.oef,
                             folder.path= NULL){
  if(period == 'day'){
    period.length = 1
    fore.name = 'day'
  }
  if(period == 'week'){
    period.length = 7
    fore.name = 'week'
  }
  # load list of synthetic catalogs
  fore = load(paste0('fore.', date.string, '.', fore.name, '.Rds'))
  fore.cat.list = get(fore)
  
  # convert string to Date obj
  #date.string.2 <- gsub('T', ' ', date.string)
  start.fore.date <- as.POSIXct(date.string, format = '%Y-%m-%d %H:%M:%OS')
  end.fore.date <- start.fore.date + period.length*24*60*60
  # find events to be selected and convert days differences in dates
  list.idx.sel <- lapply(fore.cat.list, function(x){
    transf.time <- start.cat.date + x$ts*60*60*24
    idx.time.sel = (transf.time > start.fore.date) & (transf.time < end.fore.date)
    idx.mag.sel = x$mags > M.min
    data.frame(idx.sel = idx.time.sel & idx.mag.sel,
               time = transf.time)
  })
  # select simulated events after date.string
  list.fore.sel <- lapply(1:length(fore.cat.list), function(x){
    df.out <- fore.cat.list[[x]][list.idx.sel[[x]]$idx.sel, ]
    df.out$time.date = list.idx.sel[[x]]$time[list.idx.sel[[x]]$idx.sel]
    df.out
  })
  # select events in area.oef region and merge all catalogs in unique data.frame
  final.forecast <- foreach(i = 1:length(list.fore.sel), .combine = rbind) %do% {
    fore.cat <- list.fore.sel[[i]]
    if(is.null(nrow(fore.cat))){
      fore.cat <- data.frame(Lon = 1,
                             Lat = 1,
                             Time = 1,
                             Mag = 1,
                             Idx.cat = 1)
      fore.cat <- fore.cat[-1,]
    }
    else{
      if(nrow(fore.cat) > 0){
        fore.cat.sp <- fore.cat
        coordinates(fore.cat.sp) <- c('x', 'y')
        proj4string(fore.cat.sp) <- crs.obj_km
        fore.cat.sp <- spTransform(fore.cat.sp, proj4string(area.oef))
        fore.cat.sp <- crop(fore.cat.sp, area.oef)
        if(is.null(fore.cat.sp)){
          fore.cat <- data.frame(Lon = 1,
                                 Lat = 1,
                                 Time = 1,
                                 Mag = 1,
                                 Idx.cat = 1)
          fore.cat <- fore.cat[-1,]
        } else {
          fore.cat <- data.frame(Lon = fore.cat.sp@coords[,1],
                                 Lat = fore.cat.sp@coords[,2],
                                 Time = gsub(' ', 'T', fore.cat.sp$time.date),
                                 Mag = fore.cat.sp$mags,
                                 Idx.cat = fore.cat.sp$idx.cat)  
        }
      }
      else{
        fore.cat <- data.frame(Lon = 1,
                               Lat = 1,
                               Time = 1,
                               Mag = 1,
                               Idx.cat = 1)
        fore.cat <- fore.cat[-1,]
      }
      rownames(fore.cat) <- NULL
      fore.cat
    }
  }
  date.name <- gsub(' ', 'T', date.string)
  write.table(final.forecast, 
              file = paste0(folder.path, 'forecast.', date.name, period,'.txt'),
              row.names = FALSE,
              sep = ',')
}





sample.tapGR <- function(n, beta.p, Mc, M0){
  b = beta.p/log(10)
  beta.tapGR <- (2/3)*b
  Mom0 <- 10^( (3/2)*M0 + 9.1) 
  Momc <- 10^( (3/2)*Mc + 9.1) 
  R1 = runif(n, 0, 1)
  R2 = runif(n, 0, 1)
  
  M1 = Mom0*(R1^(-1/beta.tapGR))
  M2 = Mom0 - Momc*log(R2)
  
  moments = pmin(M1 , M2)
  magnitudes = 2/3 * (log10(moments) - 9.1)
  return(magnitudes)
}


sample.magnitudes <- function(n, beta.p, M0, Mc = NULL, distro = 'GR'){
  if(distro == 'GR'){
    magnitudes = rexp(n, beta.p) + M0
  }
  else if(distro == 'Tap-GR'){
    if(is.null(Mc)){
      stop('Corner magnitude (Mc) is missing with no default')
    }
    else{
      magnitudes = sample.tapGR(n = n, beta.p = beta.p, Mc = Mc, M0 = M0)
    }
  } else {
    stop("Unknown distribution please choose between 'GR' and 'Tap-GR")
  }
  return(magnitudes)
}














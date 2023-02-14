library(inlabru)
library(INLA)
library(sp)
library(sf)
library(raster)
library(rgeos)
library(RColorBrewer)
library(tmap)
library(ggplot2)
library(maps)
library(rgdal)
library(dplyr)
library(viridis)
source("utils/inlabru_import_functions.R")
source('utils/code_for_etas_FINAL.R')

#####################################
## LOAD DATA AND EXPLORATIVE PLOTS ##
#####################################

list.input <- input.file.to.list('user_input_Aquila.txt')

exp.plot <- explorative.plots(list.input)
exp.plot$time.mag
exp.plot$time.hist
exp.plot$space

##########################
## BACKGROUND MODEL FIT ##
##########################

list.output.bkg <- background.model(list.input)

# to plot the estimated mean background field spatial variation
mesh.pix <- list.output.bkg$mesh
pix <- pixels(mesh.pix, nx = 50, ny = 50)

# check if reasonable
loglambda.plot <- predict(list.output.bkg$fit.bkg, 
                          pix, 
                          ~ Intercept + Smooth)


ggplot() + gg(loglambda.plot['mean']) + 
  gg(list.output.bkg$catalog.bru.km, size = 0.1) + 
  gg(list.output.bkg$bdy.km, color = 'orange') + 
  scale_fill_viridis()


####################
## ETAS MODEL FIT ##
####################

start.time <- Sys.time()
fit.etas.aquila <- isotropic.ETAS.model(list.input = list.output.bkg)
comp.time <- difftime(Sys.time(), start.time, units = 'mins')
comp.time
save(fit.etas.aquila, file = 'fit.etas.aquila.Rds')
load('fit.etas.aquila.Rds')


# retrieve posterior distributions of parameters
post.df.tot <- get_posterior(list.output.bkg$link.functions, 
                             fit.etas.aquila, c('mu', 'K', 'alpha', 'c', 'p', 'sigma')) %>%
  mutate(case = 'total')
list.output.bkg$link.functions$sigma2 <- list.output.bkg$link.functions$sigma
post.df.tot2 <- get_posterior(list.output.bkg$link.functions, 
                             fit.etas.ama2, c('mu', 'K', 'alpha', 'c', 'p', 'sigma',
                                              'sigma2')) %>%
  mutate(case = 'total')

ggplot(post.df.tot, aes(x,y)) +
  geom_line() + 
  facet_wrap(facets = vars(param), scales = 'free')



# get productivity per different magnitudes
prod.mag.total <- get_productivity_mag(list.input = list.output.bkg, 
                           model.fit = fit.etas.aquila) 
  
ggplot(prod.mag.total, aes(x = magnitudes, y = mean)) + geom_line() + 
  geom_ribbon(aes(x = magnitudes, ymin = q0.025, ymax = q0.975), alpha= 0.2) + 
  geom_hline(yintercept = 1, color = 'red')


################
## SIMULATION ##
################

# simulate parameters value
set.seed(123)
post.samp <- get_posterior_sample(list.input = list.output.bkg,
                                  model.fit = fit.etas.aquila, 
                                  n.samp = 10000,
                                  scale = 'ETAS')
save(post.samp, file = 'post.samp.Aquila.Rds')
load('post.samp.Aquila.Rds')

# things needed to be added to produce forecasts
list.output.bkg$catalog.bru.km$mags <- list.output.bkg$catalog.bru.km$magnitudes
# background rate at the mesh points
list.output.bkg$loglambda.bkg.at.mesh <- predict(list.output.bkg$fit.bkg,
                                                     list.output.bkg$mesh, 
                                                     ~ Intercept + Smooth)$mean
# beta parameter of GR law
list.output.bkg$beta.p <- 1/mean(list.output.bkg$catalog$Mw - list.output.bkg$M0 - 0.05)

# crs objects 
crs.obj2 <- CRS(SRS_string=paste0('EPSG:', 7794))
crs.obj_km <- fm_crs_set_lengthunit(crs.obj2, "km")


# read dates to be forecasts
pred.dates <- read.csv('earthquake_catalogue.csv', sep = ',', 
                       header = TRUE)
time_date <- as.POSIXct(gsub('T', ' ', pred.dates$time_string), 
                        format = '%Y-%m-%d %H:%M')

cat_rows <- head(list.output.bkg$catalog[order(list.output.bkg$catalog$Mw,
                                               decreasing = TRUE),], 9)

pred.dates.corr <- foreach(i = 1:length(time_date), .combine = rbind) %do% {
    cat_rows[which.min(abs(difftime(time_date[i], cat_rows$time_date))),]
}

## prediction are for 1 second after the event 
dates_to_pred <- as.character(pred.dates.corr$time_date + 1)
additional_dates <- c('2009-04-06 00:00:00', '2009-04-07 00:00:00',
                      '2009-04-09 00:00:00', '2009-04-13 00:00:00')

##############
## FORECAST ##
##############

# this is to produce daily forecasts starting from the dates provided
dates_to_pred <- c(additional_dates, dates_to_pred)
for(i in 2:length(dates_to_pred)){
  print(dates_to_pred[i])
  single.forecast.week <- produce_single_forecast(post.samples = post.samp[1:10000,], 
                                                  start.date.cat = list.output.bkg$time.int[1],
                                                  start.date.fore = dates_to_pred[i], 
                                                  increase.start.fore.sec = 0,
                                                  period.len = 1,
                                                  list.input = list.output.bkg, 
                                                  T.retro = 3,
                                                  crs_obj = crs.obj_km,
                                                  Mc = 7,
                                                  mag.distro = 'Tap-GR')
  save(single.forecast.week, file = paste0('fore.',dates_to_pred[i],'.day.TGR.Rds'))
}


###########################
# READ AND WRITE TXT ######
###########################

# area of interest for the forecasts
area.oef <- readOGR('shape_Area/rectangle_for_OEF.shp')

preds_dates <- dates_to_pred#, additional_pred)

# write txt files from forecast lists
for(i in 1:length(preds_dates)){
  print(preds_dates[i])
  convert.forecast(date.string = preds_dates[i],
                   start.cat.date = list.output.bkg$time.int[1],
                   M.min = 3.99,
                   period = 'day',
                   area.oef = area.oef,
                   folder.path = 'daily_fore_aquila/')
}













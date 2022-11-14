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

##########################
## BACKGROUND MODEL FIT ##
##########################

list.output.bkg <- background.model(list.input)

# check if reasonable
loglambda.plot <- predict(list.output.bkg$fit.bkg, 
                          pixels(list.output.bkg$mesh, nx = 50, ny = 50), 
                          ~ Intercept + Smooth)


ggplot() + gg(loglambda.plot['mean']) + 
  gg(list.output.bkg$catalog.bru.km, size = 0.1) + 
  gg(list.output.bkg$bdy.km, color = 'orange') + 
  scale_fill_viridis()


####################
## ETAS MODEL FIT ##
####################

# takes around 30 mins
start.time <- Sys.time()
fit.etas <- isotropic.ETAS.model(list.input = list.output.bkg)
comp.time <- difftime(Sys.time(), start.time, units = 'mins')
comp.time

# retrieve posterior distributions of parameters
post.df.tot <- get_posterior(list.output.bkg$link.functions, 
                             fit.etas, c('mu', 'K', 'alpha', 'c', 'p', 'sigma')) %>%
  mutate(case = 'total')

# check them
ggplot(post.df.tot, aes(x,y,color = case)) +
  geom_line() + 
  facet_wrap(facets = vars(param), scales = 'free')


# get productivity per different magnitudes
prod.mag.total <- get_productivity_mag(list.input = list.output.bkg, 
                           model.fit = fit.etas) 
  
ggplot(prod.mag.total, aes(x = magnitudes, y = mean)) + geom_line() + 
  geom_ribbon(aes(x = magnitudes, ymin = q0.025, ymax = q0.975), alpha= 0.2) + 
  geom_hline(yintercept = 1, color = 'red')


################
## SIMULATION ##
################

# simulate parameters value
post.samp <- get_posterior_sample(list.input = list.output.bkg,
                                  model.fit = fit.etas, 
                                  n.samp = 10000)

# things needed to be added to produce forecasts
# this will be soon removed
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

##############################################################################
## FIRST FORECAST ############################################################
##############################################################################


single.forecast <- produce_single_forecast(post.samples = post.samp[1:10000,], 
                                                    start.date.cat = list.output.bkg$time.int[1],
                                                    start.date.fore = '2009-04-06 01:32:00', 
                                                    increase.start.fore.sec = 60,
                                                    period.len = 1,
                                                    list.input = list.output.bkg, 
                                                    T.retro = 5,
                                                    crs_obj = crs.obj_km)
save(single.forecast, file = 'fore.2009-04-06T01:32:00.day.Rds')

## weekly
single.forecast.week <- produce_single_forecast(post.samples = post.samp[1:10000,], 
                                           start.date.cat = list.output.bkg$time.int[1],
                                           start.date.fore = '2009-04-06 01:32:00', 
                                           increase.start.fore.sec = 60,
                                           period.len = 7,
                                           list.input = list.output.bkg, 
                                           T.retro = 5,
                                           crs_obj = crs.obj_km)
save(single.forecast.week, file = 'fore.2009-04-06T01:32:00.week.Rds')
N.sim.week <- vapply(single.forecast.week, nrow, 0)
quantile(N.sim.week, c(0.025, 0.975))
summary(N.sim.week)


###############################################################################
## SECOND FORECAST ############################################################
###############################################################################

single.forecast <- produce_single_forecast(post.samples = post.samp[1:10000,], 
                                           start.date.cat = list.output.bkg$time.int[1],
                                           start.date.fore = '2009-04-06 02:37:00', 
                                           increase.start.fore.sec = 60,
                                           period.len = 1,
                                           list.input = list.output.bkg, 
                                           T.retro = 5,
                                           crs_obj = crs.obj_km)
save(single.forecast, file = 'fore.2009-04-06T02:37:00.day.Rds')
N.sim <- vapply(single.forecast, nrow, 0)
quantile(N.sim, c(0.025, 0.975))
summary(N.sim)
ggplot(single.forecast[[which.max(N.sim)]], aes(x, y)) + geom_point() + gg(list.output.bkg$bdy.km) +
  coord_equal()


## weekly
single.forecast.week <- produce_single_forecast(post.samples = post.samp[1:10000,], 
                                                start.date.cat = list.output.bkg$time.int[1],
                                                start.date.fore = '2009-04-06 02:37:00', 
                                                increase.start.fore.sec = 60,
                                                period.len = 7,
                                                list.input = list.output.bkg, 
                                                T.retro = 5,
                                                crs_obj = crs.obj_km)
save(single.forecast.week, file = 'fore.2009-04-06T02:37:00.week.Rds')




###############################################################################
## Third FORECAST ############################################################
###############################################################################

single.forecast <- produce_single_forecast(post.samples = post.samp[1:10000,], 
                                           start.date.cat = list.output.bkg$time.int[1],
                                           start.date.fore = '2009-04-06 23:15:00', 
                                           increase.start.fore.sec = 60,
                                           period.len = 1,
                                           list.input = list.output.bkg, 
                                           T.retro = 5,
                                           crs_obj = crs.obj_km)
single.forecast2 <- produce_single_forecast(post.samples = post.samp[5888:10000,], 
                                            start.date.cat = list.output.bkg$time.int[1],
                                            start.date.fore = '2009-04-06 23:15:00', 
                                            increase.start.fore.sec = 60,
                                            period.len = 1,
                                            list.input = list.output.bkg, 
                                            T.retro = 5,
                                            crs_obj = crs.obj_km)
save(single.forecast2, file = 'fore.2009-04-06T23:15:00_2.day.Rds')
N.sim <- vapply(single.forecast, nrow, 0)
quantile(N.sim, c(0.025, 0.975))
summary(N.sim)
ggplot(single.forecast[[which.max(N.sim)]], aes(x, y)) + geom_point() + gg(list.output.bkg$bdy.km) +
  coord_equal()


## weekly
single.forecast.week <- produce_single_forecast(post.samples = post.samp[1:10000,], 
                                                start.date.cat = list.output.bkg$time.int[1],
                                                start.date.fore = '2009-04-06 23:15:00', 
                                                increase.start.fore.sec = 60,
                                                period.len = 7,
                                                list.input = list.output.bkg, 
                                                T.retro = 5,
                                                crs_obj = crs.obj_km)
save(single.forecast.week, file = 'fore.2009-04-06T23:15:00.week.Rds')

###############################################################################
## Fourth FORECAST ############################################################
###############################################################################

single.forecast <- produce_single_forecast(post.samples = post.samp[1:10000,], 
                                           start.date.cat = list.output.bkg$time.int[1],
                                           start.date.fore = '2009-04-07 09:26:00', 
                                           increase.start.fore.sec = 60,
                                           period.len = 1,
                                           list.input = list.output.bkg, 
                                           T.retro = 5,
                                           crs_obj = crs.obj_km)
save(single.forecast, file = 'fore.2009-04-07T09:26:00.day.Rds')
N.sim <- vapply(single.forecast, nrow, 0)
quantile(N.sim, c(0.025, 0.975))
summary(N.sim)
ggplot(single.forecast[[which.max(N.sim)]], aes(x, y)) + geom_point() + gg(list.output.bkg$bdy.km) +
  coord_equal()


## weekly
single.forecast.week <- produce_single_forecast(post.samples = post.samp[1:10000,], 
                                                start.date.cat = list.output.bkg$time.int[1],
                                                start.date.fore = '2009-04-07 09:26:00', 
                                                increase.start.fore.sec = 60,
                                                period.len = 7,
                                                list.input = list.output.bkg, 
                                                T.retro = 5,
                                                crs_obj = crs.obj_km)
save(single.forecast.week, file = 'fore.2009-04-07T09:26:00.week.Rds')


###############################################################################
## Fifth FORECAST ############################################################
###############################################################################

## weekly
single.forecast.week <- produce_single_forecast(post.samples = post.samp[1:10000,], 
                                                start.date.cat = list.output.bkg$time.int[1],
                                                start.date.fore = '2009-04-07 17:47:00', 
                                                increase.start.fore.sec = 60,
                                                period.len = 7,
                                                list.input = list.output.bkg, 
                                                T.retro = 5,
                                                crs_obj = crs.obj_km)
save(single.forecast.week, file = 'fore.2009-04-07T17:47:00.week.Rds')

e4409368f4dd478aaee57754ebbedc51d5b8af71
76c1c8bf4ab4d29477aafadfdca3a7f7dd493424

###############################################################################
## Sixth FORECAST ############################################################
###############################################################################

## weekly
single.forecast.week <- produce_single_forecast(post.samples = post.samp[1:10000,], 
                                                start.date.cat = list.output.bkg$time.int[1],
                                                start.date.fore = '2009-04-09 00:52:00', 
                                                increase.start.fore.sec = 60,
                                                period.len = 7,
                                                list.input = list.output.bkg, 
                                                T.retro = 5,
                                                crs_obj = crs.obj_km)
save(single.forecast.week, file = 'fore.2009-04-09T00:52:00.week.Rds')



###########################
# READ AND WRITE TXT ######
###########################
# area of interest for the forecasts
area.oef <- readOGR('shape_Area/rectangle_for_OEF.shp')

date.forecasted <- c('2009-04-06T01:32:00', '2009-04-06T02:37:00',
                     '2009-04-06T23:15:00', '2009-04-07T09:26:00')

# rerun to remove rownames.
for(i in 1:length(date.forecasted)){
  print(date.forecasted[i])
  convert.forecast(date.string =  date.forecasted[i],
                   period = 'day',
                   area.oef = area.oef)
  
  convert.forecast(date.string =  date.forecasted[i],
                   period = 'week',
                   area.oef = area.oef)
}




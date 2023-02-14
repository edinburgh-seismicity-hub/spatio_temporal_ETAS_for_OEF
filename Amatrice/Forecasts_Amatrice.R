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

list.input <- input.file.to.list('user_input_Amatrice.txt')

exp.plot <- explorative.plots(list.input)
exp.plot$time.mag
exp.plot$time.hist
exp.plot$space
nrow(list.input$catalog)

##########################
## BACKGROUND MODEL FIT ##
##########################
#source('utils/code_for_ETAS_server.R')

list.output.bkg <- background.model(list.input)

# debugonce(pixels)
mesh.pix <- list.output.bkg$mesh
# mesh.pix$crs <- NULL
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
#qlnorm(c(0.025, 0.5, 0.975), meanlog = -4, sdlog = 0.5)
list.output.bkg$link.functions$K <- \(x) {unif.t(x, 0, 20)}
# takes around 30 mins
debugonce(ETAS.fit.isotropic)

start.time <- Sys.time()
fit.etas.ama <- isotropic.ETAS.model(list.input = list.output.bkg)
comp.time <- difftime(Sys.time(), start.time, units = 'mins')
comp.time
save(fit.etas.ama, file = 'fit.etas.Amatrice.Rds')
load('fit.etas.Amatrice.Rds')

start.time <- Sys.time()
fit.etas.ama2 <- spatial.ETAS.model(list.input = list.output.bkg)
comp.time <- difftime(Sys.time(), start.time, units = 'mins')
comp.time
save(fit.etas.ama2, file = 'fit.etas.Amatrice2.Rds')

start.time <- Sys.time()
fit.etas.ama3 <- cov.ETAS.model(list.input = list.output.bkg)
comp.time <- difftime(Sys.time(), start.time, units = 'mins')
comp.time
save(fit.etas.ama3, file = 'fit.etas.Amatrice3.Rds')


#load('fit.etas.Rds')
# retrieve posterior distributions of parameters
post.df.tot <- get_posterior(list.output.bkg$link.functions, 
                             fit.etas.ama, c('mu', 'K', 'alpha', 'c', 'p', 'sigma')) %>%
  mutate(case = 'total')
list.output.bkg$link.functions$sigma2 <- list.output.bkg$link.functions$sigma
post.df.tot2 <- get_posterior(list.output.bkg$link.functions, 
                             fit.etas.ama2, c('mu', 'K', 'alpha', 'c', 'p', 'sigma',
                                              'sigma2')) %>%
  mutate(case = 'total')


K.med <- list.output.bkg$link.functions$K(fit.etas.ama$summary.fixed$`0.5quant`[2])
c.med <-list.output.bkg$link.functions$cc(fit.etas.ama$summary.fixed$`0.5quant`[4])
p.med <-list.output.bkg$link.functions$pp(fit.etas.ama$summary.fixed$`0.5quant`[5])
 
K.med*c.med/(p.med - 1)

post.df.tot$case <- 'iso'
post.df.tot2$case <- 'double'

# check them
ggplot(bind_rows(post.df.tot2,
                 post.df.tot), aes(x,y,color = case)) +
  geom_line() + 
  facet_wrap(facets = vars(param), scales = 'free')

ggplot(post.df.tot, aes(x,y)) +
  geom_line() + 
  facet_wrap(facets = vars(param), scales = 'free')


# model comparison?
# introduce correlation in ETAS space triggering
# ----- needs to take corr only for M > threshold
# ------- this needs Sigma taken as list and functions calculated with lapply



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
set.seed(123)
post.samp <- get_posterior_sample(list.input = list.output.bkg,
                                  model.fit = fit.etas.ama, 
                                  n.samp = 10000,
                                  scale = 'Internal')

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


#pred.dates <- read.csv('earthquake_catalogue.csv')
#time_string <- gsub('T', ' ', pred.dates$time_string)
#time_date <- as.POSIXct(time_string)

pred.dates <- read.csv('data/Amatrice_Norcia_earthquakes.csv', sep = ',', 
                       header = TRUE)
time_date <- as.POSIXct(pred.dates$ITACA_event_time, 
                        format = '%d/%m/%Y %H:%M')

cat_rows <- head(list.output.bkg$catalog[order(list.output.bkg$catalog$Mw,
                                               decreasing = TRUE),], 9)

pred.dates.corr <- foreach(i = 1:9, .combine = rbind) %do% {
    cat_rows[which.min(abs(difftime(time_date[i], cat_rows$time_date))),]
}

## prediction are for 1 second after the event 
dates_to_pred <- as.character(pred.dates.corr$time_date + 1)

##############################################################################
## FIRST FORECAST ############################################################
##############################################################################

source('utils/code_for_etas_FINAL.R')

## weekly
additional_pred <- c("2016-08-24 00:00:00", "2016-10-26 00:00:00", "2016-10-30 00:00:00", "2017-01-18 00:00:00")
# dates_to_pred
dates_to_pred <- c(additional_pred, dates_to_pred)
for(i in 9:length(dates_to_pred)){
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



start.date <- '2009-04-06 00:00:00'
convert.forecast(date.string = start.date,
                 start.cat.date = list.output.bkg$time.int[1],
                 M.min = 3.99,
                 period = 'week',
                 area.oef = area.oef)
###########################
# READ AND WRITE TXT ######
###########################
# area of interest for the forecasts
area.oef <- readOGR('shape_Area/rectangle_for_OEF.shp')

preds_dates <- c(dates_to_pred)#, additional_pred)
debugonce(convert.forecast)
# write txt files from forecast lists

for(i in 9:length(dates_to_pred)){
  print(dates_to_pred[i])
  convert.forecast(date.string = dates_to_pred[i],
                   start.cat.date = list.output.bkg$time.int[1],
                   M.min = 3.99,
                   period = 'day',
                   area.oef = area.oef,
                   folder.path = 'ufficial_daily_amatrice/')
}

aa <- read.table('ufficial_daily_amatrice/forecast.2016-10-26T19:18:08day.txt',
                 header = T, sep = ',')

hist(aa$Mag, breaks = seq(3.99, 8, by = 0.2))

NN <- vapply(1:10000, function(x) sum(aa$Idx.cat == x), 0)

plot(table(NN)/10000, xlim = c(0,10))




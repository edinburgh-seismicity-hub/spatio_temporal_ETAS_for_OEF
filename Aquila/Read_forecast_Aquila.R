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

list.input <- input.file.to.list('Aquila/user_input_Aquila.txt')
list.output.bkg <- background.model(list.input)

## read dates of prediction
pred.dates <- read.csv('Aquila/earthquake_catalogue.csv', sep = ',', 
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
dates_to_pred <- c(dates_to_pred, additional_dates)

# set up area
area.oef <- readOGR('shape_Area/rectangle_for_OEF.shp')
inner.bdy <- square_poly_from_bbox(matrix(c(12.9, 13.75, 42.1, 43.25), byrow = TRUE, ncol = 2),
                                   area.oef@proj4string, 0)
inner.bdy.km <- spTransform(inner.bdy, list.output.bkg$bdy.km@proj4string)

pix.mesh <- pixels(list.output.bkg$mesh, nx = 200, ny = 200, mask = inner.bdy.km) 
pix.mesh$ID = 1:length(pix.mesh)

fore.to.pix.list <- foreach(date.pred = dates_to_pred) %do% {
  print(date.pred)
  forecast.to.pix.df(start.fore.string = date.pred, 
                     pix.mesh = pix.mesh, 
                     list.input = list.output.bkg,
                     fore.path = 'Aquila/Forecasts/')
}  

summary.table <- bind_rows(
  lapply(1:length(fore.to.pix.list), \(idx) 
         data.frame(start.date = dates_to_pred[idx],
                    obs = sum(fore.to.pix.list[[idx]]$count),
                    mean = sum(fore.to.pix.list[[idx]]$mean),
                    q0.025 = sum(fore.to.pix.list[[idx]]$q0.025),
                    q0.5 = sum(fore.to.pix.list[[idx]]$q0.5),
                    q0.975 = sum(fore.to.pix.list[[idx]]$q0.975))) )

summary.table[order(as.POSIXct(summary.table$start.date)),]

fore.df <- bind_rows(fore.to.pix.list)
ggplot(fore.df, aes(x,y,fill = log(q0.5))) + 
  geom_tile() + 
  scale_fill_viridis() + 
  xlim(list.output.bkg$bdy.km@bbox[1,]) +
  ylim(list.output.bkg$bdy.km@bbox[2,]) + 
  facet_wrap(facets = vars(date))

fore.df.for.grid <- rbind(data.frame(x = fore.df$x,
                                     y = fore.df$y,
                                     date.fore = fore.df$date,
                                     value = fore.df$count,
                                     quantity = 'observed'),
                          data.frame(x = fore.df$x,
                                     y = fore.df$y,
                                     date.fore = fore.df$date,
                                     value = fore.df$mean,
                                     quantity = 'mean'),
                          data.frame(x = fore.df$x,
                                     y = fore.df$y,
                                     date.fore = fore.df$date,
                                     value = fore.df$q0.025,
                                     quantity = 'q0.025'),
                          data.frame(x = fore.df$x,
                                     y = fore.df$y,
                                     date.fore = fore.df$date,
                                     value = fore.df$q0.5,
                                     quantity = 'median'),
                          data.frame(x = fore.df$x,
                                     y = fore.df$y,
                                     date.fore = fore.df$date,
                                     value = fore.df$q0.975,
                                     quantity = 'q0.975'))

fore.df.for.grid$quantity <- factor(fore.df.for.grid$quantity,
                                    levels = c('observed', 'mean',
                                               'q0.025', 'median','q0.975'))

obs.list <- foreach(date.pred = dates_to_pred) %do% {
  print(date.pred)
  obs.to.forecast(start.fore.string = date.pred, 
                  list.input = list.output.bkg)
}  

library(ggstar)
obs.df <- bind_rows(obs.list)
#obs.df$date.fore <- as.POSIXct(obs.df$date.fore)
pdf('fore.grid.Aquila.pdf', width = 480/25, height = 480/35)
ggplot() + 
  geom_tile(data = fore.df.for.grid, mapping = aes(x,y,fill = log(value))) + 
  geom_star(data = obs.df[obs.df$magnitudes > 5,],
             mapping = aes(x,y), fill = 'red', size = 2) + 
  scale_fill_viridis() + 
  xlim(inner.bdy.km@bbox[1,]) +
  ylim(inner.bdy.km@bbox[2,]) +
  labs(fill = 'log(N)') + 
  coord_equal() + 
  xlab('Easting') + 
  ylab('Northing') + 
  facet_grid(vars(quantity), vars(date.fore), switch = 'y') +
  theme_bw()
dev.off()



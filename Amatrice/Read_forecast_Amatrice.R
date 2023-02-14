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

list.input <- input.file.to.list('Amatrice/user_input_Amatrice.txt')
list.output.bkg <- background.model(list.input)

## read dates of prediction
pred.dates <- read.csv('Amatrice/Amatrice_Norcia_earthquakes.csv', sep = ',', 
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

## prediction are for 1 second after the event 
dates_to_pred <- as.character(pred.dates.corr$time_date + 1)
additional_dates <- c("2016-08-24 00:00:00", "2016-10-26 00:00:00", 
                      "2016-10-30 00:00:00", "2017-01-18 00:00:00")
dates_to_pred <- c(dates_to_pred, additional_dates)

# set up area
area.oef <- readOGR('shape_Area/rectangle_for_OEF.shp')
inner.bdy <- square_poly_from_bbox(matrix(c(12.9, 13.75, 42.1, 43.25), byrow = TRUE, ncol = 2),
                                   area.oef@proj4string, 0)
inner.bdy.km <- spTransform(inner.bdy, list.output.bkg$bdy.km@proj4string)

pix.mesh <- pixels(list.output.bkg$mesh, nx = 200, ny = 200, mask = inner.bdy.km) 
pix.mesh$ID = 1:length(pix.mesh)

forecast.to.pix.df <- function(start.fore.string, pix.mesh, list.input,
                               fore.path){
  start.fore <- as.POSIXct(start.fore.string, format = '%Y-%m-%d %H:%M:%OS')
  date.start.cat <- as.POSIXct(list.input$time.int[1])
  
  start.fore.days <- as.numeric(difftime(start.fore, date.start.cat, units = 'days'))
  
  # calculate end date
  end.fore <- start.fore + 60*60*24
  # select observation in fore period
  idx.obs <- list.input$catalog$time_date > start.fore & 
    list.input$catalog$time_date <= end.fore
  obs.sp <- list.input$catalog.bru.km[idx.obs,]
  # load corresponding forecast
  load(paste0(fore.path,'fore.', start.fore.string,'.day.Rds'))
  single.forecast <- single.forecast.week
  single.forecast <- lapply(single.forecast, \(x) 
                     x[x$ts > start.fore.days & 
                         x$ts <= (start.fore.days+1), ])
  #return(nrow(single.forecast[[1]]))
  #calculate observed number of events per pixel
  over_obs <- over(obs.sp, pix.mesh, byid = TRUE)
  table.pix <- table(over_obs$ID)
  pix.mesh$count <- 0
  pix.mesh$count[as.numeric(names(table.pix))] = table.pix
  pix.mesh$logcount <- log(pix.mesh$count)
  # calculate number of events per pixel per synthetic catalogue
  fore.count.matrix = matrix(0, nrow = length(pix.mesh$ID), ncol = length(single.forecast))
  counter.zeros = 0
  for(j in 1:length(single.forecast)){
    single.forecast.cat <- single.forecast[[j]]
    if(is.null(single.forecast.cat)){
      counter.zeros = counter.zeros + 1
      next
    }
    if(nrow(single.forecast.cat) == 0){
      counter.zeros = counter.zeros + 1
      next
    }
    coordinates(single.forecast.cat) <- c('x', 'y')
    proj4string(single.forecast.cat) <- pix.mesh@proj4string
    single.forecast.cat <- spTransform(single.forecast.cat, pix.mesh@proj4string)
    over_fore <- over(single.forecast.cat, pix.mesh, byid = TRUE)
    table.pix.fore <- table(over_fore$ID)
    fore.count.matrix[as.numeric(names(table.pix.fore)),j] = table.pix.fore
  }
  if(counter.zeros == length(single.forecast)){
    print('no events')
    output <- data.frame(x = 0,
                         y = 0,
                         count = 0,
                         mean = 0,
                         sd = 0,
                         q0.025 = 0,
                         q0.5 = 0,
                         q0.975 = 0,
                         date = 0)
    return(output[-1,])
  } else{
    return(data.frame(x = pix.mesh@coords[,1],
                      y = pix.mesh@coords[,2],
                      count = pix.mesh$count,
                      mean = apply(fore.count.matrix, 1, mean),
                      sd = apply(fore.count.matrix, 1, sd),
                      q0.025 =  apply(fore.count.matrix, 1, \(x) quantile(x, 0.025)),
                      q0.5 = apply(fore.count.matrix, 1, \(x) quantile(x, 0.5)),
                      q0.975 = apply(fore.count.matrix, 1, \(x)quantile(x, 0.975)),
                      date = as.character(start.fore)))
  }
}


N.test.stat <- function(start.fore.string, list.input,
                               fore.path){
  start.fore <- as.POSIXct(start.fore.string, format = '%Y-%m-%d %H:%M:%OS')
  date.start.cat <- as.POSIXct(list.input$time.int[1])
  
  start.fore.days <- as.numeric(difftime(start.fore, date.start.cat, units = 'days'))
  
  # calculate end date
  end.fore <- start.fore + 60*60*24
  # select observation in fore period
  idx.obs <- list.input$catalog$time_date > start.fore & 
    list.input$catalog$time_date <= end.fore
  obs.sp <- list.input$catalog.bru.km[idx.obs,]
  # load corresponding forecast
  load(paste0(fore.path,'fore.', start.fore.string,'.day.Rds'))
  single.forecast <- single.forecast.week
  single.forecast <- lapply(single.forecast, \(x) 
                            x[x$ts > start.fore.days & 
                                x$ts <= (start.fore.days+1), ])
  single.forecast.df <- bind_rows(single.forecast)
  N.fore <- vapply(1:10000, \(x) sum(single.forecast.df$idx.cat == x),0)
  return(data.frame(date.fore = as.character(start.fore),
                    prob.over = mean(N.fore >= nrow(obs.sp)),
                    prob.under = mean(N.fore <= nrow(obs.sp))))
}


obs.to.forecast <- function(start.fore.string, list.input){
  start.fore <- as.POSIXct(start.fore.string, format = '%Y-%m-%d %H:%M:%OS')
  start.fore.m1 <- start.fore - 60*60*24
  
  # select observation in fore period
  idx.obs <- list.input$catalog$time_date >= start.fore.m1 & 
    list.input$catalog$time_date < start.fore
  if(sum(idx.obs) == 0){
    obs.sp <- list.input$catalog.bru.km[1,]
    obs.sp$date.fore <- as.character(start.fore)
    obs.sp <- obs.sp[-1,]
  } else{
    obs.sp <- list.input$catalog.bru.km[idx.obs,]
    obs.sp$date.fore <- as.character(start.fore)
  }
  # load corresponding forecast
  return(as.data.frame(obs.sp))
}


fore.to.pix.list <- foreach(date.pred = dates_to_pred) %do% {
  print(date.pred)
  forecast.to.pix.df(start.fore.string = date.pred, 
                     pix.mesh = pix.mesh, 
                     list.input = list.output.bkg,
                     fore.path = 'Amatrice/Forecasts/')
}  

N.test.stat.df <- foreach(date.pred = dates_to_pred, .combine = rbind) %do% {
  print(date.pred)
  N.test.stat(start.fore.string = date.pred, 
                     list.input = list.output.bkg,
                     fore.path = 'Amatrice/Forecasts/')
}  

N.test.stat.df

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
pdf('fore.grid.ama.pdf', width = 480/20, height = 480/20)
ggplot(fore.df.for.grid, aes(x,y,fill = log(value))) + 
  geom_tile() + 
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

obs.list <- foreach(date.pred = dates_to_pred) %do% {
  print(date.pred)
  obs.to.forecast(start.fore.string = date.pred, 
                  list.input = list.output.bkg)
}  
library(ggstar)
obs.df <- bind_rows(obs.list)

obs.df[obs.df$date.fore =="2016-10-26 00:00:00",]

plot(list.output.bkg$catalog$time_date, list.output.bkg$catalog$Mw, 
     xlim = as.POSIXct(c("2016-10-23 00:00:00", "2016-10-26 00:00:00")))

#obs.df$date.fore <- as.POSIXct(obs.df$date.fore)
pdf('fore.grid.Amatrice.pdf', width = 480/25, height = 480/35)
ggplot() + 
  geom_tile(data = fore.df.for.grid, mapping = aes(x,y,fill = log(value))) + 
  geom_star(data = obs.df[obs.df$magnitudes > 5,],
             mapping = aes(x,y), fill = 'red') + 
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






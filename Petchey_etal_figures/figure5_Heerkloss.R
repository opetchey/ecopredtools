## Code for figure 5 of Petchey et al
## Graphing the forecasts horizon relationship from Beninca et al 2008
## with body size and number of trophic link data
## Owen Petchey 1.12.14

rm(list=ls())
library(reshape2)
library(ggplot2)
library(stringr)
library(mgcv)
library(dplyr)

## read in and tidy the data from figure 2 of Beninca et al 2008
dd <- read.csv("~/Dropbox (Dept of Geography)/1. petchey EFH/ecopredtools/Petchey_etal_figures/data/fig5.r-squared.csv")
dd <- melt(dd, id.vars="prediction_time")
dd <- cbind(dd, do.call(rbind, str_split(dd$variable, "\\.")))
names(dd) <- c("Prediction.time", "dummy1", "Forecast.proficiency", "Group", "Model")
str(dd)
summary(dd)
## have a look at the data
ggplot(dd, aes(x=Prediction.time, y=Forecast.proficiency, col=Model)) +
  facet_wrap(~Group) + geom_point() + stat_smooth()


## we'll fit a polynomial to each data group, and use this sp to get the forecast horizon
dd.s <- split(dd, list(dd$Model, dd$Group))
rez <- list()
pred <- list()
p.over <- seq(3.35, 40, 0.1)
for(i in 1:length(dd.s)) {
  rez[[i]] <- lm(Forecast.proficiency ~ poly(Prediction.time,4), dd.s[[i]])
  pred[[i]] <- predict(rez[[i]], newdata=data.frame(Prediction.time=p.over),interval="confidence")
  pred[[i]]
}
pred[[1]]
rr <- as.data.frame(do.call("rbind", pred))
rr$Prediction.time <- rep(p.over, length(rr[,1])/length(p.over))
rr$Model.Group <- rep(names(dd.s), each=length(p.over))
rr <- cbind(rr, do.call("rbind", str_split(rr$Model.Group, "\\.")))
names(rr)[6:7] <- c("Model", "Group")
str(rr)

## Check the polynomial fits
ggplot(rr, aes(x=Prediction.time, y=fit, col=Model)) +
  facet_wrap(~Group) + 
  geom_line() +
  geom_line(aes(y=upr)) + geom_line(aes(y=lwr)) + 
  geom_point(data=dd, aes(y=Forecast.proficiency))

## a function to get the forecast horizon
get.fh <- function(y, x, threshold)
  min(x[y<threshold])
  
## get the mean, upper, and lower forecast horizons
mean.fh <- lapply(pred, function(x) get.fh(x[,"fit"], p.over,0.6))      
lwr.fh <- lapply(pred, function(x) get.fh(x[,"lwr"], p.over,0.6))      
upr.fh <- lapply(pred, function(x) get.fh(x[,"upr"], p.over,0.6))      

## put the data together
fh <- data.frame(Model.Group=names(dd.s),
                 mean=do.call(rbind, mean.fh), 
                 lwr=do.call(rbind, lwr.fh), 
                 upr=do.call(rbind, upr.fh))
fh <- cbind(fh, do.call("rbind", str_split(fh$Model.Group, "\\.")))
names(fh)[5:6] <- c("Model", "Group")

## read in the body sizes and the number of trophic links
ee <- read.csv("~/Dropbox (Dept of Geography)/1. petchey EFH/ecopredtools/Petchey_etal_figures/data/fig5.groups_real_Heerkloss.csv")[,1:7]
str(ee)

## get cell volume assuming a sphere
ee$size <- ee$biovol.ug.fresh

## calculate average size and number of links
ee <- group_by(ee, short, long) %>% summarise(mean.size=10^mean(log10(size), na.rm=T), mean.link=mean(num.links, na.rm=T))

## merge the two data sets
zz <- merge(ee, fh, by.x="short", by.y="Group")
## if this merge fails, shutdown R, restart it, and try again.


## select only the relevant variables and model
zz.oo <- subset(zz, long!="nitrogen" & long!="phosphorus" & long!="calanoids copepods" & Model=="nonlinear")
zz.oo$log.mean <- log10(zz.oo$mean.size)

## plot the size - forecast horizon relationship
quartz(width=4, height=3.5)
ggplot(zz.oo, aes(x=mean.size, y=mean, col=long)) + theme_bw() +
  geom_point(size=5) +
  geom_errorbar(aes(ymax=upr, ymin=lwr), width=0.0) +
  scale_x_log10() +
  xlab("Organism size") +
  ylab("Forecast horizon") +
  scale_color_discrete(name="Group") + theme(legend.position="none")

## plot the number of links - forecast horizon relationship
quartz(width=5.75, height=3.5)
ggplot(zz.oo, aes(x=mean.link, y=mean, col=long)) +
  geom_point(size=5) +
  geom_errorbar(aes(ymax=upr, ymin=lwr), width=0.0) +
  xlab("Number of trophic links") +
  ylab("Forecast horizon") +
  scale_color_discrete(name="Group")  + theme_bw() + theme(legend.key = element_rect(colour = "white"))


## and check the stats of the relationships
library(nlme)
m1 <- lm(mean ~ log.mean + mean.link, zz.oo)
summary(m1)

m2 <- lm(mean ~ log.mean * mean.link, zz.oo)
summary(m2)

m3 <- gls(mean ~ log.mean + mean.link, data=zz.oo, weights=varPower(), control=list(maxIter=1000, tolerance=0.01))
summary(m3)

m4 <- gls(mean ~  log.mean * mean.link, data=zz.oo, weights=varPower(), control=list(maxIter=1000, tolerance=0.01))
summary(m4)

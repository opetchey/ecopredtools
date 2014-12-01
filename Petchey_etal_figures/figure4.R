## Code for figure 4 of Petchey et al
## Graphing the results of the community dynamic model
## Model was performed in Matlab by Mikael Pontarp
## Owen Petchey 1.12.14

rm(list=ls())
library(plyr)
library(stringr)
library(reshape2)
library(ggplot2)
library(Hmisc)

## Read in and clean the mean forecast horizons
dd <- read.csv("~/Dropbox (Dept of Geography)/1. prediction concept/EFHtools/Petchey_etal_figures/data/PH_outputYdata.csv")
dd <- melt(dd)
dd <- cbind(dd, do.call("rbind", strsplit(as.character(dd$variable), "\\.")))
names(dd) <- c("Uncertainty", "Junk", "Forecast.horizon", "Variable", "Evolution")
str(dd)
dd$num.uncert <- rep(1:4, 6)

## Read in and clean the sd forecast horizons
ee <- read.csv("~/Dropbox (Dept of Geography)/1. prediction concept/EFHtools/Petchey_etal_figures/data/PH_outputERRORdata.csv")
ee <- melt(ee)
ee <- cbind(ee, do.call("rbind", strsplit(as.character(ee$variable), "\\.")))
names(ee) <- c("Uncertainty", "Junk", "Forecast.horizon", "Variable", "Evolution")
str(ee)
ee$num.uncert <- rep(1:4, 6)

## put the errors in the same data frame as the means
dd$error <- ee$Forecast.horizon
dd$upper <- dd$Forecast.horizon + dd$error
dd$lower <- dd$Forecast.horizon - dd$error

## set dodge value for ggplot
pd <- position_dodge(0.3)

## plotting
ggplot(dd, aes(x=num.uncert, y=Forecast.horizon, col=Variable, linetype=Evolution)) +
	geom_point(size=3, position=pd) +
	geom_line(position=pd) +
  	geom_errorbar(aes(ymax=upper, ymin=lower), width=0.1, position=pd) +
  	xlab("Uncertainty in environment conditions\n(1=low, 4=high)") +
  	ylab("Forecast horizon") +
  	scale_color_discrete(name="Variable",
                         breaks=c("extrate", "singlepop", "totalbiomass"),
                         labels=c("Extinction rate", "Population abundance", "Total biomass")) +
  	scale_linetype_discrete(name="Evolution",
                         breaks=c("evo", "noevo"),
                         labels=c("Yes", "No")) +
    theme_bw() + theme(legend.key = element_rect(colour = "white"))











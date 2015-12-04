## Code for figure 4 of Petchey et al
## Graphing the results of the community dynamic model
## Model was performed in Matlab by Mikael Pontarp
## Owen Petchey 1.12.14

rm(list=ls())
library(plyr)
library(dplyr)
library(stringr)
library(reshape2)
library(ggplot2)
library(Hmisc)
library(RCurl)

## Read in and clean the mean forecast horizons
dd <- read.csv(text=getURL("https://raw.githubusercontent.com/opetchey/ecopredtools/master/Petchey_etal_figures/data/fig4.PH_outputYdata.csv"))
dd <- melt(dd)
dd <- cbind(dd, do.call("rbind", strsplit(as.character(dd$variable), "\\.")))
names(dd) <- c("Uncertainty", "Junk", "Forecast.horizon", "Variable", "Evolution")
str(dd)
dd$num.uncert <- rep(1:4, 6)

## Read in and clean the sd forecast horizons
ee <- read.csv(text=getURL("https://raw.githubusercontent.com/opetchey/ecopredtools/master/Petchey_etal_figures/data/fig4.PH_outputERRORdata.csv"))
ee <- melt(ee)
ee <- cbind(ee, do.call("rbind", strsplit(as.character(ee$variable), "\\.")))
names(ee) <- c("Uncertainty", "Junk", "Forecast.horizon", "Variable", "Evolution")
str(ee)
ee$num.uncert <- rep(1:4, 6)

## put the errors in the same data frame as the means
dd$error <- ee$Forecast.horizon
dd$upper <- dd$Forecast.horizon + dd$error
dd$lower <- dd$Forecast.horizon - dd$error

## don't use the extinction rates
dd <- filter(dd, Variable!="extrate")

## set dodge value for ggplot
pd <- position_dodge(0.3)

## plotting
p <- ggplot(dd, aes(x=num.uncert, y=Forecast.horizon, col=Variable, linetype=Evolution)) +
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
p




############
## Below some code to give figures containing only parts of the data


## get the colours used:
ggplot_build(p)$data

## only evolution off or on
ggplot(filter(dd, Evolution=="evo"), aes(x=num.uncert, y=Forecast.horizon, col=Variable, linetype=Evolution)) +
  geom_point(size=3, position=pd) +
  geom_line(position=pd, linetype="dashed") +
  geom_errorbar(aes(ymax=upper, ymin=lower), width=0.1, position=pd) +
  xlab("Uncertainty in environment conditions\n(1=low, 4=high)") +
  ylab("Forecast horizon") +
  ylim(0, 100) +
  scale_color_discrete(name="Variable",
                       breaks=c("extrate", "singlepop", "totalbiomass"),
                       labels=c("Extinction rate", "Population abundance", "Total biomass")) +
  scale_linetype_manual(name="Evolution",
                          breaks=c("evo", "noevo"),
                          labels=c("Yes", "No"),
                          values=c("dashed")) +
  theme_bw() + theme(legend.key = element_rect(colour = "white"))

## only population abundance
ggplot(filter(dd, Variable=="singlepop"), aes(x=num.uncert, y=Forecast.horizon, col=Variable, linetype=Evolution)) +
  geom_point(size=3, position=pd) +
  geom_line(position=pd) +
  geom_errorbar(aes(ymax=upper, ymin=lower), width=0.1, position=pd) +
  xlab("Uncertainty in environment conditions\n(1=low, 4=high)") +
  ylab("Forecast horizon") +
  ylim(0, 100) +
  scale_color_discrete(name="Variable",
                       breaks=c("extrate", "singlepop", "totalbiomass"),
                       labels=c("Extinction rate", "Population abundance", "Total biomass")) +
  scale_linetype_manual(name="Evolution",
                        breaks=c("evo", "noevo"),
                        labels=c("Yes", "No"),
                        values=c("solid","dashed")) +
  theme_bw() + theme(legend.key = element_rect(colour = "white"))

## only total biomass
ggplot(filter(dd, Variable=="totalbiomass"), aes(x=num.uncert, y=Forecast.horizon, col=Variable, linetype=Evolution)) +
  geom_point(size=3, position=pd) +
  geom_line(position=pd) +
  geom_errorbar(aes(ymax=upper, ymin=lower), width=0.1, position=pd) +
  xlab("Uncertainty in environment conditions\n(1=low, 4=high)") +
  ylab("Forecast horizon") +
  ylim(0, 100) +
  scale_color_manual(name="Variable",
                       breaks=c("extrate", "singlepop", "totalbiomass"),
                       labels=c("Extinction rate", "Population abundance", "Total biomass"), 
                       values=c("#00BFC4", "red", "blue")) +
  scale_linetype_manual(name="Evolution",
                        breaks=c("evo", "noevo"),
                        labels=c("Yes", "No"),
                        values=c("solid","dashed")) +
  theme_bw() + theme(legend.key = element_rect(colour = "white"))







## Code for figure 3 of Petchey et al
## Owen Petchey 16.8.2014

## The same methods are used as in "figure2.r" for getting prediction skill.
## Slightly different method for getting prediction horizon (for computational speed reasons)
## Prediction horizons are then plotted against uncertainty

## Preliminaries
rm(list=ls())
library(zoo)
library(ggplot2)
library(Hmisc)
library(dplyr)
set.seed(101)

## The Ricker model and a function to iterate it
ricker  <- function(N, r, demo.stoch) {
  if(!demo.stoch)
    Nn <- N*exp(r*(1-N))
  if(demo.stoch) {
    exp.N <- N*exp(r*(1-N))
    Nn <- exp.N + rnorm(1, exp.N, sd=exp.N*0.000000)
  }
  Nn
}
iterate.ricker <- function(r, N, its, demo.stoch)
{
  Ns <- numeric(length(its)+1)
  Ns[1] <- N
  for(i in 2:its)
      Ns[i] <- ricker(Ns[i-1], r, demo.stoch)
  Ns
} 


## Set the parameters of the numerical "experiment"
## Distribution from which to choose "real" value of r
r.real.mean <- 2.9 
r.real.sd <- 0

## Distribution from which to choose "real" value of N0
N0.real.mean <- 0.8
N0.real.sd <- 0

## uncertainty in r for prediction
r.pred.sd <- c(0, 10^c(-5,-3,-1))
## uncertainty in N0 for prediction
N0.pred.sd <-  c(0, 10^c(-5,-3,-1))

pred.CV <- c(0.00001, 0.001, 0.1)
r.pred.sd <- pred.CV*r.real.mean
N0.pred.sd <- pred.CV*N0.real.mean

## switch for demographic stochasticity
demo.stoch <- c(F, T)

## replicate predicted time series
reps <- 1:1000

## set up experiment
expt <- expand.grid(r.real.mean=r.real.mean,
                    r.real.sd=r.real.sd,
                    N0.real.mean=N0.real.mean,
                    N0.real.sd=N0.real.sd,
                    r.pred.sd=r.pred.sd,
                    N0.pred.sd=N0.pred.sd,
                    demo.stoch=demo.stoch,
                    reps=reps)

## Get the real values of r and N0 (this is a bit redundant when real.r.sd and real.N0.sd = 0)                    
expt$r.real <- rnorm(length(expt[,1]), mean=expt$r.real.mean, sd=expt$r.real.sd)
expt$N0.real <- rnorm(length(expt[,1]), mean=expt$N0.real.mean, sd=expt$N0.real.sd)

## Get values of r and N0 to use for predictions
expt$r.pred <- rnorm(length(expt[,1]), mean=expt$r.real, sd=expt$r.pred.sd)
expt$N0.pred <- rnorm(length(expt[,1]), mean=expt$N0.real, sd=expt$N0.pred.sd)

## Check the experiment
str(expt)
##ggplot(expt, aes(x=N0.pred)) + geom_density()   + facet_grid(N0.pred.sd ~ r.pred.sd)

## Here, the prediction horizon is defined as when the prediction skill measure
## has been below the prediction skill threshold for a certain number of
## consecutive times

## set the prediction skill threshold
pred.skill.threshold <- 0.3
## set the prediction skill threshold duration
## i.e., the number of consecutive windows that prediction skill must
## be below prediction skill threshold
threshold.duration <- 0.5 ## in number of moving window widths
## set the maximum duration of simulations
max.its <- 128
## set the moving window size for calculating prediction skill
mov.wind.width <- 20

## run the next line only to test within the get.pred.horizon() function
##pars <- expt[1,]

## Function to get prediction horizon for each experimental replicate
get.pred.horizon <- function(pars, pred.skill.threshold, max.its, mov.wind.width, threshold.duration) {
  
  ## Some storage
  N.real <- rep(NA, length=max.its+1)
  N.sim <- rep(NA, length=max.its+1)
  N.real[1] <- pars$N0.real
  N.sim[1] <- pars$N0.pred
  ma.cor <- rep(NA, length=max.its+1-mov.wind.width+1)
  
  ## iterate Ricker just enough for one test of if prediction horizon is met
  for(i in 2:(mov.wind.width*(threshold.duration+1))) {
    N.real[i] <- ricker(N.real[i-1], pars$r.real, pars$demo.stoch)
    N.sim[i] <- ricker(N.sim[i-1], pars$r.pred, F)
  }
  ## and get the prediction skills for this initial part of the simulation
  for(j in 2:(mov.wind.width*threshold.duration+1))
    ma.cor[j] <- cor(N.real[j:(j+mov.wind.width-1)], N.sim[j:(j+mov.wind.width-1)])
  
  ## Now continue doing this until either the prediction horizon or max.its is reached
  tt=2
  while(sum(ma.cor[tt:(tt+mov.wind.width*threshold.duration-1)]>pred.skill.threshold)!=0 & i<max.its)
    {
  
    i <- i+1
    tt <- tt+1
    j <- j+1
    N.real[i] <- ricker(N.real[i-1], pars$r.real, pars$demo.stoch)
    N.sim[i] <- ricker(N.sim[i-1], pars$r.pred, F)
    ma.cor[j] <- cor(N.real[j:(j+mov.wind.width-1)], N.sim[j:(j+mov.wind.width-1)])
  }
  ## Return the prediction horizon
  tt
}

## Storage
expt$pred.horizon <- rep(NA, length(expt[,1]))
## Run the experiment
for(i in 1:length(expt[,1]))
  expt$pred.horizon[i] <- get.pred.horizon(expt[i,], pred.skill.threshold, max.its, mov.wind.width, threshold.duration)

save.image(file="~/Dropbox (Dept of Geography)/1. Petchey EFH/ecopredtools/Petchey_etal_figures/data/fig3.Rdata")


## Only run from here once a dataset is saved
rm(list=ls())
library(ggplot2)
library(scales)
library(Hmisc)
library(dplyr)
load('~/Dropbox (Dept of Geography)/1. Petchey EFH/ecopredtools/Petchey_etal_figures/data/data.fig3.Rdata')
#load("/Users/Frank/Documents/My scientific articles/2015 - Prediction horizons/ecopredtools/Petchey_etal_figures/data/fig3.Rdata")

## get average etc prediction horizons by treatments
CI <- .10
aa <- group_by(expt, N0.pred.sd, r.pred.sd, demo.stoch) %>%
  summarise(mean.pred.horizon=mean(pred.horizon),
            sd.pred.horizon=sd(pred.horizon),
            median.pred.horizon=median(pred.horizon),
            upper=smedian.hilow(pred.horizon, conf.int=CI)[3],
            lower=smedian.hilow(pred.horizon, conf.int=CI)[2])

aa$N0.pred.CV <- aa$N0.pred.sd/N0.real.mean
aa$r.pred.CV <- aa$r.pred.sd/r.real.mean

## amount to dodge by
pd <- position_dodge(0.3)

## plot medians
ggplot(aa, aes(x=N0.pred.CV, y=median.pred.horizon,
                      col=as.factor(r.pred.CV), linetype=as.factor(demo.stoch))) +
  geom_line(position=pd, alpha=0.75) +
  geom_point(size=3, position=pd, alpha=0.75) +
  labs(linetype="Demographic stochasticity",
       col="Uncertainty in r: CV(r)", x="Uncertainty in N0: CV(N0)", y="Forecast horizon") +
  ylim(c(0, 30)) +
  scale_x_continuous(breaks=c(0.1, 0.001, 0.00001), trans="log10", label=comma)+
  geom_errorbar(aes(ymax=upper, ymin=lower), width=0.8, position=pd, alpha=0.75) +
  theme_bw() + theme(legend.key = element_rect(colour = "white"))


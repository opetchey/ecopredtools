## Code for figure 2 of Petchey et al
## Making a graph showing how prediction skill decays into the future
## Discrete time Ricker model
## Owen Petchey 16.8.2014

## The prediction skill is calculated by making a "real" time series from a "real" value of r and N0.
## Then a predicted time series is made, but with a value of r and N0 with some error (uncertainty).
## Prediction skill is a moving window measure of the match between the real and predicted time series.
## For each real time series, there are many predicted, so average prediction skill can be calculated.
## This can all be repeated for specified levels of uncertainty in r and N0 used for prediction.


## Preliminaries
rm(list=ls())
library(zoo)
library(ggplot2)
library(Hmisc)
library(dplyr)
set.seed(100)


## The Ricker model and a function to iterate it
ricker  <- function(N,r) N*exp(r*(1-N)) ## copied from ecolMod package!
iterate.ricker <- function(r, N, its, demo.stoch=F)
{
  Ns <- numeric(length(its)+1)
  Ns[1] <- N
  for(i in 2:its) {
    if(!demo.stoch)
      Ns[i] <- ricker(Ns[i-1], r)
    if(demo.stoch) {
      exp.N <- ricker(Ns[i-1], r)
      Ns[i] <- exp.N + rnorm(1, exp.N, sd=exp.N*0.01)
    }
  }
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
r.pred.sd <- c(0, 0.001, 0.02)
## uncertainty in N0 for prediction
N0.pred.sd <-  c(0, 0.001, 0.02)

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

## Choose how long to run the simulations for
its <- 50

## Make the real dynanimcs
real.dyn <- lapply(1:length(expt[,1]),
                   function(x) iterate.ricker(r=expt[x, "r.real"],
                                              N=expt[x, "N0.real"],
                                              its=its,
                                              demo.stoch=expt[x, "demo.stoch"]))

## Make the real dynanimcs
pred.dyn <- lapply(1:length(expt[,1]),
                   function(x) iterate.ricker(r=expt[x, "r.pred"],
                                              N=expt[x, "N0.pred"],
                                              its=its,
                                              demo.stoch=expt[x, "demo.stoch"]))

## Choose the width of the moving window for calculating prediction skill
ma.width <- 3

## Some ugly code that will let us merge the expt info with the prediction skills
expt.long  <- expt[rep(row.names(expt), each=its-ma.width+1),]
expt.long$its <- rep(1:(its-ma.width+1), length(expt.long[,1])/(its-ma.width+1))

## Calculate the prediction skills
## Here it is the correlation between the real and predicted abundances in the window
mov.cor <- lapply(1:length(expt[,1]),
                  function(xx) rollapply(cbind(real.dyn[[xx]], pred.dyn[[xx]]),
                                         width=ma.width, function(x) cor(x)[1,2], by.column=F))

## Some housekeeping
mov.cor <- stack(as.data.frame(mov.cor))[,1]
rr <- data.frame(expt.long, mov.cor)

## A function to get empirical confidence limits
emp.CL <- function(x, percent) sort(x)[percent/100*length(x)]

## Set the prediction skill threshold... used to get the prediction horizon
pred.skill.threshold <- 0.5																

## get the average etc of the prediction skills across the replicates
rr1 <- group_by(rr, N0.pred.sd, r.pred.sd, demo.stoch, its) %>%
  summarise(median.pred.skill=median(mov.cor),
            mean.pred.skill=mean(mov.cor),
            upper.CL=emp.CL(mov.cor, 55),
            lower.CL=emp.CL(mov.cor, 45))


## Get the prediction horizon for the averaged prediction skill
rr2 <- group_by(rr1, N0.pred.sd, r.pred.sd, demo.stoch) %>%
  summarise(pred.horizon=min(its[mean.pred.skill<pred.skill.threshold]-1),
            last.pred.skill=mean.pred.skill[min(its[mean.pred.skill<pred.skill.threshold]-1)])
rr3 <- rr2[rep(rownames(rr2), each=2),]
rr3$last.pred.skill[seq(1,length(rr3$last.pred.skill), 2)] <- -0.2

save.image(file="~/Dropbox (Dept of Geography)/1. prediction concept/det_chaos/data1.Rdata")


## Only run from here once you have a dataset saved

rm(list=ls())
library(ggplot2)
load("~/Dropbox (Dept of Geography)/1. prediction concept/EFHtools/Petchey_etal_figures/data/data1.Rdata")
## Plot the loss of prediction skill through time
rr1$nice.ds <- ifelse(rr1$demo.stoch, "With demographic stochasticity", "Without demographic Stochasticity")  
rr3$nice.ds <- ifelse(rr3$demo.stoch, "With demographic stochasticity", "Without demographic Stochasticity")  
g <- ggplot(data=rr1, aes(x=its, y=mean.pred.skill, col=as.factor(N0.pred.sd), linetype=as.factor(r.pred.sd))) +
  geom_line(size=0.7) +
  geom_point(size=2) +
  labs(col="sd(N0)", linetype="sd(r)", x="Time (generations)", y="Forecast proficiency") +
  facet_grid(.~nice.ds) + 
  geom_hline(yintercept=pred.skill.threshold, col="purple", linetype=2, alpha=0.5)
g + geom_line(data=rr3, aes(x=pred.horizon,
                            y=last.pred.skill,
                            col=as.factor(N0.pred.sd),
                            linetype=as.factor(r.pred.sd))) +
  theme_bw() + theme(legend.key = element_rect(colour = "white"))


## Plot the distributions of prediction horizons
# rr4 <- group_by(rr, N0.pred.sd, r.pred.sd, demo.stoch, reps) %>%
#   summarise(pred.horizon=min(its[mov.cor<pred.skill.threshold]))
# rr4$nice.r.pred.sd <- paste("sd(r) =", rr4$r.pred.sd)	
# rr4$nice.ds <- ifelse(rr4$demo.stoch, "With demog. stoch.", "Without demog. stoch.")  
# ggplot(rr4, aes(x=pred.horizon, fill=as.factor(N0.pred.sd))) +
#   facet_grid(nice.ds~nice.r.pred.sd) +
#   geom_density(alpha=0.5, position="identity", adjust=5) +  
#   labs(x="Prediction horizon", y="Frequency", fill="sd(N0)")



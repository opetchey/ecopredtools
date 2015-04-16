## Code for figure 2b, Petchey et al. 2015 Ecology Letters.

## This figure and the underlying analysis were suggested by a reviewer.
## The reviewer's idea was to include a case study including
## a directional change in environental
## conditions. The reviewer gave some R code to illustrate; this code is below.
## Owen Petchey 27.2.15

### Start of reviewers code
#### Ricker
ricker <- function(N, r, temp)
  {
  K <- K_mean + 40*temp
  N_t <- N*exp(r*(1-(N/K)))
  return(N_t)
}

N1 <- 100
r <- 2.05
K_mean <- 100
tlimit <- 500
N <- numeric(tlimit)
N[1] <- N1
temp <- c(rnorm(tlimit/2, 12, 1), rnorm(tlimit/2, 18, 1))
for(t in 2:tlimit){
  N[t] <- ricker(N[t-1], r, temp[t])
}
par(mfrow=c(2,1))
plot(c(1:tlimit), N, type="l", xlab="time")
abline(v = 250, col="dodgerblue", lwd=4)
plot(c(1:tlimit), temp, type="l", xlab="time", ylab="environmental driver")
abline(v = 250, col="dodgerblue", lwd=4)
### End of reviewers code


### Start of Owen's code
## We'll use the code from figure 2a as a basis
rm(list=ls())
library(zoo)
library(ggplot2)
library(Hmisc)
library(dplyr)
set.seed(100)

## We will make K change through time, starting at K_start,
## and increasing per time step by K_step.
## Not a particularly elegant solution.

## The Ricker model and a function to iterate it
ricker  <- function(N, r, K)
	N*exp(r*(1-N/K)) ## copied from ecolMod package!
iterate.ricker <- function(r, N, K_step, its, demo.stoch=F)
{
  Ns <- numeric(length(its)+1)
  Ns[1] <- N
  K <- K_start
  for(i in 2:its) {
    K <- K + K_step
    if(K<0) stop("K less than zero")
    if(!demo.stoch)
      Ns[i] <- ricker(Ns[i-1], r, K)
    if(demo.stoch) {
      exp.N <- ricker(Ns[i-1], r, K)
      Ns[i] <- exp.N + rnorm(1, exp.N, sd=exp.N*0.01)
    }
  }
  Ns
} 

## Set the parameters of the numerical "experiment"

## Distribution from which to choose "real" value of K_step
K_step.real.mean <- c(10)
K_step.real.sd <- 0

## Distribution from which to choose "real" value of r
r.real.mean <- 2.9 
r.real.sd <- 0

## Distribution from which to choose "real" value of N0
N0.real.mean <- 0.8
N0.real.sd <- 0

## use standard deviation instead of sd
pred.CV <- c(0, 0.0005, 0.010)
r.pred.sd <- 0   #pred.CV*r.real.mean
N0.pred.sd <- pred.CV*N0.real.mean
K_step.pred.sd <- pred.CV*K_step.real.mean

## switch for demographic stochasticity
demo.stoch <- c(F)

## replicate predicted time series
reps <- 1:1000

## set up experiment
expt <- expand.grid(r.real.mean=r.real.mean,
                    r.real.sd=r.real.sd,
                    N0.real.mean=N0.real.mean,
                    N0.real.sd=N0.real.sd,
                    K_step.real.mean=K_step.real.mean,
                    K_step.real.sd=K_step.real.sd,                    
                    r.pred.sd=r.pred.sd,
                    N0.pred.sd=N0.pred.sd,
                    K_step.pred.sd=K_step.pred.sd,
                    demo.stoch=demo.stoch,
                    reps=reps)

## Get the real values of r, N0, and K_step (can be a bit redundant when real.???.sd = 0)                    
expt$r.real <- rnorm(length(expt[,1]), mean=expt$r.real.mean, sd=expt$r.real.sd)
expt$N0.real <- rnorm(length(expt[,1]), mean=expt$N0.real.mean, sd=expt$N0.real.sd)
expt$K_step.real <- rnorm(length(expt[,1]), mean=expt$K_step.real.mean, sd=expt$K_step.real.sd)

## Get values of r and N0 to use for predictions
expt$r.pred <- rnorm(length(expt[,1]), mean=expt$r.real, sd=expt$r.pred.sd)
expt$N0.pred <- rnorm(length(expt[,1]), mean=expt$N0.real, sd=expt$N0.pred.sd)
expt$K_step.pred <- rnorm(length(expt[,1]), mean=expt$K_step.real, sd=expt$K_step.pred.sd)

## Check the experiment
str(expt)
##ggplot(expt, aes(x=N0.pred)) + geom_density() + facet_grid(N0.pred.sd ~ r.pred.sd)

## Choose how long to run the simulations for
its <- 50

## choose the starting carrying capacity
K_start <- 100

## Make the real dynanimcs
real.dyn <- lapply(1:length(expt[,1]),
                   function(x) iterate.ricker(r=expt[x, "r.real"],
                                              N=expt[x, "N0.real"],
                                              K_step=expt[x, "K_step.real"],
                                              its=its,
                                              demo.stoch=expt[x, "demo.stoch"]))

## check the dynamics
#par(mfrow=c(2,1))                                     
#plot(real.dyn[[1]], type="l")
#plot(real.dyn[[2]], type="l")

## Make the predicted dynanimcs
pred.dyn <- lapply(1:length(expt[,1]),
                   function(x) iterate.ricker(r=expt[x, "r.pred"],
                                              N=expt[x, "N0.pred"],
                                              K_step=expt[x, "K_step.pred"],
                                              its=its,
                                              demo.stoch=expt[x, "demo.stoch"]))

## check the dynamics
#par(mfrow=c(2,1))                                     
#plot(pred.dyn[[1]], type="l")
#plot(pred.dyn[[2]], type="l")

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
rr1 <- group_by(rr, N0.pred.sd, r.pred.sd, K_step.pred.sd, demo.stoch, its) %>%
  summarise(median.pred.skill=median(mov.cor),
            mean.pred.skill=mean(mov.cor),
            upper.CL=emp.CL(mov.cor, 55),
            lower.CL=emp.CL(mov.cor, 45))

## Get the prediction horizon for the averaged prediction skill
rr2 <- group_by(rr1, N0.pred.sd, r.pred.sd, K_step.pred.sd, demo.stoch) %>%
  summarise(pred.horizon=min(its[mean.pred.skill<pred.skill.threshold]-1),
            last.pred.skill=mean.pred.skill[min(its[mean.pred.skill<pred.skill.threshold]-1)])
rr3 <- rr2[rep(as.numeric(rownames(rr2)), each=2),]
rr3$last.pred.skill[seq(1,length(rr3$last.pred.skill), 2)] <- -0.2

## Save the data if not already done.
#save.image(file="fig2b.Rdata")



## Only run from here once a dataset is saved
rm(list=ls())
library(ggplot2)
library(repmis)

## load the already saved data from github (this can take some time depending on the internet connection)
source_data("https://github.com/opetchey/ecopredtools/blob/master/Petchey_etal_figures/data/fig2b.Rdata?raw=True")


## Plot the loss of prediction skill through time
rr1$nice.ds <- ifelse(rr1$demo.stoch, "With demographic stochasticity", "Without demographic stochasticity")  
rr3$nice.ds <- ifelse(rr3$demo.stoch, "With demographic stochasticity", "Without demographic stochasticity")  

rr1$N0.pred.CV <- rr1$N0.pred.sd / N0.real.mean
rr1$K_step.pred.CV <- rr1$K_step.pred.sd / K_step.real.mean
rr3$N0.pred.CV <- rr3$N0.pred.sd / N0.real.mean
rr3$K_step.pred.CV <- rr3$K_step.pred.sd / K_step.real.mean

g <- ggplot(data=rr1, aes(x=its, y=mean.pred.skill, col=as.factor(N0.pred.CV), linetype=as.factor(K_step.pred.CV))) +
  geom_line(size=0.7, alpha=0.7) + # transparaent colours to better see when lines lay on top of eachother
  #geom_point(size=1.5, alpha=0.5) + # make graph less busy
  labs(col="CV(N0)", linetype="CV(K_step)", x="Time (generations)", y="Forecast proficiency") +
  facet_grid(.~nice.ds) + 
  geom_hline(yintercept=pred.skill.threshold, col="purple", linetype=2, alpha=0.5)

g + theme_bw() + theme(legend.key = element_rect(colour = "white"))










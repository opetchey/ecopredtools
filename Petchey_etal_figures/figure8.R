## Code for box figure of Petchey et al, explaining the Lyapunov exponent
## Owen Petchey 1.12.14

## Plot three graphs
## 1. Two time series.
## 2. The absolute difference between these two time series.
## 3. Function relating forecast horizon with Lyapunov exponent.
                  
                  
## Use the logistic map
logistic <- function(r,n) r*n*(1-n)
get.iterate.logistic <- function(r,n0,tt){
  rez <- numeric(length=its)
  rez[1] <- n0
  ### iterate tt times
  for(i in 2:tt)
    rez[i] <- logistic(r,rez[i-1])
  rez
}
r <- 3.6 ## intrinsic growth rate
n0 <- 0.01 ## starting density of time series 1
e <- 1e-5 ## difference to starting density of time series 2
its <- 50 ## number of iterations
ts1 <- get.iterate.logistic(r,n0,its)                  
ts2 <- get.iterate.logistic(r,n0+e,its)

## Function relating forecast horizon to Lyapunov exponent
Tp <- function(le, Delta, delta) 1 / le * log(Delta / delta)      
le <- seq(0, 1, length=1000) ## the Lyapunov exponent
Delta <- 0.01 ## the accuracy threshold
delta <- 1e-10 ## the accuracy of the initial conditions

layout(matrix(1:3, 1, 3))

## Plot 1
plot(1:its, ts1, type="l", xlab="Time", ylab="Population size", log="y")
lines(1:its, ts2, type="l", col="blue")
mtext(3, line=1,text="(a)", adj=0, font=2, lwd=0.1)

## Plot 2
plot(1:its, log10(abs(ts1-ts2)), type="l", ylab="Log of absolute difference", xlab="Time")      
mtext(3, line=1,text="(b)", adj=0, font=2)
                  
## Plot 3                  
plot(le, log10(Tp(le, Delta, delta)), type="l",
 ylim=c(0.8, 4.5),
 ylab="Log10(Forecast horizon)", xlab="Lyapunov exponent")  
delta <- 1e-5                  
lines(le, log10(Tp(le, Delta, delta)), type="l")  
mtext(3, line=1,text="(c)", adj=0, font=2)

                  
                  
                  
                  
                  
                  
                  
                  
                  

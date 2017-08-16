library(deSolve)
?deSolve

example(deSolve)


# http://desolve.r-forge.r-project.org/user2014/examples/rigidODE.R.txt

## The RigidODE problem

library(deSolve)

rigidode <- function(t, y, parms) {
  dy1 <- -2  * y[2] * y[3]
  dy2 <- 1.25* y[1] * y[3]
  dy3 <- -0.5* y[1] * y[2]
  list(c(dy1, dy2, dy3))
}
yini  <- c(1, 0, 0.9)
times <- seq(from = 0, to = 20, by = 0.01)
out   <- ode (times = times, y = yini, func = rigidode, parms = NULL)
head (out, n = 3)

plot(out)

library(scatterplot3d)
par(mar = c(0, 0, 0, 0))
scatterplot3d(out[,-1], type = "l", lwd = 2, xlab = "", ylab = "", zlab = "", main = "rigidode")

# http://desolve.r-forge.r-project.org/user2014/examples/vanderpol/

vdpol <- function (t, y, mu) {
  list(c(
    y[2],
    mu * (1 - y[1]^2) * y[2] - y[1]
  ))
}

yini <- c(y1 = 2, y2 = 0)

stiff <- ode(y = yini, func = vdpol, times = 0:3000, parms = 1000)

nonstiff <- ode(y = yini, func = vdpol, times = seq(0, 30, by = 0.01), parms = 1)

head(stiff, n = 3)

plot(stiff, type = "l", which = "y1",lwd = 2, ylab = "y",main = "IVP ODE, stiff")

plot(nonstiff, type = "l", which = "y1",lwd = 2, ylab = "y",main = "IVP ODE, nonstiff")

system.time(
  stiff <- ode(y = yini, func = vdpol, times = 0:3000, parms = 1000,
               method = "bdf")
)

system.time(
  nonstiff <- ode(y = yini, func = vdpol,times = seq(0, 30, by = 0.01), parms = 1,
                  method = "bdf")
)

system.time(
  stiff <- ode(y = yini, func = vdpol, times = 0:3000, parms = 1000,
               method = "ode23")
)

system.time(
  nonstiff <- ode(y = yini, func = vdpol,times = seq(0, 30, by = 0.01), parms = 1,
                  method = "ode23")
)



# http://desolve.r-forge.r-project.org/user2014/examples/FME/

library(deSolve)
library(FME)
library(Cairo)

## A two compartment pharmacokinetic model
twocomp <- function (time, y, parms, ...) {
  with(as.list(c(parms, y)), {
    dCL <- kFL * CF - kLF * CL - ke * CL  # concentration in liver
    dCF <-    kLF * CL  - kFL * CF        # concentration in fat
    list(c(dCL, dCF))
  })
}
parms <- c(ke = 0.2,  kFL = 0.1,	kLF = 0.05)
times <- seq(0, 40, length=200)
y0      <-  c(CL = 1, CF = 0)
out <- ode(y0, times, twocomp, parms)


## Data
dat <- data.frame(
  time = seq(0, 28, 4),
  CL = c(1.31,  0.61,  0.49, 0.41,  0.20,  0.12,  0.16,  0.21),
  CF = c(1e-03, 0.041, 0.05, 0.039, 0.031, 0.025, 0.017, 0.012)
)

## simulate and plot the model and the data
out <- ode(y0, times, twocomp, parms)
plot(out, obs=dat, obspar=list(pch=16, col="red"))


## === fit the model =================================================

## define a cost function
cost <- function(p) {
  out  <-  ode(y0, times, twocomp, p)
  modCost(out, dat, weight = "none") # try "std" od "mean" for weight
}

## Note: initial parameters taken from above, may be adjusted manually
fit <- modFit(f = cost, p = parms)
summary(fit)

## Now plot original and fitted models and data
out1 <- ode(y0, times, twocomp, parms)
out2 <- ode(y0, times, twocomp, coef(fit))
plot(out1, out2, obs=dat, obspar=list(pch=16, col="red"))



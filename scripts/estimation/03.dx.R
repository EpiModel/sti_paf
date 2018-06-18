
# STI guidelines Diagnostics

suppressMessages(library("EpiModelHIV"))
rm(list = ls())

load("est/fit.rda")
data(est)

# Main model diagnostics

dx.m <- netdx(est[[1]], nsims = 10, nsteps = 500, dynamic = TRUE,
              set.control.ergm = control.simulate.ergm(MCMC.burnin = 1e6))

print(dx.m)
plot(dx.m, qnts.alpha = 0.9)

# Simulation to check cumulative main partners per year
nsims <- 5
means <- rep(NA, nsims)
for (JJ in 1:nsims) {
  dx.m <- netdx(est[[1]], nsims = 1, nsteps = 52*20, dynamic = TRUE, keep.tedgelist = TRUE,
                set.control.ergm = control.simulate.ergm(MCMC.burnin = 1e6), verbose = FALSE)
  ids <- c(dx.m$tedgelist[[1]]$head, dx.m$tedgelist[[1]]$tail)
  means[JJ] <- mean(tabulate(ids))
  cat(JJ)
}
summary(means)

# average per year
mean(means)/20

# expectation based on mean degree and duration
# Prev (cuml p number per week) = Incid (mean deg) * Duration (of Ps in weeks),
# so Incid = Prev/Duration
md.m <- 2*(est[[1]]$target.stats[1]/10000)
cuml.m.1y <- 52*(md.m/est[[1]]$coef.diss$duration)
cuml.m.1y

# 0.27 vs 0.26: good!

# Casual model diagnostics

dx.c <- netdx(est[[2]], nsims = 10, nsteps = 500, dynamic = TRUE,
              set.control.ergm = control.simulate.ergm(MCMC.burnin = 1e6))

print(dx.c)
plot(dx.c, qnts.alpha = 0.8)

nsims <- 5
means <- rep(NA, nsims)
for (JJ in 1:nsims) {
  dx.c <- netdx(est[[2]], nsims = 1, nsteps = 52*20, dynamic = TRUE, keep.tedgelist = TRUE,
                set.control.ergm = control.simulate.ergm(MCMC.burnin = 1e6), verbose = FALSE)
  ids <- c(dx.c$tedgelist[[1]]$head, dx.c$tedgelist[[1]]$tail)
  means[JJ] <- mean(tabulate(ids))
  cat(JJ)
}
summary(means)
mean(means)/20

md.c <- 2*(est[[2]]$target.stats[1]/10000)
cuml.c.1y <- 52*(md.c/est[[2]]$coef.diss$duration)
cuml.c.1y

# 0.91 vs 0.89: good!

# One-off model diagnostics

dx.i <- netdx(est[[3]], nsims = 10000, dynamic = FALSE)

print(dx.i)
# plot(dx.i)

sim <- simulate(est[[3]]$fit, nsim = 1000, monitor = ~nodefactor("riskg", base = 0), statsonly = TRUE)
stats <- sim[, 14:18]

mds.i <- c(0.0000,
           0.0010,
           0.0054,
           0.0102,
           0.0315)

52*colMeans(stats)/2000

est[[3]]
cuml.i.1y <- 365*mds.i
cuml.i.1y

52*colMeans(stats)/2000

tot.cuml.1y <- cuml.m.1y + cuml.c.1y + cuml.i.1y

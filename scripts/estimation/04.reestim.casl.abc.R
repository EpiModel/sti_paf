
library(methods)
suppressMessages(library(EpiModelHIV))
suppressMessages(library(doParallel))
suppressMessages(library(foreach))

## Environmental Arguments
simno <- as.numeric(Sys.getenv("PBS_ARRAYID"))

load("est/fit.rda")

# Pull casual network
est <- est[[2]]

# calibration function
f <- function(est) {

  ncoef <- length(est$coef.form) - 3
  for (jj in 1:ncoef) {
    est$coef.form[jj] <- est$coef.form[jj] + runif(1, -0.15, 0.15)
  }

  dx <- netdx(est, nsims = 1, nsteps = 500, verbose = FALSE,
              set.control.ergm = control.simulate.ergm(MCMC.burnin = 1e6))

  stats <- list()
  for (kk in 1:ncoef) {
    stats[[kk]] <- mean(tail(dx$stats[[1]][, kk], 100))
  }

  out <- as.numeric(c(est$coef.form[1:ncoef], stats))

  return(out)
}

# Run parallel on each node
registerDoParallel(16)
nsims <- 1000
sout <- foreach(s = 1:nsims) %dopar% {
  f(est)
}

sim <- as.data.frame(do.call("rbind", sout))
save(sim, file = paste0("data/sim.", simno, ".rda"))

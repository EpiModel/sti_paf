
# devtools::install_github("statnet/tergm")
library(EpiModelHIV)
sessionInfo()

# post-processing of data
fn <- list.files("data/", pattern = "sim", full.names = TRUE)

load(fn[1])
fsim <- sim
for (i in 2:length(fn)) {
  load(fn[i])
  fsim <- rbind(fsim, sim)
  cat("*")
}
dim(fsim)

load("est/fit.rda")
est <- est[[2]]
target.stats <- est$target.stats[1:4]

# rejection algorithm, weighted threshold
rejection <- function(sim,
                      target.stat = target.stats,
                      threshold = 0.05) {

  p <- sim[, 1:4]
  dat <- sim[, 5:8]

  diffs <- list()
  for (jj in 1:length(target.stats)) {
    diffs[[jj]] <- abs(dat[, jj] - target.stat[jj])
  }
  diffs <- as.data.frame(diffs)
  names(diffs) <- paste0("v", 1:length(target.stats))

  rdiff <- rowSums(diffs)
  cutoff <- quantile(rdiff, threshold)

  in.threshold <- which(rdiff <= cutoff)

  post <- sim[in.threshold, ]
  out <- list()
  out$param <- post[, 1:4]
  out$stats <- post[, 5:8]
  return(out)
}

post <- rejection(fsim, threshold = 0.001)
str(post)

# Accepted adjusted coefficients
selected.param <- colMeans(post$param)
selected.stats <- colMeans(post$stats)

cbind(target.stats, selected.stats)



# Test it
est$coef.form[1:length(target.stats)] <- selected.param

dx <- netdx(est, nsteps = 500, nsims = 20, ncores = 1, dynamic = TRUE,
            set.control.ergm = control.simulate.ergm(MCMC.burnin = 2e6))

print(dx)
plot(dx, qnts.alpha = 0.9)

# Write out to coefficients
load("est/fit.rda")
est[[2]]$coef.form[1:length(target.stats)] <- selected.param
save(est, file = "est/fit.rda")

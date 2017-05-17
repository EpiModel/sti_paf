
rm(list = ls())
suppressMessages(library("EpiModelHIV"))
devtools::load_all("~/Dropbox/Dev/EpiModelHIV/EpiModelHIV")

# Main Test Script ----------------------------------------------------

data(st)

param <- param_msm(nwstats = st,
                   prep.start = 26,
                   prep.coverage = 0.5)
init <- init_msm(nwstats = st)

control <- control_msm(simno = 1,
                       nsteps = 100,
                       nsims = 1,
                       ncores = 1)

data(est)
sim <- netsim(est, param, init, control)

df <- as.data.frame(sim)
names(df)
df$incid.gc
df$ir100.gc
df$prev.gc.dual
df$prev.ct.dual

df$times.rgc
df$times.ugc
df$ir100.sti.prep

plot(sim, y = "i.prev")
plot(sim, y = c("prev.rgcct", "prev.ugcct"),
     mean.col = 1:2, leg = TRUE)
plot(sim, y = c("ir100.gc", "ir100.ct"))


# Testing/Timing ------------------------------------------------------

control$bi.mods

undebug(sti_recov)
undebug(sti_tx)
undebug(prevalence_msm)
undebug(sti_trans)

dat <- initialize_msm(est, param, init, control, s = 1)

for (at in 2:100) {
  dat <- aging_msm(dat, at)       ## <1 ms
  dat <- deaths_msm(dat, at)      ## 4 ms
  dat <- births_msm(dat, at)      ## 6 ms
  dat <- test_msm(dat, at)        ## 2 ms
  dat <- tx_msm(dat, at)          ## 3 ms
  dat <- prep_msm(dat, at)        ## 2 ms
  dat <- progress_msm(dat, at)    ## 2 ms
  dat <- vl_msm(dat, at)          ## 3 ms
  dat <- simnet_msm(dat, at)      ## 53 ms
  dat <- disclose_msm(dat, at)    ## 1 ms
  dat <- acts_msm(dat, at)        ## 1 ms
  dat <- condoms_msm(dat, at)     ## 2 ms
  dat <- riskhist_msm(dat, at)    ## 4 ms
  dat <- position_msm(dat, at)    ## 1 ms
  dat <- trans_msm(dat, at)       ## 1 ms
  dat <- sti_trans(dat, at)       ## 4 ms
  dat <- sti_recov(dat, at)       ## 3 ms
  dat <- sti_tx(dat, at)          ## 2 ms
  dat <- prevalence_msm(dat, at)  ## 1 ms
  cat(at, ".", sep = "")
}

library(microbenchmark)

res <- microbenchmark(simnet_msm(dat, at = 2), times = 100)
summary(res, unit = "ms")


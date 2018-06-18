
library("methods")
suppressMessages(library("EpiModelHIV"))
suppressMessages(library("doParallel"))
suppressMessages(library("foreach"))
suppressMessages(library("EasyABC"))

f <- function(x) {

  set.seed(x[1])

  suppressMessages(library("EpiModelHIV"))

  data(st)

    param <- param_msm(nwstats = st,

                     ai.scale = 1.04,
                     ai.scale.pospos = 1.04,

                     tst.rect.sti.rr = 1,

                     rgc.tprob = x[2],
                     ugc.tprob = x[3],
                     rct.tprob = x[4],
                     uct.tprob = x[5],
                     syph.tprob = x[6],

                     # HIV acquisition
                     hiv.rgc.rr = 1.71,
                     hiv.ugc.rr = 1.24,
                     hiv.rct.rr = 1.71,
                     hiv.uct.rr = 1.24,
                     hiv.syph.rr = 1.61,

                     syph.incub.sympt.prob = 0,
                     syph.prim.sympt.prob = 0.82,
                     syph.seco.sympt.prob = 0.90,
                     syph.earlat.sympt.prob = 0,
                     syph.latelat.sympt.prob = 0,
                     syph.tert.sympt.prob = 1.0,

                     syph.prim.sympt.prob.tx = 0.85,
                     syph.seco.sympt.prob.tx = 0.85,
                     syph.earlat.sympt.prob.tx = 0.10,
                     syph.latelat.sympt.prob.tx = 0.10,
                     syph.tert.sympt.prob.tx = 1.0,

                     ept.coverage = 0.5,
                     stianntest.gc.hivneg.coverage = 0.44,
                     stianntest.ct.hivneg.coverage = 0.44,
                     stianntest.syph.hivneg.coverage = 0.44, # 0.45
                     stihighrisktest.gc.hivneg.coverage = 0.1,
                     stihighrisktest.ct.hivneg.coverage = 0.1,
                     stihighrisktest.syph.hivneg.coverage = 0.1,
                     stianntest.gc.hivpos.coverage = 0.61,
                     stianntest.ct.hivpos.coverage = 0.61,
                     stianntest.syph.hivpos.coverage = 0.65, #0.67
                     stihighrisktest.gc.hivpos.coverage = 0.1,
                     stihighrisktest.ct.hivpos.coverage = 0.1,
                     stihighrisktest.syph.hivpos.coverage = 0.1,

                     prep.start = 7000,
                     stitest.start = 5201,
                     ept.start = 5201,

                     #partlist.start = 1,
                     stitest.active.int = 364,
                     sti.highrisktest.int = 182,
                     ept.risk.int = 60)

    init <- init_msm(nwstats = st,
                     prev.ugc = 0.01,
                     prev.rgc = 0.01,
                     prev.uct = 0.01,
                     prev.rct = 0.01,
                     prev.syph.B = 0.01,
                     prev.syph.W = 0.01)

  control <- control_msm(simno = 1,
                         nsteps = 5200,
                         nsims = 1, ncores = 1,
                         verbose = FALSE)

  data(est)
  sim <- netsim(est, param, init, control)

  df <- tail(as.data.frame(sim), 52)

  gc.incid <- mean(tail(df$ir100.gc, 10))
  ct.incid <- mean(tail(df$ir100.ct, 10))
  hiv.prev <- mean(tail(df$i.prev, 10))
  syph.incid <- mean(tail(df$ir100.syph, 10))
  syph.prev <- mean(tail(df$prev.syph, 10))
  pssyph.prev <- mean(tail(df$prev.primsecosyph, 10))

  gcslope <- mean(df$ir100.gc[52] - df$ir100.gc[47])
  ctslope <- mean(df$ir100.ct[52] - df$ir100.ct[47])
  syphslope <- mean(df$ir100.syph[52] - df$ir100.syph[47])
  hivslope <- mean(df$ir100[52] - df$ir100[47])
  hivprevslope <- mean(df$i.prev[52] - df$i.prev[47])
  syphprevslope <- mean(df$prev.syph[52] - df$prev.syph[47])

  out <- c(gc.incid, ct.incid, hiv.prev, syph.incid,
           syph.prev, pssyph.prev,
           gcslope, ctslope, syphslope, hivslope,
           hivprevslope, syphprevslope)

  return(out)
}

priors <- list(c("unif", 0.511, 0.514), #rgc.tprob
               c("unif", 0.431, 0.433), #ugc.tprob
               c("unif", 0.279, 0.280), #rct.tprob
               c("unif", 0.216, 0.217), #uct.tprob
               c("unif", 0.120, 0.123)) #syph.tprob

# NG inc, CT inc, HIV prev, syph inc,
# syph prev, PS syph prev,
# GC inc slope, ct inc slope, syph inc slope, HIV inc slope,
# HIV prev slope, Syph prev slope
targets <- c(3.5, 5.6, 0.15, 2.6,
             0.012, 0.006,
             0, 0, 0, 0,
             0, 0)

# NHBS NG/CT testing (Hoots 2011 and 2014 self-report data):
# 2014: NG: 46.2% HIV-MSM, 64.1%  HIV+ MSM
# 2014: CT: 45.8% HIV-MSM, 62.8%  HIV+ MSM
# NHBS syphilis testing (2014 self-report data Qian An): 45% HIV- MSM, 68% HIV+ MSM
# Flagg STD 2015 (MMP - 2008-2010 data among MSM): 54% syphilis, 20% NG, 20% CT
# Mattson CID 2017 (MMP- 2009-2013 data among MSM):
# Syphilis (2009-2013): 54, 57, 58, 60, 66
# NG/CT (2009-2013): 18, 22, 26, 32, 39

( nsim <- as.numeric(Sys.getenv("NSIM")) )
( pacc <- as.numeric(Sys.getenv("PACC")) )

a <- ABC_sequential(method = "Lenormand",
                    model = f,
                    prior = priors,
                    nb_simul = nsim,
                    summary_stat_target = targets,
                    p_acc_min = pacc,
                    progress_bar = TRUE,
                    n_cluster = 16,
                    use_seed = TRUE,
                    verbose = FALSE)

fn <- paste0("data/smc4.", pacc*100, "pct.", nsim, "sim.rda")
save(a, file = fn)

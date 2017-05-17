
library("methods")
suppressMessages(library("EpiModelHIV"))
suppressMessages(library("doParallel"))
suppressMessages(library("foreach"))
suppressMessages(library("EasyABC"))

f <- function(x) {

  set.seed(x[1])

  suppressMessages(library("EpiModelHIV"))

  param <- param_msm(nwstats = st,

                     rgc.tprob = x[2],
                     ugc.tprob = x[3],
                     rct.tprob = x[4],
                     uct.tprob = x[5],

                     rgc.sympt.prob = x[6],
                     ugc.sympt.prob = x[7],
                     rct.sympt.prob = x[8],
                     uct.sympt.prob = x[9],

                     rgc.dur.asympt = x[10],
                     ugc.dur.asympt = x[11],
                     rct.dur.asympt = x[12],
                     uct.dur.asympt = x[13],

                     gc.prob.cease = 0,
                     ct.prob.cease = 0,

                     hiv.rgc.rr = x[14],
                     hiv.ugc.rr = x[15],
                     hiv.rct.rr = x[14],
                     hiv.uct.rr = x[15])

  init <- init_msm(nwstats = st, prev.ugc = 0.01, prev.rgc = 0.01,
                   prev.uct = 0.01, prev.rct = 0.01)

  control <- control_msm(simno = 1,
                         nsteps = 1820,
                         nsims = 1, ncores = 1,
                         verbose = FALSE)

  data(est)
  sim <- netsim(est, param, init, control)

  df <- tail(as.data.frame(sim), 52)

  gc.incid <- mean(df$ir100.gc)
  ct.incid <- mean(df$ir100.ct)
  hiv.incid <- mean(df$ir100)
  hiv.prev <- mean(df$i.prev)

  out <- c(gc.incid, ct.incid, hiv.incid, hiv.prev)

  return(out)
}

priors <- list(c("unif", 0.30, 0.60),
               c("unif", 0.20, 0.50),
               c("unif", 0.30, 0.60),
               c("unif", 0.20, 0.50),
               c("unif", 0.01, 0.15),
               c("unif", 0.60, 0.95),
               c("unif", 0.01, 0.15),
               c("unif", 0.60, 0.95),
               c("unif", 26, 52),
               c("unif", 26, 52),
               c("unif", 39, 65),
               c("unif", 39, 65),
               c("unif", 2, 3),
               c("unif", 1, 2))

targets <- c(4.2, 6.6, 3.8, 0.26)

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

fn <- paste0("data/smc3.", pacc*100, "pct.", nsim, "sim.rda")
save(a, file = fn)

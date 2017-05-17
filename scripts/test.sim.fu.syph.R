
rm(list=ls())
library("methods")
suppressMessages(library("EpiModelHIV"))
# devtools::load_all("~/Dropbox/Dev/EpiModelHIV/EpiModelHIV")

load("est/nwstats.rda")
param <- param_msm(nwstats = st,
                   ai.scale = 1.12, # 1.11889726, # was 1.13
                   
                   rsyph.tprob = 0.057,
                   usyph.tprob = 0.047,
                   
                   hiv.rsyph.rr = 2.7, 
                   hiv.usyph.rr = 1.7,
                   syph.rhiv.rr = 4.00,
                   syph.uhiv.rr = 3.00,
                   
                   syph.earlat.rr = 0.5, #2/3, 0
                   incu.syph.int = 27,
                   prim.syph.int = 60,
                   seco.syph.int = 120,
                   earlat.syph.int = 365 - 27 - 60 - 120,
                   latelat.syph.int = 9 * 52 * 7,
                   latelatelat.syph.int = 20 * 52 * 7,
                   tert.syph.int = 20 * 52 * 7,
                   # immune.syph.int = 5 * 52 * 7,
                   syph.tert.prog.prob = 0.15 / (52 * 7 * 20),
                   
                   rgc.tprob = 0.42, 
                   ugc.tprob = 0.315, 
                   rct.tprob = 0.21, 
                   uct.tprob = 0.15,
                   
                   hiv.rgc.rr = 2.7, #2.780673,
                   hiv.ugc.rr = 1.7, #1.732363,
                   hiv.rct.rr = 2.7, #2.780673,
                   hiv.uct.rr = 1.7, #1.732363,
                   
                   # adjust prim and seco from 0.1385 each
                   stage.syph.B.prob = c(0.00, 0.20, 0.077, 0.277, 0.22, 0.22, 0.006),
                   stage.syph.W.prob = c(0.00, 0.20, 0.077, 0.277, 0.22, 0.22, 0.006),
                   
                   syph.prim.sympt.prob.tx = 0.35, # Tuite PLoS One 2014, Bissessor AIDS 2010, Kourbatova STD 2008 use 0.45
                   syph.prim.asympt.prob.tx = 0.00,
                   syph.seco.sympt.prob.tx = 0.60, # Tuite PLoS One 2014, Bissessor AIDS 2010, Kourbatova STD 2008
                   syph.seco.asympt.prob.tx = 0.00,
                   syph.earlat.prob.tx = 0.15, # Tuite PLoS One 2014, Bissessor AIDS 2010, Kourbatova STD 2008
                   syph.latelat.prob.tx = 0.10,
                   syph.tert.sympt.prob.tx = 0.90,
                   syph.tert.asympt.prob.tx = 0.00,
                   
                   hivdx.syph.sympt.tx.rr = 2.00,
                   
                   prep.coverage = 0.4,
                   ept.coverage = 0,
                   
                   prep.start = 2601,
                   stitest.start = 5000,
                   
                   stitest.active.int = 364,
                   sti.highrisktest.int = 182) # adjustable for 3 or 6 months

init <- init_msm(st)

control <- control_msm(simno = 1,
                       start = 2601, nsteps = 2700,
                       nsims = 1, ncores = 1,
                       initialize.FUN = reinit_msm)

load("est/stimod.burnin.rda")

sim2 <- netsim(sim, param, init, control)


dat <- reinit_msm(sim, param, init, control, s = 1)

for (at in 2601:2650) {
    dat <- aging_msm(dat, at)
    dat <- deaths_msm(dat, at)
    dat <- births_msm(dat, at)
    dat <- test_msm(dat, at)
    dat <- tx_msm(dat, at)
    dat <- prep_msm(dat, at)
    dat <- progress_msm(dat, at)
    dat <- vl_msm(dat, at)
    dat <- simnet_msm(dat, at)
    dat <- disclose_msm(dat, at)
    dat <- acts_msm(dat, at)
    dat <- condoms_msm(dat, at)
    dat <- riskhist_msm(dat, at)
    dat <- position_msm(dat, at)
    dat <- trans_msm(dat, at)
    dat <- sti_trans(dat, at)
    dat <- sti_recov(dat, at)
    dat <- sti_tx(dat, at)
    dat <- prevalence_msm(dat, at)
    cat(at, ".", sep = "")
}
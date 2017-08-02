rm(list=ls())
library("methods")
suppressMessages(library("EpiModelHIV"))
suppressMessages(library("EpiModelHPC"))
# devtools::load_all("~/Dropbox/Dev/EpiModelHIV/EpiModelHIV")


gctrans <- 1.0
cttrans <- 1.0
syphtrans <- 1.0
gccttrans <- 1.0
gcsyphtrans <- 1.0
ctsyphtrans <- 1.0
allstitrans <- 1.0



load("est/nwstats.rda")
load("est/stimod.burnin.rda")
param <- param_msm(nwstats = st,
                   ai.scale = 1.03,
                   
                   syph.earlat.rr = 0.5,
                   incu.syph.int = 27,
                   prim.syph.int = 60,
                   seco.syph.int = 120,
                   earlat.syph.int = 365 - 27 - 60 - 120,
                   latelat.syph.int = 9 * 52 * 7,
                   latelatelat.syph.int = 20 * 52 * 7,
                   tert.syph.int = 20 * 52 * 7,
                   syph.tert.prog.prob = 0.00015625599,
                   
                   # STI acquisition
                   rgc.tprob = 0.447,
                   ugc.tprob = 0.337,
                   rct.tprob = 0.2025,
                   uct.tprob = 0.1825,
                   rsyph.tprob = 0.1526957,
                   usyph.tprob = 0.1326838,
                   
                   # HIV acquisition
                   hiv.rgc.rr = 1.80292790,
                   hiv.ugc.rr = 1.1989083,
                   hiv.rct.rr = 1.80292790,
                   hiv.uct.rr = 1.1989083,
                   hiv.rsyph.rr = 1.80292790,
                   hiv.usyph.rr = 1.1989083,
                   
                   # HIV transmission
                   hiv.trans.gc.rr = gctrans,
                   hiv.trans.ct.rr = cttrans,
                   hiv.trans.syph.rr = syphtrans,
                   hiv.trans.gc.ct.rr = gccttrans,
                   hiv.trans.gc.syph.rr = gcsyphtrans,
                   hiv.trans.ct.syph.rr = ctsyphtrans,
                   hiv.trans.allsti.rr = allstitrans,
                   
                   syph.prim.sympt.prob.tx = 0.60,
                   syph.seco.sympt.prob.tx = 0.688235,
                   syph.earlat.sympt.prob.tx = 0.10,
                   syph.latelat.sympt.prob.tx = 0.10,
                   syph.tert.sympt.prob.tx = 1.0,
                   
                   syph.prim.asympt.prob.tx = 1.0,
                   syph.seco.asympt.prob.tx = 1.0,
                   syph.earlat.asympt.prob.tx = 1.0,
                   syph.latelat.asympt.prob.tx = 1.0,
                   syph.tert.asympt.prob.tx = 1.0,
                   gc.asympt.prob.tx = 1.0,
                   ct.asympt.prob.tx = 1.0,
                   
                   hivdx.syph.sympt.tx.rr = 1.5,
                   
                   partnercut = 1,
                   stianntest.coverage = 0.1,
                   stihighrisktest.coverage = 0.0,
                   prep.coverage = 0,
                   ept.coverage = 0,
                   
                   prep.start = 7000,
                   stitest.start = 5201,
                   ept.start = 7000,
                   
                   stitest.active.int = 364,
                   sti.highrisktest.int = 182) # adjustable for 3 or 6 months


init <- init_msm(st)

control <- control_msm(simno = 1,
                       start = 5201,
                       nsteps = 5202,
                       nsims = 1,
                       ncores = 1,
                       initialize.FUN = reinit_msm,
                       verbose = FALSE)

# sim2 <- netsim(sim, param, init, control)

# colMeans(sim2$epi$prop.CT.asympt.tx, na.rm = TRUE)
# colMeans(sim2$epi$prop.GC.asympt.tx, na.rm = TRUE)
# colMeans(sim2$epi$prop.rGC.tx, na.rm = TRUE)
# colMeans(sim2$epi$prop.rCT.tx, na.rm = TRUE)
# 
# 
# # Testing/Timing ------------------------------------------------------
# 
# control$bi.mods
# debug(test_sti_msm)
# debug(sti_tx)
# # debug(sti_recov)
# # debug(prevalence_msm)
# 
# # load("est/stimod.burnin.rda")

at <- 5200
dat <- reinit_msm(sim, param, init, control, s = 1)

at <- at + 1

dat <- aging_msm(dat, at)
dat <- deaths_msm(dat, at)
dat <- births_msm(dat, at)
dat <- hiv_test_msm(dat, at)
dat <- sti_test_msm(dat, at)
dat <- hiv_tx_msm(dat, at)
dat <- prep_msm(dat, at)
dat <- sti_ept_msm(dat, at)
dat <- hiv_progress_msm(dat, at)
dat <- syph_progress_msm(dat, at)
dat <- hiv_vl_msm(dat, at)
dat <- simnet_msm(dat, at)
dat <- hiv_disclose_msm(dat, at)
dat <- part_msm(dat, at)
dat <- acts_msm(dat, at)
dat <- condoms_msm(dat, at)
dat <- riskhist_prep_msm(dat, at)
dat <- riskhist_stitest_msm(dat, at)
dat <- position_msm(dat, at)

sum_GC <- rep(NA,length(1:1000))
sum_CT <- rep(NA,length(1:1000))
sum_syph <- rep(NA,length(1:1000))
sum_urethral <- rep(NA,length(1:1000))
sum_rectal <- rep(NA,length(1:1000))
cell1_gc <- rep(NA,length(1:1000))
cell2_gc <- rep(NA,length(1:1000))
cell3_gc <- rep(NA,length(1:1000))
cell4_gc <- rep(NA,length(1:1000))
cell1_ct <- rep(NA,length(1:1000))
cell2_ct <- rep(NA,length(1:1000))
cell3_ct <- rep(NA,length(1:1000))
cell4_ct <- rep(NA,length(1:1000))
cell1_syph <- rep(NA,length(1:1000))
cell2_syph <- rep(NA,length(1:1000))
cell3_syph <- rep(NA,length(1:1000))
cell4_syph <- rep(NA,length(1:1000))
cell1_sti <- rep(NA,length(1:1000))
cell2_sti <- rep(NA,length(1:1000))
cell3_sti <- rep(NA,length(1:1000))
cell4_sti <- rep(NA,length(1:1000))
# Summary Output
incid <- rep(NA,length(1:1000))
ir100 <- rep(NA,length(1:1000))
trans.main <- rep(NA,length(1:1000))
trans.pers <- rep(NA,length(1:1000))
trans.inst <- rep(NA,length(1:1000))

df.direct <- data.frame(sum_GC,sum_CT,sum_syph,sum_urethral,sum_rectal,
                        cell1_gc,cell2_gc,cell3_gc,cell4_gc,
                        cell1_ct,cell2_ct,cell3_ct,cell4_ct,
                        cell1_syph,cell2_syph,cell3_syph,cell4_syph,
                        cell1_sti,cell2_sti,cell3_sti,cell4_sti,
                        incid,ir100,trans.main,trans.pers,trans.inst)

for (sim in 1:1000) {
  hiv_trans_msm(dat,at)
}

# 
# for (at in 2601:2650) {
#   dat <- aging_msm(dat, at)
#   dat <- deaths_msm(dat, at)
#   dat <- births_msm(dat, at)
#   dat <- hiv_test_msm(dat, at)
#   dat <- sti_test_msm(dat, at)
#   dat <- hiv_tx_msm(dat, at)
#   dat <- prep_msm(dat, at)
#   dat <- sti_ept_msm(dat, at)
#   dat <- hiv_progress_msm(dat, at)
#   dat <- syph_progress_msm(dat, at)
#   dat <- hiv_vl_msm(dat, at)
#   dat <- simnet_msm(dat, at)
#   dat <- hiv_disclose_msm(dat, at)
#   dat <- part_msm(dat, at)
#   dat <- acts_msm(dat, at)
#   dat <- condoms_msm(dat, at)
#   dat <- riskhist_msm(dat, at)
#   dat <- position_msm(dat, at)
#   dat <- hiv_trans_msm(dat, at)
#   dat <- sti_trans_msm(dat, at)
#   dat <- sti_recov_msm(dat, at)
#   dat <- sti_tx_msm(dat, at)
#   dat <- prevalence_msm(dat, at)
#   verbose_msm(dat, type = "progress", s = 1, at)
#   cat(at, ".", sep = "")
# }
# 
# 
# 
# undebug(prep_msm)
# debug(sti_tx)

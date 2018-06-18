
## Packages
library("methods")
suppressMessages(library("EpiModelHIV"))
suppressMessages(library("EpiModelHPC"))

## Environmental Arguments
simno <- Sys.getenv("SIMNO")
jobno <- Sys.getenv("PBS_ARRAYID")
njobs <- as.numeric(Sys.getenv("NJOBS"))
fsimno <- paste(simno, jobno, sep = ".")

## Parameters
load("est/nwstats.rda")
#load("est/abc.syph.parms.rda")
#for (i in seq_along(mean.p)) {
#    assign(names(mean.p)[i], unname(mean.p[i]))
#}

param <- param_msm(nwstats = st,

                   ai.scale = 1.04,
                   ai.scale.pospos = 1.04,

                   # Correlation
                   sti.stitx.correlation = "false",
                   sti.hivdx.correlation = "false",
                   sti.correlation.time = 12,

                   # STI acquisition
                   rgc.tprob = 0.4773,
                   ugc.tprob = 0.3819,
                   rct.tprob = 0.2564,
                   uct.tprob = 0.2091,
                   syph.tprob = 0.2533,

                   # HIV acquisition
                   hiv.rgc.rr = 1.75,
                   hiv.ugc.rr = 1.26,
                   hiv.rct.rr = 1.75,
                   hiv.uct.rr = 1.26,
                   hiv.syph.rr = 1.63,

                   syph.incub.sympt.prob = 0,
                   syph.prim.sympt.prob = 0.70,
                   syph.seco.sympt.prob = 0.85,
                   syph.earlat.sympt.prob = 0,
                   syph.latelat.sympt.prob = 0,
                   syph.tert.sympt.prob = 1.0,

                   syph.prim.sympt.prob.tx = 0.80,
                   syph.seco.sympt.prob.tx = 0.80,
                   syph.earlat.sympt.prob.tx = 0.10,
                   syph.latelat.sympt.prob.tx = 0.10,
                   syph.tert.sympt.prob.tx = 1.0,

                   # HIV transmission
                   hiv.trans.gc.rr = 1.0,
                   hiv.trans.ct.rr = 1.0,
                   hiv.trans.syph.rr = 1.0,

                   ept.coverage = 0.0,
                   stianntest.gc.hivneg.coverage = 0.44,
                   stianntest.ct.hivneg.coverage = 0.44,
                   stianntest.syph.hivneg.coverage = 0.45,
                   stihighrisktest.gc.hivneg.coverage = 0.0,
                   stihighrisktest.ct.hivneg.coverage = 0.0,
                   stihighrisktest.syph.hivneg.coverage = 0.0,
                   stianntest.gc.hivpos.coverage = 0.61,
                   stianntest.ct.hivpos.coverage = 0.61,
                   stianntest.syph.hivpos.coverage = 0.67,
                   stihighrisktest.gc.hivpos.coverage = 0.0,
                   stihighrisktest.ct.hivpos.coverage = 0.0,
                   stihighrisktest.syph.hivpos.coverage = 0.0,
                   # 2014 (Hoots): NG: 46.2% HIV-MSM, 64.1%  HIV+ MSM
                   # 2014 (Hoots): CT: 45.8% HIV-MSM, 62.8%  HIV+ MSM
                   # NHBS syphilis testing (2014 self-report data Qian An):
                   # 45% HIV- MSM, 68% HIV+ MSM

                   # Balanced:
                   # HIV-negative MSM:
                   # NG: 46.2% (Hoots unpublished NHBS)
                   # CT: 45.8 % (Hoots unpublished NHBS)
                   # Syphilis: 45% (An 2017 - NHBS)

                   #   HIV-positive MSM:
                   # NG: 43% (Mattson 2017 - MMP), 47.2% (Patel 2017 – MMP)
                   # CT: 43% (Mattson 2017 - MMP), 47.2% (Patel 2017 – MMP)
                   # Syphilis: 69% (Mattson 2017 - MMP)

                   # Ratio approach
                   # NG among HIV-negative MSM: (47.2) x (46.2/64.1) = 34.0%
                   # CT among HIV-negative MSM:(47.2) x (45.8 /62.8) = 34.4%
                   # Syph among HIV-negative MSM:(45) x (45 /68) = 29.8%


                   prep.start = 7000,
                   stitest.start = 5201,
                   ept.start = 7000,

                   stitest.active.int = 364,
                   sti.highrisktest.int = 182,
                   ept.risk.int = 60)

init <- init_msm(nwstats = st,
                 prev.ugc = 0.003,
                 prev.rgc = 0.003,
                 prev.uct = 0.003,
                 prev.rct = 0.003,
                 prev.syph.B = 0.005,
                 prev.syph.W = 0.005)

control <- control_msm(simno = fsimno,
                       nsteps = 5200,
                       nsims = 16, ncores = 16,
                       verbose = FALSE)

## Simulation
netsim_hpc("est/fit.rda", param, init, control,
            save.min = FALSE, save.max = TRUE)

process_simfiles(simno = simno, min.n = njobs, compress = TRUE, outdir = "data/")

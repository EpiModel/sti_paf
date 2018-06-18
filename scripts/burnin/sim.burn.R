
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
                   ai.scale.pospos = 2.08,

                   tst.rect.sti.rr = 1,

                   # Correlation
                   sti.correlation.time = 12,

                   # STI acquisition
                   rgc.tprob = 0.5161, #0.513,
                   ugc.tprob = 0.4362, # 0.432
                   rct.tprob = 0.2813, #0.2797,
                   uct.tprob = 0.2195, # 0.2165,
                   syph.tprob = 0.1206,

                   # HIV acquisition
                   hiv.rgc.rr = 1.97, #1.75,
                   hiv.ugc.rr = 1.48, #1.27,
                   hiv.rct.rr = 1.97, #1.75,
                   hiv.uct.rr = 1.48, #1.27,
                   hiv.syph.rr = 1.64,

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

                   ept.coverage = 0.0,
                   stianntest.gc.hivneg.coverage = 0.44,
                   stianntest.ct.hivneg.coverage = 0.44,
                   stianntest.syph.hivneg.coverage = 0, #0.44, # 0.45
                   stihighrisktest.gc.hivneg.coverage = 0.0,
                   stihighrisktest.ct.hivneg.coverage = 0.0,
                   stihighrisktest.syph.hivneg.coverage = 0.0,
                   stianntest.gc.hivpos.coverage = 0.61,
                   stianntest.ct.hivpos.coverage = 0.61,
                   stianntest.syph.hivpos.coverage = 0, #0.65, #0.67
                   stihighrisktest.gc.hivpos.coverage = 0.0,
                   stihighrisktest.ct.hivpos.coverage = 0.0,
                   stihighrisktest.syph.hivpos.coverage = 0.0,

                   prep.start = 7000,
                   stitest.start = 5201,
                   ept.start = 5201,

                   #partlist.start = 1,
                   stitest.active.int = 364,
                   sti.highrisktest.int = 182,
                   ept.risk.int = 60)

init <- init_msm(nwstats = st,
                 prev.ugc = 0.0015,
                 prev.rgc = 0.0015,
                 prev.uct = 0.0015,
                 prev.rct = 0.0015, # 0.03
                 prev.syph.B = 0.01, # 0.03
                 prev.syph.W = 0.01) # 0.03

control <- control_msm(simno = fsimno,
                       nsteps = 5200,
                       nsims = 16, ncores = 16,
                       verbose = FALSE)

## Simulation
netsim_hpc("est/fit.rda", param, init, control,
            save.min = FALSE, save.max = TRUE)

process_simfiles(simno = simno, min.n = njobs, compress = TRUE, outdir = "data/")

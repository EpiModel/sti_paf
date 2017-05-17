
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
                   ai.scale = 1.11,
                   
                   rsyph.tprob = 0.07, #0.04668348,
                   usyph.tprob = 0.05, #0.03598792,
                   
                   hiv.rsyph.rr = 2.80, #2.98876572, 
                   hiv.usyph.rr = 2.10, #1.7456618,
                   syph.rhiv.rr = 1, #6.54189295,
                   syph.uhiv.rr = 1, #5.09641658,
                   
                   syph.earlat.rr = 0.5,
                   incu.syph.int = 27,
                   prim.syph.int = 60,
                   seco.syph.int = 120,
                   earlat.syph.int = 365 - 27 - 60 - 120,
                   latelat.syph.int = 9 * 52 * 7,
                   latelatelat.syph.int = 20 * 52 * 7,
                   tert.syph.int = 20 * 52 * 7,
                   syph.tert.prog.prob = 0.00015625599,
                   
                   rgc.tprob = 0.422, #0.4133300,
                   ugc.tprob = 0.310, #0.31404720,
                   rct.tprob = 0.195,
                   uct.tprob = 0.1655,
                   
                   hiv.rgc.rr = 2.55, #2.35,
                   hiv.ugc.rr = 1.855, #1.35,
                   hiv.rct.rr = 2.55, #2.35,
                   hiv.uct.rr = 1.855, #1.35,
                   
                   syph.prim.sympt.prob.tx = 0.35, # Tuite PLoS One 2014, Bissessor AIDS 2010, Kourbatova STD 2008 use 0.45
                   syph.seco.sympt.prob.tx = 0.60, # Tuite PLoS One 2014, Bissessor AIDS 2010, Kourbatova STD 2008
                   syph.earlat.sympt.prob.tx = 0.15, # Tuite PLoS One 2014, Bissessor AIDS 2010, Kourbatova STD 2008
                   syph.latelat.sympt.prob.tx = 0.10,
                   syph.tert.sympt.prob.tx = 0.90,
                   
                   syph.prim.asympt.prob.tx = 1,
                   syph.seco.asympt.prob.tx = 1,
                   syph.earlat.asympt.prob.tx = 1, 
                   syph.latelat.asympt.prob.tx = 1,
                   syph.tert.asympt.prob.tx = 1,
                   
                   hivdx.syph.sympt.tx.rr = 1.45,
                   
                   prep.coverage = 0.0,
                   ept.coverage = 0.0,
                   stianntest.coverage = 0.5,
                   stihighrisktest.coverage = 0.8,

                   prep.start = 5000,
                   stitest.start = 2601,
                   ept.start = 5000,
                   
                   stitest.elig.model = "sti",
                   
                   stitest.active.int = 364,
                   sti.highrisktest.int = 182) # adjustable for 3 or 6 months

init <- init_msm(nwstats = st, 
                 prev.B = 0.10, 
                 prev.W = 0.10,
                 prev.ugc = 0.015,
                 prev.rgc = 0.015,
                 prev.uct = 0.015,
                 prev.rct = 0.015,
                 prev.syph.B = 0.020,
                 prev.syph.W = 0.020,
                 
                 # adjust prim and seco from 0.1385 each
                 # Incubating, primary, secondary, early latent, late latent, late late latent, tertiary
                 stage.syph.B.prob = c(0.00, 0.20, 0.077, 0.277, 0.22, 0.22, 0.006),
                 stage.syph.W.prob = c(0.00, 0.20, 0.077, 0.277, 0.22, 0.22, 0.006))

control <- control_msm(simno = fsimno,
                       nsteps = 2600,
                       nsims = 16, ncores = 16,
                       verbose = FALSE)

## Simulation
netsim_hpc("est/fit.rda", param, init, control,
            save.min = TRUE, save.max = TRUE)

process_simfiles(simno = simno, min.n = njobs, compress = TRUE, outdir = "data/")

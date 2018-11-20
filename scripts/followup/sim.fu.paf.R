
## Packages
library("methods")
suppressMessages(library("EpiModelHIV"))
suppressMessages(library("EpiModelHPC"))
library("EpiModel")

## Environmental Arguments
simno <- as.numeric(Sys.getenv("SIMNO"))
jobno <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
njobs <- as.numeric(Sys.getenv("NJOBS"))
fsimno <- paste(simno, jobno, sep = ".")
#aipospos <- as.numeric(Sys.getenv("AIPOSPOS"))
#syphacq <- as.numeric(Sys.getenv("SYPHACQ"))
rectacq <- as.numeric(Sys.getenv("RECTACQ"))
urethacq <- as.numeric(Sys.getenv("URETHACQ"))
#syphtrans <- as.numeric(Sys.getenv("SYPHTRANS"))
gctrans <- as.numeric(Sys.getenv("GCTRANS"))
cttrans <- as.numeric(Sys.getenv("CTTRANS"))
gccttrans <- as.numeric(Sys.getenv("GCCTTRANS"))
# ctsyphtrans <- as.numeric(Sys.getenv("CTSYPHTRANS"))
# gcsyphtrans <- as.numeric(Sys.getenv("GCSYPHTRANS"))
allstitrans <- as.numeric(Sys.getenv("ALLSTITRANS"))

cat("Array number is ", jobno)
cat("\n fsimno is ", fsimno)

## Parameters
load("est/nwstats.rda")

param <- param_msm(nwstats = st,

                   ai.scale = 1.04,
                   #ai.scale.pospos = 2.08,

                   # Correlation
                   sti.correlation.time = 12,

                   # STI acquisition
                   rgc.tprob = 0.5161, #0.513,
                   ugc.tprob = 0.4362, # 0.432
                   rct.tprob = 0.2813, #0.2797,
                   uct.tprob = 0.2195, # 0.2165,
                   syph.tprob = 0, #0.1206,

                   # HIV acquisition
                   hiv.rgc.rr = rectacq,
                   hiv.ugc.rr = urethacq,
                   hiv.rct.rr = rectacq,
                   hiv.uct.rr = urethacq,
                   #hiv.syph.rr = syphacq,

                   # HIV transmission
                   hiv.trans.gc.rr = gctrans,
                   hiv.trans.ct.rr = cttrans,
                   #hiv.trans.syph.rr = syphtrans,
                   hiv.trans.gc.ct.rr = gccttrans,
                   #hiv.trans.gc.syph.rr = gcsyphtrans,
                   #hiv.trans.ct.syph.rr = ctsyphtrans,
                   hiv.trans.allsti.rr = allstitrans,

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

                   syph.prim.asympt.prob.tx = 1.0,
                   syph.seco.asympt.prob.tx = 1.0,
                   syph.earlat.asympt.prob.tx = 1.0,
                   syph.latelat.asympt.prob.tx = 1.0,
                   syph.tert.asympt.prob.tx = 1.0,
                   gc.asympt.prob.tx = 1.0,
                   ct.asympt.prob.tx = 1.0,

                   partnercut = 1,
                   ept.coverage = 0.0,
                   stianntest.gc.hivneg.coverage = 0.44,
                   stianntest.ct.hivneg.coverage = 0.44,
                   stianntest.syph.hivneg.coverage = 0.0, #0.44, # 0.45
                   stihighrisktest.gc.hivneg.coverage = 0.0,
                   stihighrisktest.ct.hivneg.coverage = 0.0,
                   stihighrisktest.syph.hivneg.coverage = 0.0,
                   stianntest.gc.hivpos.coverage = 0.61,
                   stianntest.ct.hivpos.coverage = 0.61,
                   stianntest.syph.hivpos.coverage = 0.0, #0.65, #0.67
                   stihighrisktest.gc.hivpos.coverage = 0.0,
                   stihighrisktest.ct.hivpos.coverage = 0.0,
                   stihighrisktest.syph.hivpos.coverage = 0.0,
                   prep.coverage = 0,

                   prep.start = 7000,
                   stitest.start = 7000,
                   ept.start = 7000,

                   stitest.active.int = 364,
                   sti.highrisktest.int = 182) # adjustable for 3 or 6 months


init <- init_msm(st)

control <- control_msm(simno = fsimno,
                       start = 5201,
                       nsteps = 5720,
                       nsims = 16,
                       ncores = 16,
                       initialize.FUN = reinit_msm,
                       verbose = FALSE)

# Simulation
netsim_hpc("est/stimod.burnin.rda", param, init, control,
           compress = TRUE, verbose = FALSE)

process_simfiles(simno = simno, min.n = njobs,
                 outdir = "data/", compress = TRUE, delete.sub = TRUE,
                 truncate.at = 5200,
                 vars =
                   c("num", "ir100", "incid", "ir100.gc", "incid.gc", "incid.gcct",
                     "ir100.ct", "incid.ct", "ir100.syph", "incid.syph", "incid.sti",
                     "ir100.rct", "ir100.uct", "ir100.rgc", "ir100.ugc",
                     "ir100.sti", "ir100.gcct",
                     "incid.gc.hivneg", "incid.gc.hivpos",
                     "incid.ct.hivneg", "incid.ct.hivpos",
                     "incid.syph.hivneg", "incid.syph.hivpos",
                     "ir100.gc.hivneg", "ir100.gc.hivpos",
                     "ir100.ct.hivneg", "ir100.ct.hivpos",
                     "ir100.syph.hivneg", "ir100.syph.hivpos",
                     "prop.edges.negneg", "prop.edges.negpos", "prop.edges.pospos",
                     "num.acts.negneg", "num.acts.negpos", "num.acts.pospos",
                     "prop.uai.negneg", "prop.uai.negpos", "prop.uai.pospos",
                     "prop.acts.negneg", "prop.acts.negpos", "prop.acts.pospos",
                     "prop.main.edges.negneg", "prop.main.edges.negpos",
                     "prop.main.edges.pospos", "prop.cas.edges.negneg",
                     "prop.cas.edges.negpos", "prop.cas.edges.pospos",
                     "prop.inst.edges.negneg", "prop.inst.edges.negpos", "prop.inst.edges.pospos",
                     "hivtests.nprep", "hivtests.pos", "hivtests.prep",
                     'test.gc.12mo', 'test.gc.12mo.hivpos', 'test.gc.12mo.hivneg',
                     'test.ct.12mo', 'test.ct.12mo.hivpos', 'test.ct.12mo.hivneg',
                     "prev.gc", "prev.rgc", "prev.ugc",
                     "prev.ct", "prev.rct", "prev.uct", "prev.sti",
                     "prev.rgc.hivneg", "prev.ugc.hivneg",
                     "prev.rct.hivneg", "prev.uct.hivneg",
                     "prev.rgc.hivpos","prev.ugc.hivpos",
                     "prev.rct.hivpos", "prev.uct.hivpos",
                     "prev.hivposmultsti", "prev.hivnegmultsti",
                     "txGC", "txCT",
                     "txGC_asympt", "txCT_asympt", "txsyph_asympt", "txSTI", "txSTI_asympt",
                     "sum_GC", "sum_CT", "sum_syph", "sum_urethral", "sum_rectal",
                     "cell1_rectureth", "cell2_rectureth", "cell3_rectureth", "cell4_rectureth",
                     "cell1_newinf", "cell2_newinf", "cell3_newinf", "cell4_newinf",
                     "cell1_gc", "cell2_gc", "cell3_gc", "cell4_gc",
                     "cell1_ct", "cell2_ct", "cell3_ct", "cell4_ct",
                     "cell1_sti", "cell2_sti", "cell3_sti", "cell4_sti",
                     "stiactiveind.prop", "stiactiveind",
                     "recentpartners", "recentpartners.prop",
                     "rGC_hiv", "uGC_hiv", "rCT_hiv", "uCT_hiv", "rGC_hiv_acq",
                     "uGC_hiv_acq", "rCT_hiv_acq", "uCT_hiv_acq", "rGC_hiv_trans",
                     "uGC_hiv_trans", "rCT_hiv_trans", "uCT_hiv_trans",
                     "rGC_hiv_trans_events_perperson",
                     "uGC_hiv_trans_events_perperson",
                     "rCT_hiv_trans_events_perperson",
                     "uCT_hiv_trans_events_perperson"))

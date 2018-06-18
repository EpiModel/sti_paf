
## Packages
library("methods")
suppressMessages(library("EpiModelHIV"))
suppressMessages(library("EpiModelHPC"))

## Environmental Arguments
simno <- as.numeric(Sys.getenv("SIMNO"))
jobno <- as.numeric(Sys.getenv("PBS_ARRAYID"))
njobs <- as.numeric(Sys.getenv("NJOBS"))
fsimno <- paste(simno, jobno, sep = ".")
syphacq <- as.numeric(Sys.getenv("SYPHACQ"))
rectacq <- as.numeric(Sys.getenv("RECTACQ"))
urethacq <- as.numeric(Sys.getenv("URETHACQ"))
syphtrans <- as.numeric(Sys.getenv("SYPHTRANS"))
gctrans <- as.numeric(Sys.getenv("GCTRANS"))
cttrans <- as.numeric(Sys.getenv("CTTRANS"))
gccttrans <- as.numeric(Sys.getenv("GCCTTRANS"))
ctsyphtrans <- as.numeric(Sys.getenv("CTSYPHTRANS"))
gcsyphtrans <- as.numeric(Sys.getenv("GCSYPHTRANS"))
allstitrans <- as.numeric(Sys.getenv("ALLSTITRANS"))
anngcnegcov <- as.numeric(Sys.getenv("ANNGCNEGCOV"))
anngcposcov <- as.numeric(Sys.getenv("ANNGCPOSCOV"))
annctnegcov <- as.numeric(Sys.getenv("ANNCTNEGCOV"))
annctposcov <- as.numeric(Sys.getenv("ANNCTPOSCOV"))
annsyphnegcov <- as.numeric(Sys.getenv("ANNSYPHNEGCOV"))
annsyphposcov <- as.numeric(Sys.getenv("ANNSYPHPOSCOV"))
hrgcnegcov <- as.numeric(Sys.getenv("HRGCNEGCOV"))
hrgcposcov <- as.numeric(Sys.getenv("HRGCPOSCOV"))
hrctnegcov <- as.numeric(Sys.getenv("HRCTNEGCOV"))
hrctposcov <- as.numeric(Sys.getenv("HRCTPOSCOV"))
hrsyphnegcov <- as.numeric(Sys.getenv("HRSYPHNEGCOV"))
hrsyphposcov <- as.numeric(Sys.getenv("HRSYPHPOSCOV"))

## Parameters
load("est/nwstats.rda")

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
                   syph.tprob = 0.2553,

                   # HIV acquisition
                   hiv.rgc.rr = rectacq,
                   hiv.ugc.rr = urethacq,
                   hiv.rct.rr = rectacq,
                   hiv.uct.rr = urethacq,
                   hiv.syph.rr = syphacq,

                   # HIV transmission
                   hiv.trans.gc.rr = gctrans,
                   hiv.trans.ct.rr = cttrans,
                   hiv.trans.syph.rr = syphtrans,
                   hiv.trans.gc.ct.rr = gccttrans,
                   hiv.trans.gc.syph.rr = gcsyphtrans,
                   hiv.trans.ct.syph.rr = ctsyphtrans,
                   hiv.trans.allsti.rr = allstitrans,

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

                   syph.prim.asympt.prob.tx = 1.0,
                   syph.seco.asympt.prob.tx = 1.0,
                   syph.earlat.asympt.prob.tx = 1.0,
                   syph.latelat.asympt.prob.tx = 1.0,
                   syph.tert.asympt.prob.tx = 1.0,
                   gc.asympt.prob.tx = 1.0,
                   ct.asympt.prob.tx = 1.0,

                   partnercut = 1,
                   ept.coverage = 0.0,
                   stianntest.gc.hivneg.coverage = anngcnegcov,
                   stianntest.ct.hivneg.coverage = annctnegcov,
                   stianntest.syph.hivneg.coverage = annsyphnegcov,
                   stihighrisktest.gc.hivneg.coverage = hrgcnegcov,
                   stihighrisktest.ct.hivneg.coverage = hrctnegcov,
                   stihighrisktest.syph.hivneg.coverage = hrsyphnegcov,
                   stianntest.gc.hivpos.coverage = anngcposcov,
                   stianntest.ct.hivpos.coverage = annctposcov,
                   stianntest.syph.hivpos.coverage = annsyphposcov,
                   stihighrisktest.gc.hivpos.coverage = hrgcposcov,
                   stihighrisktest.ct.hivpos.coverage = hrctposcov,
                   stihighrisktest.syph.hivpos.coverage = hrsyphposcov,
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

## Simulation
netsim_hpc("est/stimod.burnin.rda", param, init, control,
           compress = TRUE, verbose = FALSE)

process_simfiles(simno = simno, min.n = njobs,
                 outdir = "data/", compress = TRUE, delete.sub = TRUE,
                 truncate.at = 5200,
                 vars =
                   c("num", "ir100", "incid", "ir100.gc", "incid.gc", "incid.gcct",
                     "ir100.ct", "incid.ct", "ir100.syph", "incid.syph", "incid.sti",
                     "ir100.rct", "ir100.uct", "ir100.rgc", "ir100.ugc",
                     "ir100.sti", "ir100.sti.prep", "ir100.gcct",
                     "ir100.sti.tttraj1", "ir100.sti.tttraj2",
                     "ir100.gc.tttraj1", "ir100.gc.tttraj2",
                     "ir100.ct.tttraj1", "ir100.ct.tttraj2",
                     "ir100.syph.tttraj1", "ir100.syph.tttraj2",
                     "ir100.gcct.tttraj1", "ir100.gcct.tttraj2",
                     "incid.sti.tttraj1", "incid.sti.tttraj2",
                     "incid.gc.tttraj1", "incid.gc.tttraj2",
                     "incid.ct.tttraj1", "incid.ct.tttraj2",
                     "incid.syph.tttraj1", "incid.syph.tttraj2",
                     "incid.gcct.tttraj1", "incid.gcct.tttraj2",
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
                     "GCasympttests", "uGCasympttests", "rGCasympttests",
                     "CTasympttests", "uCTasympttests", "rCTasympttests",
                     "syphasympttests", "stiasympttests",
                     "GCasympttests.pos", "uGCasympttests.pos", "rGCasympttests.pos",
                     "CTasympttests.pos", "uCTasympttests.pos", "rCTasympttests.pos",
                     "syphasympttests.pos", "syphearlyasympttests.pos",
                     "syphlateasympttests.pos", "stiasympttests.pos",
                     # "GCasympttests.hivneg", "uGCasympttests.hivneg", "rGCasympttests.hivneg",
                     # "CTasympttests.hivneg", "uCTasympttests.hivneg", "rCTasympttests.hivneg",
                     # "syphasympttests.hivneg", "stiasympttests.hivneg",
                     #"GCasympttests.pos.hivneg", "uGCasympttests.pos.hivneg", "rGCasympttests.pos.hivneg",
                     #"CTasympttests.pos.hivneg", "uCTasympttests.pos.hivneg", "rCTasympttests.pos.hivneg",
                     #"syphasympttests.pos.hivneg", "syphearlyasympttests.pos.hivneg",
                     #"syphlateasympttests.pos.hivneg", "stiasympttests.pos.hivneg",
                     # "GCasympttests.hivpos", "uGCasympttests.hivpos", "rGCasympttests.hivpos",
                     # "CTasympttests.hivpos", "uCTasympttests.hivpos", "rCTasympttests.hivpos",
                     #"syphasympttests.hivpos", "stiasympttests.hivpos",
                     #"GCasympttests.pos.hivpos", "uGCasympttests.pos.hivpos", "rGCasympttests.pos.hivpos",
                     #"CTasympttests.pos.hivpos", "uCTasympttests.pos.hivpos", "rCTasympttests.pos.hivpos",
                     #"syphasympttests.pos.hivpos", "syphearlyasympttests.pos.hivpos",
                     #"syphlateasympttests.pos.hivpos", "stiasympttests.pos.hivpos",
                     "tt.traj.syph1", "tt.traj.gc1", "tt.traj.ct1",
                     "tt.traj.syph2", "tt.traj.gc2", "tt.traj.ct2",
                     "tt.traj.sti1", "tt.traj.sti2",
                     "GCsympttests", "uGCsympttests", "rGCsympttests",
                     "CTsympttests", "uCTsympttests", "rCTsympttests",
                     "syphsympttests", "stisympttests",
                     "rCTasympttests.tttraj1", "rCTasympttests.tttraj2",
                     "uCTasympttests.tttraj1", "uCTasympttests.tttraj2",
                     "CTasympttests.tttraj1", "CTasympttests.tttraj2",
                     "rGCasympttests.tttraj1", "rGCasympttests.tttraj2",
                     "uGCasympttests.tttraj1", "uGCasympttests.tttraj2",
                     "GCasympttests.tttraj1", "GCasympttests.tttraj2",
                     "syphasympttests.tttraj1", "syphasympttests.tttraj2",
                     "stiasympttests.tttraj1", "stiasympttests.tttraj2",
                     "rCTsympttests.tttraj1", "rCTsympttests.tttraj2",
                     "uCTsympttests.tttraj1", "uCTsympttests.tttraj2",
                     "CTsympttests.tttraj1", "CTsympttests.tttraj2",
                     "rGCsympttests.tttraj1", "rGCsympttests.tttraj2",
                     "uGCsympttests.tttraj1", "uGCsympttests.tttraj2",
                     "GCsympttests.tttraj1", "GCsympttests.tttraj2",
                     "syphsympttests.tttraj1", "syphsympttests.tttraj2",
                     "stisympttests.tttraj1", "stisympttests.tttraj2",
                     "txGC.tttraj1", "txGC_asympt.tttraj1",
                     "txGC.tttraj2", "txGC_asympt.tttraj2",
                     "txCT.tttraj1", "txCT_asympt.tttraj1",
                     "txCT.tttraj2", "txCT_asympt.tttraj2",
                     "txsyph.tttraj1", "txsyph_asympt.tttraj1",
                     "txsyph.tttraj2", "txsyph_asympt.tttraj2",
                     "txearlysyph.tttraj1", "txlatesyph.tttraj1",
                     "txearlysyph.tttraj2", "txlatesyph.tttraj2",
                     "txSTI.tttraj1", "txSTI.tttraj2",
                     "txSTI_asympt.tttraj1", "txSTI_asympt.tttraj2",
                     "hivtests.nprep", "hivtests.pos", "hivtests.prep",
                     'test.gc.12mo', 'test.gc.12mo.hivpos', 'test.gc.12mo.hivneg',
                     'test.ct.12mo', 'test.ct.12mo.hivpos', 'test.ct.12mo.hivneg',
                     'test.syph.12mo', 'test.syph.12mo.hivpos', 'test.syph.12mo.hivneg',
                     "i.prev", "prev.primsecosyph", "prev.syph", "prev.earlysyph", "prev.latesyph",
                     "prev.gc", "prev.rgc", "prev.ugc",
                     "prev.ct", "prev.rct", "prev.uct", "prev.sti",
                     "prev.dxhiv.dxipssyph",
                     "prev.primsecosyph.tttraj1", "prev.primsecosyph.tttraj2",
                     "prev.syph.tttraj1", "prev.syph.tttraj2",
                     "prev.gcct.tttraj1", "prev.gcct.tttraj2",
                     "prev.gc.tttraj1", "prev.gc.tttraj2",
                     "prev.ct.tttraj1", "prev.ct.tttraj2",
                     "prev.sti.tttraj1", "prev.sti.tttraj2",
                     "prev.rgc.hivneg.only", "prev.ugc.hivneg.only",
                     "prev.rct.hivneg.only", "prev.uct.hivneg.only",
                     "prev.primsecosyph.hivneg.only",
                     "prev.rgc.hivpos.only","prev.ugc.hivpos.only",
                     "prev.rct.hivpos.only", "prev.uct.hivpos.only",
                     "prev.primsecosyph.hivpos.only",
                     "prev.hivposmultsti", "prev.hivnegmultsti",
                     "txearlysyph", "txlatesyph", "txsyph", "txGC", "txCT",
                     "txGC_asympt", "txCT_asympt", "txsyph_asympt", "txSTI", "txSTI_asympt",
                     "tx.gc.prop", "tx.ct.prop", "tx.gcct.prop", "tx.syph.prop",
                     "gc.infect.dur", "ct.infect.dur", "gcct.infect.dur", "syph.infect.dur",
                     "sum_GC", "sum_CT", "sum_syph", "sum_urethral", "sum_rectal",
                     "cell1_gc", "cell2_gc", "cell3_gc", "cell4_gc",
                     "cell1_ct", "cell2_ct", "cell3_ct", "cell4_ct",
                     "cell1_syph", "cell2_syph", "cell3_syph", "cell4_syph",
                     "cell1_sti", "cell2_sti", "cell3_sti", "cell4_sti",
                     "deathage", "stiactiveind.prop", "stiactiveind",
                     "recentpartners", "recentpartners.prop",
                     "num.newearlydiagsyph",
                     "num.newlatediagsyph", "early.late.syphratio",
                     "early.late.diagsyphratio",
                     "num.asympt.tx", "num.asympt.cases", "num.rect.cases", "num.rect.tx",
                     "time.hivneg",
                     "eptCov", "eptpartelig", "eptpartprovided", "eptpartuptake",
                     "eptTx", "propindexeptElig", "eptprop_provided", "eptprop_tx",
                     "eptuninfectedprovided","eptuninfecteduptake","eptgcinfectsti",
                     "eptctinfectsti","eptgcinfecthiv", "eptctinfecthiv",
                     "stage.time.ar.ndx","stage.time.ar.dx", "stage.time.af.ndx",
                     "stage.time.af.dx", "stage.time.early.chronic.ndx",
                     "stage.time.early.chronic.dx.yrone",
                     "stage.time.early.chronic.dx.yrstwotolate",
                     "stage.time.early.chronic.art", "stage.time.late.chronic.ndx",
                     "stage.time.late.chronic.dx", "stage.time.late.chronic.art",
                     "stage.time.aids.ndx", "stage.time.aids.dx","stage.time.aids.art"))

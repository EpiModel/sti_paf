## Process STI Testing Guidelines Data

rm(list = ls())
library("EpiModelHIV")
library("EpiModelHPC")
library("dplyr")
source("analysis/fx.R")

# ( fn <- list.files("data/followup", full.names = TRUE) )
# # fn <- list.files(pattern = "data/followup/n[3-4][0-9][0-9][0-9].rda")
# for (i in fn) {
#     load(i)
#     sim <- truncate_sim(sim, at = 2600)
#     save(sim, file = i, compress = TRUE)
#     cat("*")
# }

# truncate and limit sim 3000/4000 files on Hyak
# cd /gscratch/csde/kweiss2/sti/data
# module load r_3.2.4
# R
fn <- list.files(pattern = "n[3-4][0-9][0-9][0-9].rda")
for (i in fn) {
  load(i)
  sim <- truncate_sim(sim, at = 2600)
  vars.needed <- c("num", "ir100", "incid", "ir100.gc", "incid.gc",
                   "ir100.ct", "incid.ct", "ir100.syph", "incid.syph",
                   "ir100.rct", "ir100.uct", "ir100.rgc", "ir100.ugc",
                   "ir100.sti", "ir100.sti.prep", 
                   "GCasympttests", "uGCasympttests", "rGCasympttests", 
                   "CTasympttests", "uCTasympttests", "rCTasympttests", 
                   "syphasympttests", "stiasympttests",
                   "GCasympttests.pos", "uGCasympttests.pos", "rGCasympttests.pos", 
                   "CTasympttests.pos", "uCTasympttests.pos", "rCTasympttests.pos", 
                   "syphasympttests.pos", "stiasympttests.pos",
                   "GCsympttests", "uGCsympttests", "rGCsympttests", 
                   "CTsympttests", "uCTsympttests", "rCTsympttests", 
                   "syphsympttests", "stisympttests",
                   "hivtests", "hivtests.pos", "hivtests.prep",
                   "totalGCasympttests", "totaluGCasympttests", "totalrGCasympttests", 
                   "totalCTasympttests", "totaluCTasympttests", "totalrCTasympttests", 
                   "totalsyphasympttests", "totalstiasympttests",
                   "totalGCasympttests.pos", "totaluGCasympttests.pos", "totalrGCasympttests.pos", 
                   "totalCTasympttests.pos", "totaluCTasympttests.pos", "totalrCTasympttests.pos", 
                   "totalsyphasympttests.pos", "totalstiasympttests.pos",
                   "totalGCsympttests", "totaluGCsympttests", "totalrGCsympttests", 
                   "totalCTsympttests", "totaluCTsympttests", "totalrCTsympttests", 
                   "totalsyphsympttests", "totalstisympttests",
                   "stiactiveind", "recentpartners","recentSTI","newpartner", "concurrpart", "partnersti", "uai.nmain","uai.any",
                   "i.prev", "prev.primsecosyph", "prev.syph",
                   "prev.gc", "prev.rgc", "prev.ugc",
                   "prev.ct", "prev.rct", "prev.uct",
                   "txearlysyph", "txlatesyph", "txsyph", "txGC", "txCT",
                   "hiv_sum", "sti_hiv_sum", "sti_u_hiv_sum", "sti_r_hiv_sum", "sti_syph_hiv_sum", 
                   "sti_paf", "sti_u_paf", "sti_r_paf", "sti_syph_paf", 
                   "sti_u_paf", "sti_u_sympt_paf", "sti_u_asympt_paf","sti_r_paf", 
                   "sti_r_sympt_paf", "sti_r_asympt_paf", "sti_syph_paf", "sti_syph_sympt_paf", "sti_syph_asympt_paf",
                   "totalhivtests", "totalhivtests.pos", "totalhivtests.prep", "time.hivneg",
                   "stage.time.ar.ndx", "stage.time.af.ndx", "stage.time.chronic.ndx", "stage.time.aids.ndx",
                   "stage.time.ar.dx", "stage.time.af.dx", "stage.time.chronic.dx", "stage.time.aids.dx",
                   "stage.time.ar.art", "stage.time.af.art", "stage.time.chronic.art", "stage.time.aids.art")
  
  i.vars <- which(names(sim$epi) %in% vars.needed)
  sim$epi <- sim$epi[i.vars]
  out.fn <- paste0("followup/", i)
  save(sim, file = out.fn, compress = "gzip")
  #file.remove(i)
  cat(i, "\n")
}


### 1 by 1 processing on Hyak - cd /gscratch/csde/kweiss2/sti/data
rm(list = ls())
library("EpiModelHIV")
library("EpiModelHPC")
library("dplyr")
load("sim.n2000.rda")
sim <- truncate_sim(sim, at = 2600)
vars.needed <- c("num", "ir100", "incid", "ir100.gc", "incid.gc",
                 "ir100.ct", "incid.ct", "ir100.syph", "incid.syph",
                 "ir100.rct", "ir100.uct", "ir100.rgc", "ir100.ugc",
                 "ir100.sti", "ir100.sti.prep", 
                 "GCasympttests", "uGCasympttests", "rGCasympttests", 
                 "CTasympttests", "uCTasympttests", "rCTasympttests", 
                 "syphasympttests", "stiasympttests",
                 "GCasympttests.pos", "uGCasympttests.pos", "rGCasympttests.pos", 
                 "CTasympttests.pos", "uCTasympttests.pos", "rCTasympttests.pos", 
                 "syphasympttests.pos", "stiasympttests.pos",
                 "GCsympttests", "uGCsympttests", "rGCsympttests", 
                 "CTsympttests", "uCTsympttests", "rCTsympttests", 
                 "syphsympttests", "stisympttests",
                 "hivtests", "hivtests.pos", "hivtests.prep",
                 "totalGCasympttests", "totaluGCasympttests", "totalrGCasympttests", 
                 "totalCTasympttests", "totaluCTasympttests", "totalrCTasympttests", 
                 "totalsyphasympttests", "totalstiasympttests",
                 "totalGCasympttests.pos", "totaluGCasympttests.pos", "totalrGCasympttests.pos", 
                 "totalCTasympttests.pos", "totaluCTasympttests.pos", "totalrCTasympttests.pos", 
                 "totalsyphasympttests.pos", "totalstiasympttests.pos",
                 "totalGCsympttests", "totaluGCsympttests", "totalrGCsympttests", 
                 "totalCTsympttests", "totaluCTsympttests", "totalrCTsympttests", 
                 "totalsyphsympttests", "totalstisympttests",
                 "stiactiveind", "recentpartners","recentSTI","newpartner", "concurrpart", "partnersti", "uai.nmain","uai.any",
                 "i.prev", "prev.primsecosyph", "prev.syph",
                 "prev.gc", "prev.rgc", "prev.ugc",
                 "prev.ct", "prev.rct", "prev.uct",
                 "txearlysyph", "txlatesyph", "txsyph", "txGC", "txCT",
                 "hiv_sum", "sti_hiv_sum", "sti_u_hiv_sum", "sti_r_hiv_sum", "sti_syph_hiv_sum", 
                 "sti_paf", "sti_u_paf", "sti_r_paf", "sti_syph_paf", 
                 "sti_u_paf", "sti_u_sympt_paf", "sti_u_asympt_paf","sti_r_paf", 
                 "sti_r_sympt_paf", "sti_r_asympt_paf", "sti_syph_paf", "sti_syph_sympt_paf", "sti_syph_asympt_paf",
                 "totalhivtests", "totalhivtests.pos", "totalhivtests.prep", "time.hivneg",
                 "stage.time.ar.ndx", "stage.time.af.ndx", "stage.time.chronic.ndx", "stage.time.aids.ndx",
                 "stage.time.ar.dx", "stage.time.af.dx", "stage.time.chronic.dx", "stage.time.aids.dx",
                 "stage.time.ar.art", "stage.time.af.art", "stage.time.chronic.art", "stage.time.aids.art")

i.vars <- which(names(sim$epi) %in% vars.needed)
sim$epi <- sim$epi[i.vars]
save(sim, file = "followup/sim.n2000.rda", compress = "gzip")


## Locally merge files
sim <- merge_simfiles(3118, indir = "data/", ftype = "min")
sim <- truncate_sim(sim, at = 2600)
vars.needed <- c("num", "ir100", "incid", "ir100.gc", "incid.gc",
                 "ir100.ct", "incid.ct", "ir100.syph", "incid.syph",
                 "ir100.rct", "ir100.uct", "ir100.rgc", "ir100.ugc",
                 "ir100.sti", "ir100.sti.prep", 
                 "GCasympttests", "uGCasympttests", "rGCasympttests", 
                 "CTasympttests", "uCTasympttests", "rCTasympttests", 
                 "syphasympttests", "stiasympttests",
                 "GCasympttests.pos", "uGCasympttests.pos", "rGCasympttests.pos", 
                 "CTasympttests.pos", "uCTasympttests.pos", "rCTasympttests.pos", 
                 "syphasympttests.pos", "stiasympttests.pos",
                 "GCsympttests", "uGCsympttests", "rGCsympttests", 
                 "CTsympttests", "uCTsympttests", "rCTsympttests", 
                 "syphsympttests", "stisympttests",
                 "hivtests", "hivtests.pos", "hivtests.prep",
                 "totalGCasympttests", "totaluGCasympttests", "totalrGCasympttests", 
                 "totalCTasympttests", "totaluCTasympttests", "totalrCTasympttests", 
                 "totalsyphasympttests", "totalstiasympttests",
                 "totalGCasympttests.pos", "totaluGCasympttests.pos", "totalrGCasympttests.pos", 
                 "totalCTasympttests.pos", "totaluCTasympttests.pos", "totalrCTasympttests.pos", 
                 "totalsyphasympttests.pos", "totalstiasympttests.pos",
                 "totalGCsympttests", "totaluGCsympttests", "totalrGCsympttests", 
                 "totalCTsympttests", "totaluCTsympttests", "totalrCTsympttests", 
                 "totalsyphsympttests", "totalstisympttests",
                 "stiactiveind", "recentpartners","recentSTI","newpartner", "concurrpart", "partnersti", "uai.nmain","uai.any",
                 "i.prev", "prev.primsecosyph", "prev.syph",
                 "prev.gc", "prev.rgc", "prev.ugc",
                 "prev.ct", "prev.rct", "prev.uct",
                 "txearlysyph", "txlatesyph", "txsyph", "txGC", "txCT",
                 "hiv_sum", "sti_hiv_sum", "sti_u_hiv_sum", "sti_r_hiv_sum", "sti_syph_hiv_sum", 
                 "sti_paf", "sti_u_paf", "sti_r_paf", "sti_syph_paf", 
                 "sti_u_paf", "sti_u_sympt_paf", "sti_u_asympt_paf","sti_r_paf", 
                 "sti_r_sympt_paf", "sti_r_asympt_paf", "sti_syph_paf", "sti_syph_sympt_paf", "sti_syph_asympt_paf",
                 "totalhivtests", "totalhivtests.pos", "totalhivtests.prep", "time.hivneg",
                 "stage.time.ar.ndx", "stage.time.af.ndx", "stage.time.chronic.ndx", "stage.time.aids.ndx",
                 "stage.time.ar.dx", "stage.time.af.dx", "stage.time.chronic.dx", "stage.time.aids.dx",
                 "stage.time.ar.art", "stage.time.af.art", "stage.time.chronic.art", "stage.time.aids.art")

i.vars <- which(names(sim$epi) %in% vars.needed)
sim$epi <- sim$epi[i.vars]
save(sim, file = "data/followup/sim.3118.rda", compress = "gzip")
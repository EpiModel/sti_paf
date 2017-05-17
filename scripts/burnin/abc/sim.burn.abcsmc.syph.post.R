
## abc model analysis

library("EpiModelHIV")
library("EasyABC")

# system("scp scripts/burnin/abc/*.abcsmc.[Rs]* hyak:/gscratch/csde/sjenness/sti2")
# system("scp scripts/burnin/abc/*.abcsmc2.[Rs]* hyak:/gscratch/csde/sjenness/sti2")
# system("scp scripts/burnin/abc/*.abcsmc3.[Rs]* hyak:/gscratch/csde/sjenness/sti2")
# system("scp scripts/burnin/abc/*.abcsmc4.[Rs]* hyak:/gscratch/csde/sjenness/stia")
# 
# system("scp hyak:/gscratch/csde/sjenness/stia/data/*.rda scripts/burnin/abc/")
# 
# 
# ## averaged ATL/demo fits
# load("scripts/burnin/abc/smc.5pct.100sim.rda")
# load("scripts/burnin/abc/smc2.5pct.100sim.rda")
# load("scripts/burnin/abc/smc3.10pct.100sim.rda")
# load("scripts/burnin/abc/smc3.5pct.100sim.rda")
# load("scripts/burnin/abc/smc3.1pct.100sim.rda")


p <- as.data.frame(a$param)
s <- as.data.frame(a$stats)
w <- a$weights

# names(p) <- c("rgc.tprob", "ugc.tprob", "rct.tprob", "uct.tprob",
#               "rgc.sympt.prob", "ugc.sympt.prob", "rct.sympt.prob", "uct.sympt.prob",
#               "rgc.dur.asympt", "ugc.dur.asympt", "rct.dur.asympt", "uct.dur.asympt",
#               "hiv.rect.rr", "hiv.ureth.rr")

# names(p) <- c("ai.scale", "rsyph.tprob", "usyph.tprob", "hiv.rsyph.rr", "hiv.usyph.rr", "syph.hiv.rr", "rgc.tprob",
#               "ugc.tprob", "rct.tprob", "uct.tprob", "hiv.rct.rr", "hiv.uct.rr", "hiv.rgc.rr", "hiv.ugc.rr", 
#               "syph.prim.sympt.prob.tx", "syph.seco.sympt.prob.tx", "syph.earlat.prob.tx",
#               "syph.latelat.prob.tx")

# names(p) <- c("rsyph.tprob", "usyph.tprob", "rectalsti.rr", "urethralsti.rr", "syph.rhiv.rr", "syph.uhiv.rr", 
#               "rgc.tprob", "ugc.tprob", "rct.tprob","uct.tprob", "syph.prim.sympt.prob.tx", 
#               "syph.seco.sympt.prob.tx", "syph.earlat.prob.tx", "syph.latelat.prob.tx")

names(p) <- c("rsyph.tprob", "usyph.tprob", "hiv.rsyph.rr", "hiv.usyph.rr", "rgc.tprob","ugc.tprob", "rct.tprob","uct.tprob", "hiv.rsti.rr", "hiv.usti.rr") 
              # 
              # "ugc.tprob", 
              #"rct.tprob",
              # "uct.tprob"#, "hivdx.syph.sympt.tx.rr")

# names(s) <- c("gc.incid", "ct.incid", "hiv.incid", #syph.incid,
#              "hiv.prev", "prev.primsecosyph.hivpos", "prev.primsecosyph.hivneg", "prev.primsecosyph",
#              "prev.hiv.primsecosyphpos", "gcslope", "ctslope", "hivslope", "hivprevslope", "syphprevslope")#, "prev.earlysyph", "prev.latesyph")

# names(s) <- c("gc.incid", "ct.incid", "hiv.incid", "hiv.prev", "prev.primsecosyph",
#               "gcslope", "ctslope", "hivslope", "hivprevslope", "syphprevslope")#, "prev.earlysyph", "prev.latesyph")

names(s) <- c("gc.incid", "ct.incid", "hiv.incid", "syph.incid","hiv.prev")#,
              #"gcslope", "ctslope", "syphslope", "hivslope")
              #, "hivprevslope", "syphprevslope")

View(s)
View(p)

( mean.s <- apply(s, 2, function(x) sum(x * w)) )
( mean.p <- apply(p, 2, function(x) sum(x * w)) )

# hist(s$prev.primsecosyph.hivpos / s$prev.primsecosyph.hivneg)

#tar.syph <- c(4.2, 6.6, 3.8, 0.26, 0.046, 0, 0, 0, 0, 0) 


# tar.syph <- c(4.2, 6.6, 3.8, #0.9, 
#               0.26, 0.103, 0.026, 0.046, 0.498, 0.554, 0.446, 0, 0, 0, 0, 0) #, 0.1385, 0.1385, 0.2770, 0.2000, 0.2000, 0.0460)

data.frame(mean.s, tar.syph)

mean.p

par(mar = c(3,3,1,1), mgp = c(2,1,0), mfrow = c(3,2))
for (i in 1:ncol(s)) {
  hist(s[, i], col = "bisque2", border = "white", main = names(s)[i])
  abline(v = tar[i], lwd = 2, col = "red")
}

par(mar = c(3,3,1,1), mgp = c(2,1,0), mfrow = c(4,4))
for (i in 1:ncol(p)) {
  hist(p[, i], col = "bisque2", border = "white", main = names(p)[i])
}

for (i in seq_along(mean.p)) {
    assign(names(mean.p)[i], unname(mean.p[i]))
}

save(mean.p, file = "scripts/burnin/abc/abc.avg.parms.syph.100sim.1pct.rda")
save(mean.p, file = "est/abc.syph.parms.rda")
#save(mean.p, file = "est/meta.parms.rda")

for (i in seq_along(mean.p)) {
  assign(names(mean.p)[i], unname(mean.p[i]))
}

mean.p <- c(0.0245, 2.40, 2.30)
names(mean.p) <- c("syph.tprob", "hiv.syph.rr", "syph.hiv.rr")
save(mean.p, file = "est/abc.syph.parms.rda")

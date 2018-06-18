
## abc model analysis

library("EpiModelHIV")
library("EasyABC")

p <- as.data.frame(a$param)
s <- as.data.frame(a$stats)
w <- a$weights

names(p) <- c("rgc.tprob","ugc.tprob", "rct.tprob","uct.tprob",
              "syph.tprob", "hiv.rect.rr", "hiv.ureth.rr",
              "hiv.syph.rr", "syph.prim.sympt.prob",
              "syph.seco.sympt.prob", "syph.prim.sympt.prob.tx",
              "syph.seco.sympt.prob.tx")#,
              # "stianntest.gc.hivneg.coverage", "stianntest.gc.hivpos.coverage",
              # "stianntest.ct.hivneg.coverage", "stianntest.ct.hivpos.coverage",
              # "stianntest.syph.hivneg.coverage", "stianntest.syph.hivpos.coverage")
names(s) <- c("gc.incid", "ct.incid", "hiv.prev", "syph.incid",
              "syph.prev", "pssyph.prev",
              "gctest.hivneg", "gctest.hivpos", "cttest.hivneg",
              "cttest.hivpos", "syphtest.hivneg", "syphtest.hivpos",
              "gcslope", "ctslope", "syphslope", "hivslope",
              "hivprevslope", "syphprevslope")

comb <- cbind(s, p)
#View(comb)

( mean.s <- apply(s, 2, function(x) sum(x * w)) )
( mean.p <- apply(p, 2, function(x) sum(x * w)) )

# hist(s$prev.primsecosyph.hivpos / s$prev.primsecosyph.hivneg)

tar <- c(3.5, 5.6, 0.15, 1.5,
         0.02, 0.01,
         0.34, 0.472,
         0.344, 0.472,
         0.45, 0.69,
         0, 0, 0, 0, 0, 0)

par(mar = c(3,3,1,1), mgp = c(2,1,0), mfrow = c(3,3))
for (i in 1:ncol(s)) {
  hist(s[, i], col = "bisque2", border = "white", main = names(s)[i])
  abline(v = tar[i], lwd = 2, col = "red")
}

par(mar = c(3,3,1,1), mgp = c(2,1,0), mfrow = c(3,3))
for (i in 1:ncol(p)) {
  hist(p[, i], col = "bisque2", border = "white", main = names(p)[i])
  abline(v = mean.p[i], lwd = 2, col = "red")
}

for (i in seq_along(mean.p)) {
    assign(names(mean.p)[i], unname(mean.p[i]))
}

save(mean.p, file = "scripts/burnin/abc/abc.avg.parms.syph.100sim.1pct.rda")
save(mean.p, file = "est/abc.syph.parms.rda")

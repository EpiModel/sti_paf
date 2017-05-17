
## Analysis script for HIV risk compensation paper

rm(list = ls())
library("EpiModelHIV")
library("EpiModelHPC")
library("dplyr")
source("analysis/fx.R")

# Figures -------------------------------------------------------------

load("data/sim.n1300.rda")
sim.base <- truncate_sim(sim, at = 2600)

list.files("data/")
fn <- list.files("data/", full.names = TRUE)[2:5]
pal <- viridis(n = length(fn), option = "D")

pdf(file = "analysis/GC-CT-PrevCurves-byPrEPCov.pdf", h = 5, w = 12)
par(mar = c(3,3,1,1), mgp = c(2,1,0), mfrow = c(1,3))

for (i in seq_along(fn)) {
  load(fn[i])
  sim <- truncate_sim(sim, at = 2600)
  plot(sim, y = "ir100", mean.col = pal[i], qnts = FALSE, qnts.col = pal[i],
       add = i > 1, ylim = c(0, 4), main = "HIV Incidence by HIV PrEP Coverage")
}
legend("topright", col = pal, legend = c("5%", "10%", "20%", "40%"),
       lwd = 3, bty = "n")

for (i in seq_along(fn)) {
  load(fn[i])
  sim <- truncate_sim(sim, at = 2600)
  plot(sim, y = "ir100.gc", mean.col = pal[i], qnts = FALSE, qnts.col = pal[i],
       add = i > 1, ylim = c(0, 6), main = "GC Incidence by HIV PrEP Coverage")
}
legend("topright", col = pal, legend = c("5%", "10%", "20%", "40%"),
       lwd = 3, bty = "n")

pal <- viridis(n = length(fn), option = "D")
for (i in seq_along(fn)) {
  load(fn[i])
  sim <- truncate_sim(sim, at = 2600)
  plot(sim, y = "ir100.ct", mean.col = pal[i], qnts = FALSE, qnts.col = pal[i],
       add = i > 1, ylim = c(0, 8), main = "CT Incidence by HIV PrEP Coverage")
}
legend("topright", col = pal, legend = c("5%", "10%", "20%", "40%"),
       lwd = 3, bty = "n")
dev.off()


fn <- list.files("data/", full.names = TRUE)[c(1, 4, 6, 8)]

for (i in seq_along(fn)) {
  load(fn[i])

  nsim <- sim$control$nsims
  cov <- rep(sim$param$prep.coverage, nsim)
  scint <- rep(sim$param$prep.sti.screen.int, nsim)
  rc <- rep(sim$param$rcomp.prob, nsim)
  dft <- data.frame(cov, scint, rc)

  epi <- truncate_sim(sim, at = 2600)$epi

  dft <- within(dft, {
    ir100.gc <- unname(colMeans(tail(epi$ir100.gc, 52)))
    ir100.ct <- unname(colMeans(tail(epi$ir100.ct, 52)))
    ir100 <- unname(colMeans(tail(epi$ir100, 52)))
    tx.gc <- unname(colMeans(tail(epi$txGC, 52)))
    tx.ct <- unname(colMeans(tail(epi$txCT, 52)))
    hiv.prev <- as.numeric(tail(epi$i.prev, 1))
  })

  if (i == 1) {
    df <- dft
  } else {
    df <- rbind(df, dft)
  }
  cat("*")
}


df$scint <- ifelse(df$scint == 26 & df$cov == 0, 0, df$scint)

table(df$scint)

par(mfrow=c(1,1), mar = c(3,3,2,1), mgp = c(2,1,0))
boxplot(ir100.gc ~ scint, data = df, outline = FALSE, col = viridis(5, option = "D")[2:5],
        ylab = "IR per 100 PYAR", xlab = "PrEP-related Screening Interval",
        main = "GC Incidence by PrEP-related STI Screening Interval")
boxplot(ir100.ct ~ scint, data = df, outline = FALSE, col = viridis(5, option = "D")[2:5],
        ylab = "IR per 100 PYAR", xlab = "PrEP-related Screening Interval",
        main = "CT Incidence by PrEP-related STI Screening Interval")



load("data/sim.n1300.rda")
sim.base <- sim
load("data/sim.n1305.rda")
sim.comp <- sim
mo3 <- gc_pia_vec(sim.base, sim.comp)

load("data/sim.n1303.rda")
sim.comp <- sim
mo6 <- gc_pia_vec(sim.base, sim.comp)

load("data/sim.n1307.rda")
sim.comp <- sim
mo12 <- gc_pia_vec(sim.base, sim.comp)

df <- data.frame(mo3, mo6, mo12)
boxplot(df, outline = FALSE, col = )

epi_stats(sim.base, sim.comp)



# Treatment stats -----------------------------------------------------

system("scp hyak:/gscratch/csde/kweiss2/sti/data/sim.n1403.rda data/")

load("data/sim.n1402.rda")
sim.base <- truncate_sim(sim, at = 2600)
names(sim.base$epi)
head(sim.base$epi$prop.GC.asympt.tx)
tail(sim.base$epi$prop.GC.asympt.tx)
rowMeans(sim.base$epi$prop.GC.asympt.tx, na.rm = TRUE)

head(sim.base$epi$prop.CT.asympt.tx)
tail(sim.base$epi$prop.CT.asympt.tx)

head(sim.base$epi$prop.rGC.tx)
rowMeans(sim.base$epi$prop.rGC.tx)
tail(sim.base$epi$prop.rGC.tx)

is.nan(sim.base$epi$prop.rGC.tx$sim94[181])



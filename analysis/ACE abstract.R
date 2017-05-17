## STI PrEP Table 2
# Varying coverage of annual and high-risk testing

rm(list = ls())
library("EpiModelHIV")
library("EpiModelHPC")
library("dplyr")
source("analysis/fx.R")

# Base - No annual or high-risk
load("data/followup/sim.n3000.rda")
sim.base <- sim
epi_stats(sim.base, at = 520, qnt.low = 0.25, qnt.high = 0.75)

# Newer way
# Base - No annual or high-risk
haz <- as.numeric(colMeans(tail(sim.base$epi$ir100, 52)))
ir.base <- unname(colMeans(sim.base$epi$ir100)) * 1000
incid.base <- unname(colSums(sim.base$epi$incid))

haz.gc <- as.numeric(colMeans(tail(sim.base$epi$ir100.gc, 52)))
ir.base.gc <- unname(colMeans(sim.base$epi$ir100.gc)) * 1000
incid.base.gc <- unname(colSums(sim.base$epi$incid.gc))

haz.ct <- as.numeric(colMeans(tail(sim.base$epi$ir100.ct, 52)))
ir.base.ct <- unname(colMeans(sim.base$epi$ir100.ct)) * 1000
incid.base.ct <- unname(colSums(sim.base$epi$incid.ct))

haz.syph <- as.numeric(colMeans(tail(sim.base$epi$ir100.syph, 52)))
ir.base.syph <- unname(colMeans(sim.base$epi$ir100.syph)) * 1000
incid.base.syph <- unname(colSums(sim.base$epi$incid.syph))

sims <- c(3000, 3054, 3136, 3076, 3014)

annint <- rep(NA, length(sims))
hrint <- rep(NA, length(sims))

anncov <- rep(NA, length(sims))
hrcov <- rep(NA, length(sims))

gc.pia <- rep(NA, length(sims))
gc.nnt <- rep(NA, length(sims))

ct.pia <- rep(NA, length(sims))
ct.nnt <- rep(NA, length(sims))

syph.pia <- rep(NA, length(sims))
syph.nnt <- rep(NA, length(sims))

# add sims to data frame as an object?
df <- data.frame(annint, hrint, anncov, hrcov, gc.pia, ct.pia, syph.pia, gc.nnt, ct.nnt, syph.nnt)

for (i in seq_along(sims)) {
    
    fn <- list.files("data/followup/", pattern = as.character(sims[i]), full.names = TRUE)
    load(fn)
    
    #sim <- truncate_sim(sim, at = 2600)
    mn <- as.data.frame(sim)
    
    df$annint[i] <- sim$param$stitest.active.int
    df$hrint[i] <- sim$param$sti.highrisktest.int
    df$anncov[i] <- sim$param$stianntest.coverage
    df$hrcov[i] <- sim$param$stihighrisktest.coverage
    
    # Incidence Rate
    ir <- (colSums(sim$epi$incid, na.rm = TRUE)) /
        sum((1 - mn$i.prev)  * mn$num) * 52 * 1e5
    ir.gc <- (colSums(sim$epi$incid.gc, na.rm = TRUE)) /
        sum((1 - mn$prev.gc)  * mn$num) * 52 * 1e5
    ir.ct <- (colSums(sim$epi$incid.ct, na.rm = TRUE)) /
        sum((1 - mn$prev.ct)  * mn$num) * 52 * 1e5
    ir.syph <- (colSums(sim$epi$incid.syph, na.rm = TRUE)) /
        sum((1 - mn$prev.syph)  * mn$num) * 52 * 1e5
    
    vec.ir.hiv <- unname(colMeans(tail(sim$epi$ir100, 52)))
    vec.ir.gc <- unname(colMeans(tail(sim$epi$ir100.gc, 52)))
    vec.ir.ct <- unname(colMeans(tail(sim$epi$ir100.ct, 52)))
    vec.ir.syph <- unname(colMeans(tail(sim$epi$ir100.syph, 52)))
    
    ir.comp.gc <- unname(colMeans(sim$epi$ir100.gc)) * 1000
    vec.nia.gc <- round(ir.base.gc - ir.comp.gc, 1)
    vec.pia.gc <- vec.nia.gc/ir.base.gc
    vec.pia.gc <- vec.pia.gc[vec.pia.gc > -Inf]
    df$gc.pia[i] <- quantile(vec.pia.gc, probs = 0.50, na.rm = TRUE, names = FALSE)
    
    ir.comp.ct <- unname(colMeans(sim$epi$ir100.ct)) * 1000
    vec.nia.ct <- round(ir.base.ct - ir.comp.ct, 1)
    vec.pia.ct <- vec.nia.ct/ir.base.ct
    vec.pia.ct <- vec.pia.ct[vec.pia.ct > -Inf]
    df$ct.pia[i] <- quantile(vec.pia.ct, probs = 0.50, na.rm = TRUE, names = FALSE)
    
    ir.comp.syph <- unname(colMeans(sim$epi$ir100.syph)) * 1000
    vec.nia.syph <- round(ir.base.syph - ir.comp.syph, 1)
    vec.pia.syph <- vec.nia.syph/ir.base.syph
    vec.pia.syph <- vec.pia.syph[vec.pia.syph > -Inf]
    df$syph.pia[i] <- quantile(vec.pia.syph, probs = 0.50, na.rm = TRUE, names = FALSE)
    
    #NNT
    total.hiv.tests <- unname(tail(sim$epi$totalhivtests, 1))
    gc.asympt.tests <- unname(colMeans(tail(sim$epi$totalGCasympttests, 1)))
    ct.asympt.tests <- unname(colMeans(tail(sim$epi$totalCTasympttests, 1)))
    syph.asympt.tests <- unname(colMeans(tail(sim$epi$totalsyphasympttests, 1)))
    total.asympt.tests <- unname(colMeans(tail(sim$epi$totalstiasympttests, 1)))
    
    # HIV could be total HIV tests or total sti tests
    vec.hiv.nnt <- total.asympt.tests / (incid.base - unname(colSums(sim$epi$incid)))
    vec.gc.nnt <- gc.asympt.tests / (median(incid.base.gc) - unname(colSums(sim$epi$incid.gc)))
    vec.ct.nnt <- ct.asympt.tests / (median(incid.base.ct) - unname(colSums(sim$epi$incid.ct)))
    vec.syph.nnt <- syph.asympt.tests / (median(incid.base.syph) - unname(colSums(sim$epi$incid.syph)))
    
    df$gc.nnt[i] <- quantile(vec.gc.nnt, probs = 0.50, na.rm = TRUE, names = FALSE)
    df$ct.nnt[i] <- quantile(vec.ct.nnt, probs = 0.50, na.rm = TRUE, names = FALSE)
    df$syph.nnt[i] <- quantile(vec.syph.nnt, probs = 0.50, na.rm = TRUE, names = FALSE)
    
    cat("*")
    
}
View(df)

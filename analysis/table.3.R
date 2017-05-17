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

## Base STI lower-risk testing interval (364 days): n3054
## Varying STI lower-risk testing interval: n3131-n3141
## Base STI higher-risk testing interval: n3014
## Varying STI higher-risk testing interval: n3153 - n3173
sims <- c(3000, 3131:3141, 3054, 3153, 3156, 3158, 3160, 3162, 3164, 3166, 3168, 3170, 3172, 3014)

qnt.low <- 0.25
qnt.high <- 0.75

annint <- rep(NA, length(sims))
hrint <- rep(NA, length(sims))

hiv.incid <- rep(NA, length(sims))
hiv.incid.low <- rep(NA, length(sims))
hiv.incid.high <- rep(NA, length(sims))

hiv.hr <- rep(NA, length(sims))
hiv.hr.low <- rep(NA, length(sims))
hiv.hr.high <- rep(NA, length(sims))

hiv.pia <- rep(NA, length(sims))
hiv.pia.low <- rep(NA, length(sims))
hiv.pia.high <- rep(NA, length(sims))

hiv.nnt <- rep(NA, length(sims))
hiv.nnt.low <- rep(NA, length(sims))
hiv.nnt.high <- rep(NA, length(sims))

gc.incid <- rep(NA, length(sims))
gc.incid.low <- rep(NA, length(sims))
gc.incid.high <- rep(NA, length(sims))

gc.hr <- rep(NA, length(sims))
gc.hr.low <- rep(NA, length(sims))
gc.hr.high <- rep(NA, length(sims))

gc.pia <- rep(NA, length(sims))
gc.pia.low <- rep(NA, length(sims))
gc.pia.high <- rep(NA, length(sims))

gc.nnt <- rep(NA, length(sims))
gc.nnt.low <- rep(NA, length(sims))
gc.nnt.high <- rep(NA, length(sims))

ct.incid <- rep(NA, length(sims))
ct.incid.low <- rep(NA, length(sims))
ct.incid.high <- rep(NA, length(sims))

ct.hr <- rep(NA, length(sims))
ct.hr.low <- rep(NA, length(sims))
ct.hr.high <- rep(NA, length(sims))

ct.pia <- rep(NA, length(sims))
ct.pia.low <- rep(NA, length(sims))
ct.pia.high <- rep(NA, length(sims))

ct.nnt <- rep(NA, length(sims))
ct.nnt.low <- rep(NA, length(sims))
ct.nnt.high <- rep(NA, length(sims))

syph.incid <- rep(NA, length(sims))
syph.incid.low <- rep(NA, length(sims))
syph.incid.high <- rep(NA, length(sims))

syph.hr <- rep(NA, length(sims))
syph.hr.low <- rep(NA, length(sims))
syph.hr.high <- rep(NA, length(sims))

syph.pia <- rep(NA, length(sims))
syph.pia.low <- rep(NA, length(sims))
syph.pia.high <- rep(NA, length(sims))

syph.nnt <- rep(NA, length(sims))
syph.nnt.low <- rep(NA, length(sims))
syph.nnt.high <- rep(NA, length(sims))

# add sims to data frame as an object?
df <- data.frame(annint, hrint, hiv.incid.low, hiv.incid, hiv.incid.high, hiv.hr.low, hiv.hr, hiv.hr.high, 
                 hiv.pia.low, hiv.pia, hiv.pia.high, hiv.nnt.low, hiv.nnt, hiv.nnt.high,
                 gc.incid.low, gc.incid, gc.incid.high, gc.hr.low, gc.hr, gc.hr.high, 
                 gc.pia.low, gc.pia, gc.pia.high, gc.nnt.low, gc.nnt, gc.nnt.high,
                 ct.incid.low, ct.incid, ct.incid.high, ct.hr.low, ct.hr, ct.hr.high, 
                 ct.pia.low, ct.pia, ct.pia.high, ct.nnt.low, ct.nnt, ct.nnt.high,
                 syph.incid.low, syph.incid, syph.incid.high, syph.hr.low, syph.hr, syph.hr.high, 
                 syph.pia.low, syph.pia, syph.pia.high, syph.nnt.low, syph.nnt, syph.nnt.high)

for (i in seq_along(sims)) {
    
    fn <- list.files("data/followup/", pattern = as.character(sims[i]), full.names = TRUE)
    load(fn)
    
    #sim <- truncate_sim(sim, at = 2600)
    mn <- as.data.frame(sim)
    
    df$annint[i] <- sim$param$stitest.active.int
    df$hrint[i] <- sim$param$sti.highrisktest.int
    
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
    
    df$hiv.incid.low[i] <- quantile(vec.ir.hiv, probs = qnt.low, na.rm = TRUE, names = FALSE)
    df$hiv.incid[i] <- quantile(vec.ir.hiv, probs = 0.50, na.rm = TRUE, names = FALSE)
    df$hiv.incid.high[i] <- quantile(vec.ir.hiv, probs = qnt.high, na.rm = TRUE, names = FALSE)
    
    df$gc.incid.low[i] <- quantile(vec.ir.gc, probs = qnt.low, na.rm = TRUE, names = FALSE)
    df$gc.incid[i] <- quantile(vec.ir.gc, probs = 0.50, na.rm = TRUE, names = FALSE)
    df$gc.incid.high[i] <- quantile(vec.ir.gc, probs = qnt.high, na.rm = TRUE, names = FALSE)
    
    df$ct.incid.low[i] <- quantile(vec.ir.ct, probs = qnt.low, na.rm = TRUE, names = FALSE)
    df$ct.incid[i] <- quantile(vec.ir.ct, probs = 0.50, na.rm = TRUE, names = FALSE)
    df$ct.incid.high[i] <- quantile(vec.ir.ct, probs = qnt.high, na.rm = TRUE, names = FALSE)
        
    df$syph.incid.low[i] <- quantile(vec.ir.syph, probs = qnt.low, na.rm = TRUE, names = FALSE)
    df$syph.incid[i] <- quantile(vec.ir.syph, probs = 0.50, na.rm = TRUE, names = FALSE)
    df$syph.incid.high[i] <- quantile(vec.ir.syph, probs = qnt.high, na.rm = TRUE, names = FALSE)
    
    # HR
    num.hiv <- unname(colMeans(tail(sim$epi$ir100, 52)))
    denom.hiv <- unname(colMeans(tail(sim.base$epi$ir100, 52)))
    hr.vec.hiv <- num.hiv/denom.hiv
    hr.vec.hiv <- hr.vec.hiv[hr.vec.hiv < Inf]
    df$hiv.hr.low[i] <- quantile(hr.vec.hiv, probs = qnt.low, na.rm = TRUE, names = FALSE)
    df$hiv.hr[i] <- quantile(hr.vec.hiv, probs = 0.50, na.rm = TRUE, names = FALSE)
    df$hiv.hr.high[i] <- quantile(hr.vec.hiv, probs = qnt.high, na.rm = TRUE, names = FALSE)
    
    num.gc <- unname(colMeans(tail(sim$epi$ir100.gc, 52)))
    denom.gc <- unname(colMeans(tail(sim.base$epi$ir100.gc, 52)))
    hr.vec.gc <- num.gc/denom.gc
    hr.vec.gc <- hr.vec.gc[hr.vec.gc < Inf]
    df$gc.hr.low[i] <- quantile(hr.vec.gc, probs = qnt.low, na.rm = TRUE, names = FALSE)
    df$gc.hr[i] <- quantile(hr.vec.gc, probs = 0.50, na.rm = TRUE, names = FALSE)
    df$gc.hr.high[i] <- quantile(hr.vec.gc, probs = qnt.high, na.rm = TRUE, names = FALSE)
    
    num.ct <- unname(colMeans(tail(sim$epi$ir100.ct, 52)))
    denom.ct <- unname(colMeans(tail(sim.base$epi$ir100.ct, 52)))
    hr.vec.ct <- num.ct/denom.ct
    hr.vec.ct <- hr.vec.ct[hr.vec.ct < Inf]
    df$ct.hr.low[i] <- quantile(hr.vec.ct, probs = qnt.low, na.rm = TRUE, names = FALSE)
    df$ct.hr[i] <- quantile(hr.vec.ct, probs = 0.50, na.rm = TRUE, names = FALSE)
    df$ct.hr.high[i] <- quantile(hr.vec.ct, probs = qnt.high, na.rm = TRUE, names = FALSE)
    
    num.syph <- unname(colMeans(tail(sim$epi$ir100.syph, 52)))
    denom.syph <- unname(colMeans(tail(sim.base$epi$ir100.syph, 52)))
    hr.vec.syph <- num.syph/denom.syph
    hr.vec.syph <- hr.vec.syph[hr.vec.syph < Inf]
    df$syph.hr.low[i] <- quantile(hr.vec.syph, probs = qnt.low, na.rm = TRUE, names = FALSE)
    df$syph.hr[i] <- quantile(hr.vec.syph, probs = 0.50, na.rm = TRUE, names = FALSE)
    df$syph.hr.high[i] <- quantile(hr.vec.syph, probs = qnt.high, na.rm = TRUE, names = FALSE)
    
    #PIA
    ir.comp <- unname(colMeans(sim$epi$ir100)) * 1000
    vec.nia.hiv <- round(ir.base - ir.comp, 1)
    vec.pia.hiv <- vec.nia.hiv/ir.base
    vec.pia.hiv <- vec.pia.hiv[vec.pia.hiv > -Inf]
    df$hiv.pia.low[i] <- quantile(vec.pia.hiv, probs = qnt.low, na.rm = TRUE, names = FALSE)
    df$hiv.pia[i] <- quantile(vec.pia.hiv, probs = 0.50, na.rm = TRUE, names = FALSE)
    df$hiv.pia.high[i] <- quantile(vec.pia.hiv, probs = qnt.high, na.rm = TRUE, names = FALSE)
    
    ir.comp.gc <- unname(colMeans(sim$epi$ir100.gc)) * 1000
    vec.nia.gc <- round(ir.base.gc - ir.comp.gc, 1)
    vec.pia.gc <- vec.nia.gc/ir.base.gc
    vec.pia.gc <- vec.pia.gc[vec.pia.gc > -Inf]
    df$gc.pia.low[i] <- quantile(vec.pia.gc, probs = qnt.low, na.rm = TRUE, names = FALSE)
    df$gc.pia[i] <- quantile(vec.pia.gc, probs = 0.50, na.rm = TRUE, names = FALSE)
    df$gc.pia.high[i] <- quantile(vec.pia.gc, probs = qnt.high, na.rm = TRUE, names = FALSE)
    
    ir.comp.ct <- unname(colMeans(sim$epi$ir100.ct)) * 1000
    vec.nia.ct <- round(ir.base.ct - ir.comp.ct, 1)
    vec.pia.ct <- vec.nia.ct/ir.base.ct
    vec.pia.ct <- vec.pia.ct[vec.pia.ct > -Inf]
    df$ct.pia.low[i] <- quantile(vec.pia.ct, probs = qnt.low, na.rm = TRUE, names = FALSE)
    df$ct.pia[i] <- quantile(vec.pia.ct, probs = 0.50, na.rm = TRUE, names = FALSE)
    df$ct.pia.high[i] <- quantile(vec.pia.ct, probs = qnt.high, na.rm = TRUE, names = FALSE)
    
    ir.comp.syph <- unname(colMeans(sim$epi$ir100.syph)) * 1000
    vec.nia.syph <- round(ir.base.syph - ir.comp.syph, 1)
    vec.pia.syph <- vec.nia.syph/ir.base.syph
    vec.pia.syph <- vec.pia.syph[vec.pia.syph > -Inf]
    df$syph.pia.low[i] <- quantile(vec.pia.syph, probs = qnt.low, na.rm = TRUE, names = FALSE)
    df$syph.pia[i] <- quantile(vec.pia.syph, probs = 0.50, na.rm = TRUE, names = FALSE)
    df$syph.pia.high[i] <- quantile(vec.pia.syph, probs = qnt.high, na.rm = TRUE, names = FALSE)
    
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
    
    df$hiv.nnt.low[i] <- quantile(vec.hiv.nnt, probs = qnt.low, na.rm = TRUE, names = FALSE)
    df$hiv.nnt[i] <- quantile(vec.hiv.nnt, probs = 0.50, na.rm = TRUE, names = FALSE)
    df$hiv.nnt.high[i] <- quantile(vec.hiv.nnt, probs = qnt.high, na.rm = TRUE, names = FALSE)
    
    df$gc.nnt.low[i] <- quantile(vec.gc.nnt, probs = qnt.low, na.rm = TRUE, names = FALSE)
    df$gc.nnt[i] <- quantile(vec.gc.nnt, probs = 0.50, na.rm = TRUE, names = FALSE)
    df$gc.nnt.high[i] <- quantile(vec.gc.nnt, probs = qnt.high, na.rm = TRUE, names = FALSE)
    
    df$ct.nnt.low[i] <- quantile(vec.ct.nnt, probs = qnt.low, na.rm = TRUE, names = FALSE)
    df$ct.nnt[i] <- quantile(vec.ct.nnt, probs = 0.50, na.rm = TRUE, names = FALSE)
    df$ct.nnt.high[i] <- quantile(vec.ct.nnt, probs = qnt.high, na.rm = TRUE, names = FALSE)
    
    df$syph.nnt.low[i] <- quantile(vec.syph.nnt, probs = qnt.low, na.rm = TRUE, names = FALSE)
    df$syph.nnt[i] <- quantile(vec.syph.nnt, probs = 0.50, na.rm = TRUE, names = FALSE)
    df$syph.nnt.high[i] <- quantile(vec.syph.nnt, probs = qnt.high, na.rm = TRUE, names = FALSE)
    
    cat("*")
    
}

write.csv(df, "C:/Users/kweiss2/Documents/GitHub/stitestguidelines/analysis/Table 3.csv")



# Older way
# Varying Lower-Risk
# 3131:3141, 3054, Lower-Risk = 40%, 28 days to 364 days by ~28 days, Higher-Risk = 0%, 182 days
# 28, 63, 91, 119, 147, 182, 210, 238, 273, 301, 329, 364

load("data/followup/sim.n3131.rda")
sim$param$stitest.active.int
epi_stats(sim.base, sim, at = 520, qnt.low = 0.25, qnt.high = 0.75)

load("data/followup/sim.n3132.rda")
sim$param$stitest.active.int
epi_stats(sim.base, sim, at = 520, qnt.low = 0.25, qnt.high = 0.75)

load("data/followup/sim.n3133.rda")
sim$param$stitest.active.int
epi_stats(sim.base, sim, at = 520, qnt.low = 0.25, qnt.high = 0.75)

load("data/followup/sim.n3134.rda")
sim$param$stitest.active.int
epi_stats(sim.base, sim, at = 520, qnt.low = 0.25, qnt.high = 0.75)

load("data/followup/sim.n3135.rda")
sim$param$stitest.active.int
epi_stats(sim.base, sim, at = 520, qnt.low = 0.25, qnt.high = 0.75)

load("data/followup/sim.n3136.rda")
sim$param$stitest.active.int
epi_stats(sim.base, sim, at = 520, qnt.low = 0.25, qnt.high = 0.75)

load("data/followup/sim.n3137.rda")
sim$param$stitest.active.int
epi_stats(sim.base, sim, at = 520, qnt.low = 0.25, qnt.high = 0.75)

load("data/followup/sim.n3138.rda")
sim$param$stitest.active.int
epi_stats(sim.base, sim, at = 520, qnt.low = 0.25, qnt.high = 0.75)

load("data/followup/sim.n3139.rda")
sim$param$stitest.active.int
epi_stats(sim.base, sim, at = 520, qnt.low = 0.25, qnt.high = 0.75)

load("data/followup/sim.n3140.rda")
sim$param$stitest.active.int
epi_stats(sim.base, sim, at = 520, qnt.low = 0.25, qnt.high = 0.75)

load("data/followup/sim.n3141.rda")
sim$param$stitest.active.int
epi_stats(sim.base, sim, at = 520, qnt.low = 0.25, qnt.high = 0.75)

load("data/followup/sim.n3054.rda")
sim$param$stitest.active.int
epi_stats(sim.base, sim, at = 520, qnt.low = 0.25, qnt.high = 0.75)

# Varying Higher-Risk
# 3153, 3156, 3158, 3160, 3162, 3164, 3166, 3168, 3170, 3172, 3014, Annual = 0%, 364 days, Higher-Risk = 40%, 182 days
# 28, 49, 63, 77, 91, 105, 119, 133, 147, 161, 182 

load("data/followup/sim.n3153.rda")
sim$param$sti.highrisktest.int
epi_stats(sim.base, sim, at = 520, qnt.low = 0.25, qnt.high = 0.75)

load("data/followup/sim.n3156.rda")
sim$param$sti.highrisktest.int
epi_stats(sim.base, sim, at = 520, qnt.low = 0.25, qnt.high = 0.75)

load("data/followup/sim.n3158.rda")
sim$param$sti.highrisktest.int
epi_stats(sim.base, sim, at = 520, qnt.low = 0.25, qnt.high = 0.75)

load("data/followup/sim.n3160.rda")
sim$param$sti.highrisktest.int
epi_stats(sim.base, sim, at = 520, qnt.low = 0.25, qnt.high = 0.75)

load("data/followup/sim.n3162.rda")
sim$param$sti.highrisktest.int
epi_stats(sim.base, sim, at = 520, qnt.low = 0.25, qnt.high = 0.75)

load("data/followup/sim.n3164.rda")
sim$param$sti.highrisktest.int
epi_stats(sim.base, sim, at = 520, qnt.low = 0.25, qnt.high = 0.75)

load("data/followup/sim.n3166.rda")
sim$param$sti.highrisktest.int
epi_stats(sim.base, sim, at = 520, qnt.low = 0.25, qnt.high = 0.75)

load("data/followup/sim.n3168.rda")
sim$param$sti.highrisktest.int
epi_stats(sim.base, sim, at = 520, qnt.low = 0.25, qnt.high = 0.75)

load("data/followup/sim.n3170.rda")
sim$param$sti.highrisktest.int
epi_stats(sim.base, sim, at = 520, qnt.low = 0.25, qnt.high = 0.75)

load("data/followup/sim.n3172.rda")
sim$param$sti.highrisktest.int
epi_stats(sim.base, sim, at = 520, qnt.low = 0.25, qnt.high = 0.75)

load("data/followup/sim.n3014.rda")
sim$param$sti.highrisktest.int
epi_stats(sim.base, sim, at = 520, qnt.low = 0.25, qnt.high = 0.75)



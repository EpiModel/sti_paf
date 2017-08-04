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
# Varying Lower-Risk
# 3014, 3025, 3036, 3047, 3058, 3069, 3080, 3091, 3102, 3113, 3124: Annual = 0.1 - 1.0 by 0.1, 364 days, HR = 40%, 182 days

# Newer way:
sims <- c(3014, 3025, 3036, 3047, 3058, 3069, 3080, 3091, 3102, 3113, 3124, 3054:3064)

qnt.low <- 0.25
qnt.high <- 0.75

anncov <- rep(NA, length(sims))
hrcov <- rep(NA, length(sims))

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
df <- data.frame(anncov, hrcov, hiv.incid.low, hiv.incid, hiv.incid.high, hiv.hr.low, hiv.hr, hiv.hr.high, 
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
    denom.hiv <- unname(colMeans(sim.base$epi$ir100)) * 1000
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
    
    #HIV could be HIV tests or total STI tests
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

write.csv(df, "C:/Users/kweiss2/Documents/GitHub/stitestguidelines/analysis/Table 2.csv") 



# Older Way
load("data/followup/sim.n3014.rda")
sim$param$stianntest.coverage
sim$param$stihighrisktest.coverage
epi_stats(sim.base, sim, at = 520, qnt.low = 0.25, qnt.high = 0.75)

load("data/followup/sim.n3025.rda")
sim$param$stianntest.coverage
sim$param$stihighrisktest.coverage
epi_stats(sim.base, sim, at = 520, qnt.low = 0.25, qnt.high = 0.75)

load("data/followup/sim.n3036.rda")
sim$param$stianntest.coverage
sim$param$stihighrisktest.coverage
epi_stats(sim.base, sim, at = 520, qnt.low = 0.25, qnt.high = 0.75)

load("data/followup/sim.n3047.rda")
sim$param$stianntest.coverage
sim$param$stihighrisktest.coverage
epi_stats(sim.base, sim, at = 520, qnt.low = 0.25, qnt.high = 0.75)

load("data/followup/sim.n3058.rda")
sim$param$stianntest.coverage
sim$param$stihighrisktest.coverage
epi_stats(sim.base, sim, at = 520, qnt.low = 0.25, qnt.high = 0.75)

load("data/followup/sim.n3069.rda")
sim$param$stianntest.coverage
sim$param$stihighrisktest.coverage
epi_stats(sim.base, sim, at = 520, qnt.low = 0.25, qnt.high = 0.75)

load("data/followup/sim.n3080.rda")
sim$param$stianntest.coverage
sim$param$stihighrisktest.coverage
epi_stats(sim.base, sim, at = 520, qnt.low = 0.25, qnt.high = 0.75)

load("data/followup/sim.n3091.rda")
sim$param$stianntest.coverage
sim$param$stihighrisktest.coverage
epi_stats(sim.base, sim, at = 520, qnt.low = 0.25, qnt.high = 0.75)

load("data/followup/sim.n3102.rda")
sim$param$stianntest.coverage
sim$param$stihighrisktest.coverage
epi_stats(sim.base, sim, at = 520, qnt.low = 0.25, qnt.high = 0.75)

load("data/followup/sim.n3113.rda")
sim$param$stianntest.coverage
sim$param$stihighrisktest.coverage
epi_stats(sim.base, sim, at = 520, qnt.low = 0.25, qnt.high = 0.75)

load("data/followup/sim.n3124.rda")
sim$param$stianntest.coverage
sim$param$stihighrisktest.coverage
epi_stats(sim.base, sim, at = 520, qnt.low = 0.25, qnt.high = 0.75)

# Varying Higher-Risk
# 3054 - 3064 Annual = 40%, 364 days, HR = 0.0 - 1.0 by 0.1, 182 days
load("data/followup/sim.n3054.rda")
sim$param$stianntest.coverage
sim$param$stihighrisktest.coverage
epi_stats(sim.base, sim, at = 520, qnt.low = 0.25, qnt.high = 0.75)

load("data/followup/sim.n3055.rda")
sim$param$stianntest.coverage
sim$param$stihighrisktest.coverage
epi_stats(sim.base, sim, at = 520, qnt.low = 0.25, qnt.high = 0.75)

load("data/followup/sim.n3056.rda")
sim$param$stianntest.coverage
sim$param$stihighrisktest.coverage
epi_stats(sim.base, sim, at = 520, qnt.low = 0.25, qnt.high = 0.75)

load("data/followup/sim.n3057.rda")
sim$param$stianntest.coverage
sim$param$stihighrisktest.coverage
epi_stats(sim.base, sim, at = 520, qnt.low = 0.25, qnt.high = 0.75)

## Same as above scenario
load("data/followup/sim.n3058.rda")
sim$param$stianntest.coverage
sim$param$stihighrisktest.coverage
epi_stats(sim.base, sim, at = 520, qnt.low = 0.25, qnt.high = 0.75)

load("data/followup/sim.n3059.rda")
sim$param$stianntest.coverage
sim$param$stihighrisktest.coverage
epi_stats(sim.base, sim, at = 520, qnt.low = 0.25, qnt.high = 0.75)

load("data/followup/sim.n3060.rda")
sim$param$stianntest.coverage
sim$param$stihighrisktest.coverage
epi_stats(sim.base, sim, at = 520, qnt.low = 0.25, qnt.high = 0.75)

load("data/followup/sim.n3061.rda")
sim$param$stianntest.coverage
sim$param$stihighrisktest.coverage
epi_stats(sim.base, sim, at = 520, qnt.low = 0.25, qnt.high = 0.75)

load("data/followup/sim.n3062.rda")
sim$param$stianntest.coverage
sim$param$stihighrisktest.coverage
epi_stats(sim.base, sim, at = 520, qnt.low = 0.25, qnt.high = 0.75)

load("data/followup/sim.n3063.rda")
sim$param$stianntest.coverage
sim$param$stihighrisktest.coverage
epi_stats(sim.base, sim, at = 520, qnt.low = 0.25, qnt.high = 0.75)

load("data/followup/sim.n3064.rda")
sim$param$stianntest.coverage
sim$param$stihighrisktest.coverage
epi_stats(sim.base, sim, at = 520, qnt.low = 0.25, qnt.high = 0.75)


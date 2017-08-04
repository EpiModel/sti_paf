## STI Testing Guidelines Table 1
# Varying Indications for High-Risk Testing

rm(list = ls())
library("EpiModelHIV")
library("EpiModelHPC")
library("dplyr")
source("analysis/fx.R")

# Base - No testing
load("data/followup/sim.n5001.rda")
sim.base <- sim
epi_stats(sim.base, at = 520, qnt.low = 0.25, qnt.high = 0.75)

haz <- as.numeric(colMeans(tail(sim.base$epi$ir100, 52)))
ir.base <- unname(colMeans(sim.base$epi$ir100)) * 1000
incid.base <- unname(colSums(sim.base$epi$incid))

## Varying Indications:

# Newer way:
sims <- c(5001:5030)

qnt.low <- 0.25
qnt.high <- 0.75

hiv.incid <- rep(NA, length(sims))
hiv.incid.low <- rep(NA, length(sims))
hiv.incid.high <- rep(NA, length(sims))

hiv.hr <- rep(NA, length(sims))
hiv.hr.low <- rep(NA, length(sims))
hiv.hr.high <- rep(NA, length(sims))


# add sims to data frame as an object?
df <- data.frame(hiv.incid.low, hiv.incid, hiv.incid.high, hiv.hr.low, hiv.hr, hiv.hr.high)

for (i in seq_along(sims)) {

    fn <- list.files("data/followup/", pattern = as.character(sims[i]), full.names = TRUE)
    load(fn)
    
    #sim <- truncate_sim(sim, at = 2600)
    mn <- as.data.frame(sim)
    
    # Incidence Rate
    ir.base <- (colSums(sim$epi$incid, na.rm = TRUE)) /
        sum((1 - mn$i.prev)  * mn$num) * 52 * 1e5
    
    vec.ir.hiv <- unname(colMeans(tail(sim$epi$ir100, 52)))

    df$hiv.incid.low[i] <- quantile(vec.ir.hiv, probs = qnt.low, na.rm = TRUE, names = FALSE)
    df$hiv.incid[i] <- quantile(vec.ir.hiv, probs = 0.50, na.rm = TRUE, names = FALSE)
    df$hiv.incid.high[i] <- quantile(vec.ir.hiv, probs = qnt.high, na.rm = TRUE, names = FALSE)
    
    # HR
    num.hiv <- unname(colMeans(tail(sim$epi$ir100, 52)))
    denom.hiv <- unname(colMeans(sim.base$epi$ir100, 52))
    hr.vec.hiv <- num.hiv/denom.hiv
    hr.vec.hiv <- hr.vec.hiv[hr.vec.hiv < Inf]
    df$hiv.hr.low[i] <- quantile(hr.vec.hiv, probs = qnt.low, na.rm = TRUE, names = FALSE)
    df$hiv.hr[i] <- quantile(hr.vec.hiv, probs = 0.50, na.rm = TRUE, names = FALSE)
    df$hiv.hr.high[i] <- quantile(hr.vec.hiv, probs = qnt.high, na.rm = TRUE, names = FALSE)

    cat("*")
    
}
df
names(df$hiv.incid.low) <- names(df$hiv.incid) <- names(df$hiv.incid.high) <- names(df$hiv.hr.low) <- names(df$hiv.hr) <- names(df$hiv.hr.high) <-
    c("5001",
      "5002",
      "5003",
      "5004",
      "5005",
      "5006",
      "5007",
      "5008",
      "5009",
      "5010",
      "5011",
      "5012",
      "5013",
      "5014",
      "5015",
      "5016",
      "5017",
      "5018",
      "5019",
      "5020",
      "5021",
      "5022",
      "5023",
      "5024",
      "5025",
      "5026",
      "5027",
      "5028",
      "5029",
      "5030")

df
write.csv(df, "C:/Users/Jeb/Dropbox/Hazard Ratios.csv")


# Older way:

# 3142 - STI
load("data/followup/sim.n3142.rda")
epi_stats(sim.base, sim, at = 520, qnt.low = 0.25, qnt.high = 0.75)

# 3143 - recent partners
load("data/followup/sim.n3143.rda")
epi_stats(sim.base, sim, at = 520, qnt.low = 0.25, qnt.high = 0.75)

# 3144 - new partners
load("data/followup/sim.n3144.rda")
epi_stats(sim.base, sim, at = 520, qnt.low = 0.25, qnt.high = 0.75)

# 3145 - partner who has multiple partners
load("data/followup/sim.n3145.rda")
epi_stats(sim.base, sim, at = 520, qnt.low = 0.25, qnt.high = 0.75)

# 3146 - partner with a STI
load("data/followup/sim.n3146.rda")
epi_stats(sim.base, sim, at = 520, qnt.low = 0.25, qnt.high = 0.75)

# 3147 - any CAI in a non-main
load("data/followup/sim.n3147.rda")
epi_stats(sim.base, sim, at = 520, qnt.low = 0.25, qnt.high = 0.75)

# 3148 - any CAI
load("data/followup/sim.n3148.rda")
epi_stats(sim.base, sim, at = 520, qnt.low = 0.25, qnt.high = 0.75)

# 3149 - recent or new partners
load("data/followup/sim.n3149.rda")
epi_stats(sim.base, sim, at = 520, qnt.low = 0.25, qnt.high = 0.75)

# 3150 - sti, recent, or new partners
load("data/followup/sim.n3150.rda")
epi_stats(sim.base, sim, at = 520, qnt.low = 0.25, qnt.high = 0.75)

# 3151 - CAI in non-main or any CAI
load("data/followup/sim.n3151.rda")
epi_stats(sim.base, sim, at = 520, qnt.low = 0.25, qnt.high = 0.75)

# 3152 - partner with multiple partners or with a STI
load("data/followup/sim.n3152.rda")
epi_stats(sim.base, sim, at = 520, qnt.low = 0.25, qnt.high = 0.75)
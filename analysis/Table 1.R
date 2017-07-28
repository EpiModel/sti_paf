## STI PAF Table 1
# Sensitivity analysis with varying relative risks for each STI

rm(list = ls())
library("EpiModelHIV")
library("EpiModelHPC")
library("dplyr")
source("analysis/fx.R")



#Final year of follow-up - highest prob prevails
#Obtain incidence rate and simulation interval for each of the 30 combinations in Table 1
load("data/followup/sim.n5001.rda")
sim.base <- sim
epi_stats(sim.base, at = 520, qnt.low = 0.25, qnt.high = 0.75)

haz.base <- as.numeric(colMeans(tail(sim.base$epi$ir100, 52)))
ir.base <- unname(colMeans(sim.base$epi$ir100)) * 1000
incid.base <- unname(colSums(sim.base$epi$incid))

sims <- c(5001:5015)

qnt.low <- 0.25
qnt.high <- 0.75

hiv.incid <- rep(NA, length(sims))
hiv.incid.low <- rep(NA, length(sims))
hiv.incid.high <- rep(NA, length(sims))

hiv.hr <- rep(NA, length(sims))
hiv.hr.low <- rep(NA, length(sims))
hiv.hr.high <- rep(NA, length(sims))

#add sims to dataframe as object
df <- data.frame(hiv.incid.low, hiv.incid, hiv.incid.high, hiv.hr.low, hiv.hr, hiv.hr.high)
df
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
  
  num.hiv <- unname(colMeans(tail(sim$epi$ir100, 52)))
  denom.hiv <- unname(colMeans(tail(sim.base$epi$ir100, 52)))
  hr.vec.hiv <- num.hiv/denom.hiv
  hr.vec.hiv <- hr.vec.hiv[hr.vec.hiv < Inf]
  df$hiv.hr.low[i] <- quantile(hr.vec.hiv, probs = qnt.low, na.rm = TRUE, names = FALSE)
  df$hiv.hr[i] <- quantile(hr.vec.hiv, probs = 0.50, na.rm = TRUE, names = FALSE)
  df$hiv.hr.high[i] <- quantile(hr.vec.hiv, probs = qnt.high, na.rm = TRUE, names = FALSE)
  
  cat("*")
  
}

df

names(df$hiv.incid.low) <- names(df$hiv.incid) <- names(df$hiv.incid.high) <- 
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
      "5015")

df
write.csv(df, "C:/Users/jsjone2/Desktop/Table 1 Data - Final Year - Highest Prob.csv")

#Final year of follow-up - highest prob prevails
#Obtain incidence rate and simulation interval for each of the 30 combinations in Table 1
load("data/followup/sim.n5016.rda")
sim.base <- sim
epi_stats(sim.base, at = 520, qnt.low = 0.25, qnt.high = 0.75)

haz.base <- as.numeric(colMeans(tail(sim.base$epi$ir100, 52)))
ir.base <- unname(colMeans(sim.base$epi$ir100)) * 1000
incid.base <- unname(colSums(sim.base$epi$incid))

sims <- c(5016:5030)

qnt.low <- 0.25
qnt.high <- 0.75

hiv.incid <- rep(NA, length(sims))
hiv.incid.low <- rep(NA, length(sims))
hiv.incid.high <- rep(NA, length(sims))

hiv.hr <- rep(NA, length(sims))
hiv.hr.low <- rep(NA, length(sims))
hiv.hr.high <- rep(NA, length(sims))

#add sims to dataframe as object
df <- data.frame(hiv.incid.low, hiv.incid, hiv.incid.high, hiv.hr.low, hiv.hr, hiv.hr.high)
df
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
  
  num.hiv <- unname(colMeans(tail(sim$epi$ir100, 52)))
  denom.hiv <- unname(colMeans(tail(sim.base$epi$ir100, 52)))
  hr.vec.hiv <- num.hiv/denom.hiv
  hr.vec.hiv <- hr.vec.hiv[hr.vec.hiv < Inf]
  df$hiv.hr.low[i] <- quantile(hr.vec.hiv, probs = qnt.low, na.rm = TRUE, names = FALSE)
  df$hiv.hr[i] <- quantile(hr.vec.hiv, probs = 0.50, na.rm = TRUE, names = FALSE)
  df$hiv.hr.high[i] <- quantile(hr.vec.hiv, probs = qnt.high, na.rm = TRUE, names = FALSE)
  
  cat("*")
  
}

df

names(df$hiv.incid.low) <- names(df$hiv.incid) <- names(df$hiv.incid.high) <- 
  c("5016",
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
write.csv(df, "C:/Users/jsjone2/Desktop/Table 1 Data - Final Year - Multiplicative.csv")


#entire follow-up
#Obtain incidence rate and simulation interval for each of the 30 combinations in Table 1
sims <- c(3003, 5001:5030)

qnt.low <- 0.25
qnt.high <- 0.75

hiv.incid <- rep(NA, length(sims))
hiv.incid.low <- rep(NA, length(sims))
hiv.incid.high <- rep(NA, length(sims))

#add sims to dataframe as object
df <- data.frame(hiv.incid.low, hiv.incid, hiv.incid.high)
df
for (i in seq_along(sims)) {
  
  fn <- list.files("data/followup/", pattern = as.character(sims[i]), full.names = TRUE)
  load(fn)
  
  #sim <- truncate_sim(sim, at = 2600)
  mn <- as.data.frame(sim)
  
  # Incidence Rate
  ir.base <- (colSums(sim$epi$incid, na.rm = TRUE)) /
    sum((1 - mn$i.prev)  * mn$num) * 52 * 1e5
  
  vec.ir.hiv <- unname(colMeans(tail(sim$epi$ir100, 519)))
  
  
  df$hiv.incid.low[i] <- quantile(vec.ir.hiv, probs = qnt.low, na.rm = TRUE, names = FALSE)
  df$hiv.incid[i] <- quantile(vec.ir.hiv, probs = 0.50, na.rm = TRUE, names = FALSE)
  df$hiv.incid.high[i] <- quantile(vec.ir.hiv, probs = qnt.high, na.rm = TRUE, names = FALSE)
  
  cat("*")
  
}

names(df$hiv.incid.low) <- names(df$hiv.incid) <- names(df$hiv.incid.high) <- 
  c("3003",
    "5001",
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
write.csv(df, "C:/Users/jsjone2/Desktop/Table 1 Data - Full Follow-Up.csv")
## STI PAF Table 1
# Sensitivity analysis with varying relative risks for each STI

rm(list = ls())
library("EpiModelHIV")
library("EpiModelHPC")
library("dplyr")
source("analysis/fx.R")



#Highest prob prevails
#Obtain incidence rate and simulation interval for each of the 30 combinations
sim5001 <- read.csv("data/directeffect/sim 5001 .csv")


haz.base <- as.numeric(mean(sim5001$ir100))
ir.base <- unname(mean(sim5001$ir100)) * 1000
incid.base <- unname(mean(sim5001$incid))

sims <- c(5001:5015)

qnt.low <- 0.025
qnt.high <- 0.975

scenario_num <- rep(NA, length(sims))

hiv.incid <- rep(NA, length(sims))
hiv.incid.low <- rep(NA, length(sims))
hiv.incid.high <- rep(NA, length(sims))

ir100 <- rep(NA, length(sims))
ir100.low <- rep(NA, length(sims))
ir100.high <- rep(NA, length(sims))

cell1.low <- rep(NA, length(sims))
cell1 <- rep(NA, length(sims))
cell1.high <- rep(NA, length(sims))

cell2.low <- rep(NA, length(sims))
cell2 <- rep(NA, length(sims))
cell2.high <- rep(NA, length(sims))

cell3.low <- rep(NA, length(sims))
cell3 <- rep(NA, length(sims))
cell3.high <- rep(NA, length(sims))

cell4.low <- rep(NA, length(sims))
cell4 <- rep(NA, length(sims))
cell4.high <- rep(NA, length(sims))

paf.low <- rep(NA, length(sims))
paf <- rep(NA, length(sims))
paf.high <- rep(NA, length(sims))

#add sims to dataframe as object
df <- data.frame(scenario_num, 
                 hiv.incid.low, hiv.incid, hiv.incid.high, 
                 ir100.low, ir100, ir100.high,
                 cell1.low, cell1, cell1.high, 
                 cell2.low, cell2, cell2.high,
                 cell3.low, cell3, cell3.high, 
                 cell4.low, cell4, cell4.high,
                 paf.low, paf, paf.high)
df
for (i in seq_along(sims)) {
  
  fn <- paste("data/directeffect/sim", pattern = as.character(sims[i]),".csv")
  
  
  #sim <- truncate_sim(sim, at = 2600)
  mn <- read.csv(fn)

  df$scenario_num[i] <- sims[i]
  
  #HIV diagnosis counts
  df$hiv.incid.low[i] <- quantile(mn$incid, probs = qnt.low, na.rm = TRUE, names = FALSE)
  df$hiv.incid[i] <- quantile(mn$incid, probs = 0.50, na.rm = TRUE, names = FALSE)
  df$hiv.incid.high[i] <- quantile(mn$incid, probs = qnt.high, na.rm = TRUE, names = FALSE)
    
  #Incidence Rate
  df$ir100.low[i] <- quantile(mn$ir100, probs=qnt.low, na.rm = TRUE, names = FALSE)
  df$ir100[i] <- quantile(mn$ir100, probs=0.50, na.rm = TRUE, names = FALSE)
  df$ir100.high[i] <- quantile(mn$ir100, probs=qnt.high, na.rm = TRUE, names = FALSE)
  
  #sum cell1 - cell4 for proportion
  cell_sum_low <- quantile(mn$cell1_sti, probs = qnt.low, na.rm = TRUE, names = FALSE) +
                  quantile(mn$cell2_sti, probs = qnt.low, na.rm = TRUE, names = FALSE) +
                  quantile(mn$cell3_sti, probs = qnt.low, na.rm = TRUE, names = FALSE) +
                  quantile(mn$cell4_sti, probs = qnt.low, na.rm = TRUE, names = FALSE) 
  
  cell_sum <- quantile(mn$cell1_sti, probs = 0.50, na.rm = TRUE, names = FALSE) +
              quantile(mn$cell2_sti, probs = 0.50, na.rm = TRUE, names = FALSE) +
              quantile(mn$cell3_sti, probs = 0.50, na.rm = TRUE, names = FALSE) +
              quantile(mn$cell4_sti, probs = 0.50, na.rm = TRUE, names = FALSE) 
  
  cell_sum_high <- quantile(mn$cell1_sti, probs = qnt.high, na.rm = TRUE, names = FALSE) +
                  quantile(mn$cell2_sti, probs = qnt.high, na.rm = TRUE, names = FALSE) +
                  quantile(mn$cell3_sti, probs = qnt.high, na.rm = TRUE, names = FALSE) +
                  quantile(mn$cell4_sti, probs = qnt.high, na.rm = TRUE, names = FALSE) 
  
  #Cell1 - Cell4
  df$cell1.low[i] <- (quantile(mn$cell1_sti, probs=qnt.low, na.rm = TRUE, names = FALSE) / cell_sum_low)*100
  df$cell1[i] <- (quantile(mn$cell1_sti, probs=0.50, na.rm = TRUE, names = FALSE) / cell_sum)*100
  df$cell1.high[i] <- (quantile(mn$cell1_sti, probs=qnt.high, na.rm = TRUE, names = FALSE) / cell_sum_high)*100
  
  df$cell2.low[i] <- (quantile(mn$cell2_sti, probs=qnt.low, na.rm = TRUE, names = FALSE) / cell_sum_low)*100
  df$cell2[i] <- (quantile(mn$cell2_sti, probs=0.50, na.rm = TRUE, names = FALSE) / cell_sum)*100
  df$cell2.high[i] <- (quantile(mn$cell2_sti, probs=qnt.high, na.rm = TRUE, names = FALSE) / cell_sum_high)*100
  
  df$cell3.low[i] <- (quantile(mn$cell3_sti, probs=qnt.low, na.rm = TRUE, names = FALSE) / cell_sum_low)*100
  df$cell3[i] <- (quantile(mn$cell3_sti, probs=0.50, na.rm = TRUE, names = FALSE) / cell_sum)*100
  df$cell3.high[i] <- (quantile(mn$cell3_sti, probs=qnt.high, na.rm = TRUE, names = FALSE) / cell_sum_high)*100
  
  df$cell4.low[i] <- (quantile(mn$cell4_sti, probs=qnt.low, na.rm = TRUE, names = FALSE) / cell_sum_low)*100
  df$cell4[i] <- (quantile(mn$cell4_sti, probs=0.50, na.rm = TRUE, names = FALSE) / cell_sum)*100
  df$cell4.high[i] <- (quantile(mn$cell4_sti, probs=qnt.high, na.rm = TRUE, names = FALSE) / cell_sum_high)*100
  
  #PAF
  paf.low[i] <- (quantile(mn$ir100, probs = qnt.low, na.rm = TRUE, names = FALSE) - df$ir100.low[1]) / quantile(mn$ir100, probs = qnt.low, na.rm = TRUE, names = FALSE)
  paf.low[i] <- (quantile(mn$ir100, probs = 0.50, na.rm = TRUE, names = FALSE) - df$ir100[1]) / quantile(mn$ir100, probs = 0.50, na.rm = TRUE, names = FALSE)
  paf.high[i] <- (quantile(mn$ir100, probs = qnt.high, na.rm = TRUE, names = FALSE) - df$ir100.high[1]) / quantile(mn$ir100, probs = qnt.high, na.rm = TRUE, names = FALSE)
  
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
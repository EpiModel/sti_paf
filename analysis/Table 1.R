## STI PAF Table 1
# Sensitivity analysis with varying relative risks for each STI

rm(list = ls())
library("EpiModelHIV")
library("EpiModelHPC")
library("dplyr")
source("analysis/fx.R")



#Final year of follow-up - highest prob prevails
#Obtain incidence rate and simulation interval for each of the 30 combinations in Table 1
load("data/followup/sim.n6002.rda")
sim.base <- sim
epi_stats(sim.base, at = 520, qnt.low = 0.25, qnt.high = 0.75)

haz.base <- as.numeric(colMeans(tail(sim.base$epi$ir100, 52)))
ir.base <- unname(colMeans(sim.base$epi$ir100))
ir.base
incid.base <- unname(colSums(sim.base$epi$incid))

sims <- c(6002,
          6000,
          6003,
          6004,
          6005,
          6006,
          6007,
          6008,
          6009,
          6010)

qnt.low <- 0.25
qnt.high <- 0.75

scenario_num <- rep(NA, length(sims))

hiv.incid <- rep(NA, length(sims))
hiv.incid.low <- rep(NA, length(sims))
hiv.incid.high <- rep(NA, length(sims))

ir100 <- rep(NA, length(sims))
ir100.low <- rep(NA, length(sims))
ir100.high <- rep(NA, length(sims))
ir100.mean <- rep(NA, length(sims))

cell1_rectureth.low <- rep(NA, length(sims))
cell1_rectureth <- rep(NA, length(sims))
cell1_rectureth.high <- rep(NA, length(sims))

cell2_rectureth.low <- rep(NA, length(sims))
cell2_rectureth <- rep(NA, length(sims))
cell2_rectureth.high <- rep(NA, length(sims))

cell3_rectureth.low <- rep(NA, length(sims))
cell3_rectureth <- rep(NA, length(sims))
cell3_rectureth.high <- rep(NA, length(sims))

cell4_rectureth.low <- rep(NA, length(sims))
cell4_rectureth <- rep(NA, length(sims))
cell4_rectureth.high <- rep(NA, length(sims))

cell1_gcct_newinf.low <- rep(NA, length(sims))
cell1_gcct_newinf <- rep(NA, length(sims))
cell1_gcct_newinf.high <- rep(NA, length(sims))

cell2_gcct_newinf.low <- rep(NA, length(sims))
cell2_gcct_newinf <- rep(NA, length(sims))
cell2_gcct_newinf.high <- rep(NA, length(sims))

cell3_gcct_newinf.low <- rep(NA, length(sims))
cell3_gcct_newinf <- rep(NA, length(sims))
cell3_gcct_newinf.high <- rep(NA, length(sims))

cell4_gcct_newinf.low <- rep(NA, length(sims))
cell4_gcct_newinf <- rep(NA, length(sims))
cell4_gcct_newinf.high <- rep(NA, length(sims))

paf.low <- rep(NA, length(sims))
paf <- rep(NA, length(sims))
paf.high <- rep(NA, length(sims))

#add sims to dataframe as object
df <- data.frame(scenario_num, 
                 hiv.incid.low, hiv.incid, hiv.incid.high, 
                 ir100.low, ir100, ir100.high, ir100.mean,
                 cell1_rectureth.low, cell1_rectureth, cell1_rectureth.high, 
                 cell2_rectureth.low, cell2_rectureth, cell2_rectureth.high,
                 cell3_rectureth.low, cell3_rectureth, cell3_rectureth.high, 
                 cell4_rectureth.low, cell4_rectureth, cell4_rectureth.high,
                 cell1_gcct_newinf.low, cell1_gcct_newinf, cell1_gcct_newinf.high, 
                 cell2_gcct_newinf.low, cell2_gcct_newinf, cell2_gcct_newinf.high,
                 cell3_gcct_newinf.low, cell3_gcct_newinf, cell3_gcct_newinf.high, 
                 cell4_gcct_newinf.low, cell4_gcct_newinf, cell4_gcct_newinf.high,
                 paf.low, paf, paf.high)
df
for (i in seq_along(sims)) {
  
  fn <- list.files("data/followup/", pattern = as.character(sims[i]), full.names = TRUE)
  load(fn)
  
  #sim <- truncate_sim(sim, at = 2600)
  mn <- as.data.frame(sim)
  
  df$scenario_num[i] <- sims[i]
  
  #HIV diagnosis counts
  df$hiv.incid.low[i] <- quantile(mn$incid, probs = qnt.low, na.rm = TRUE, names = FALSE)
  df$hiv.incid[i] <- quantile(mn$incid, probs = 0.50, na.rm = TRUE, names = FALSE)
  df$hiv.incid.high[i] <- quantile(mn$incid, probs = qnt.high, na.rm = TRUE, names = FALSE)
  
  #Incidence Rate
  df$ir100.low[i] <- quantile(mn$ir100, probs=qnt.low, na.rm = TRUE, names = FALSE)
  df$ir100[i] <- quantile(mn$ir100, probs=0.50, na.rm = TRUE, names = FALSE)
  df$ir100.high[i] <- quantile(mn$ir100, probs=qnt.high, na.rm = TRUE, names = FALSE)
  df$ir100.mean[i] <- mean(mn$ir100)
  
  #sum cell1 - cell4 variables for proportion
  #rectureth
  cell_rectureth_sum_low <- quantile(mn$cell1_rectureth, probs = qnt.low, na.rm = TRUE, names = FALSE) +
    quantile(mn$cell2_rectureth, probs = qnt.low, na.rm = TRUE, names = FALSE) +
    quantile(mn$cell3_rectureth, probs = qnt.low, na.rm = TRUE, names = FALSE) +
    quantile(mn$cell4_rectureth, probs = qnt.low, na.rm = TRUE, names = FALSE) 
  
  cell_rectureth_sum <- quantile(mn$cell1_rectureth, probs = 0.50, na.rm = TRUE, names = FALSE) +
    quantile(mn$cell2_rectureth, probs = 0.50, na.rm = TRUE, names = FALSE) +
    quantile(mn$cell3_rectureth, probs = 0.50, na.rm = TRUE, names = FALSE) +
    quantile(mn$cell4_rectureth, probs = 0.50, na.rm = TRUE, names = FALSE) 
  
  cell_rectureth_sum_high <- quantile(mn$cell1_rectureth, probs = qnt.high, na.rm = TRUE, names = FALSE) +
    quantile(mn$cell2_rectureth, probs = qnt.high, na.rm = TRUE, names = FALSE) +
    quantile(mn$cell3_rectureth, probs = qnt.high, na.rm = TRUE, names = FALSE) +
    quantile(mn$cell4_rectureth, probs = qnt.high, na.rm = TRUE, names = FALSE)
  
  #gcct_newinf
  cell_gcct_newinf_sum_low <- quantile(mn$cell1_gcct_newinf, probs = qnt.low, na.rm = TRUE, names = FALSE) +
    quantile(mn$cell2_gcct_newinf, probs = qnt.low, na.rm = TRUE, names = FALSE) +
    quantile(mn$cell3_gcct_newinf, probs = qnt.low, na.rm = TRUE, names = FALSE) +
    quantile(mn$cell4_gcct_newinf, probs = qnt.low, na.rm = TRUE, names = FALSE) 
  
  cell_gcct_newinf_sum <- quantile(mn$cell1_gcct_newinf, probs = 0.50, na.rm = TRUE, names = FALSE) +
    quantile(mn$cell2_gcct_newinf, probs = 0.50, na.rm = TRUE, names = FALSE) +
    quantile(mn$cell3_gcct_newinf, probs = 0.50, na.rm = TRUE, names = FALSE) +
    quantile(mn$cell4_gcct_newinf, probs = 0.50, na.rm = TRUE, names = FALSE) 
  
  cell_gcct_newinf_sum_high <- quantile(mn$cell1_gcct_newinf, probs = qnt.high, na.rm = TRUE, names = FALSE) +
    quantile(mn$cell2_gcct_newinf, probs = qnt.high, na.rm = TRUE, names = FALSE) +
    quantile(mn$cell3_gcct_newinf, probs = qnt.high, na.rm = TRUE, names = FALSE) +
    quantile(mn$cell4_gcct_newinf, probs = qnt.high, na.rm = TRUE, names = FALSE)  
  
  #Cell1 - Cell4
  #rectureth
  df$cell1_rectureth.low[i] <- (quantile(mn$cell1_rectureth, probs=qnt.low, na.rm = TRUE, names = FALSE) / cell_rectureth_sum_low)*100
  df$cell1_rectureth[i] <- (quantile(mn$cell1_rectureth, probs=0.50, na.rm = TRUE, names = FALSE) / cell_rectureth_sum)*100
  df$cell1_rectureth.high[i] <- (quantile(mn$cell1_rectureth, probs=qnt.high, na.rm = TRUE, names = FALSE) / cell_rectureth_sum_high)*100
  
  df$cell2_rectureth.low[i] <- (quantile(mn$cell2_rectureth, probs=qnt.low, na.rm = TRUE, names = FALSE) / cell_rectureth_sum_low)*100
  df$cell2_rectureth[i] <- (quantile(mn$cell2_rectureth, probs=0.50, na.rm = TRUE, names = FALSE) / cell_rectureth_sum)*100
  df$cell2_rectureth.high[i] <- (quantile(mn$cell2_rectureth, probs=qnt.high, na.rm = TRUE, names = FALSE) / cell_rectureth_sum_high)*100
  
  df$cell3_rectureth.low[i] <- (quantile(mn$cell3_rectureth, probs=qnt.low, na.rm = TRUE, names = FALSE) / cell_rectureth_sum_low)*100
  df$cell3_rectureth[i] <- (quantile(mn$cell3_rectureth, probs=0.50, na.rm = TRUE, names = FALSE) / cell_rectureth_sum)*100
  df$cell3_rectureth.high[i] <- (quantile(mn$cell3_rectureth, probs=qnt.high, na.rm = TRUE, names = FALSE) / cell_rectureth_sum_high)*100
  
  df$cell4_rectureth.low[i] <- (quantile(mn$cell4_rectureth, probs=qnt.low, na.rm = TRUE, names = FALSE) / cell_rectureth_sum_low)*100
  df$cell4_rectureth[i] <- (quantile(mn$cell4_rectureth, probs=0.50, na.rm = TRUE, names = FALSE) / cell_rectureth_sum)*100
  df$cell4_rectureth.high[i] <- (quantile(mn$cell4_rectureth, probs=qnt.high, na.rm = TRUE, names = FALSE) / cell_rectureth_sum_high)*100
  
  #gcct_newinf
  #rectureth
  df$cell1_gcct_newinf.low[i] <- (quantile(mn$cell1_gcct_newinf, probs=qnt.low, na.rm = TRUE, names = FALSE) / cell_gcct_newinf_sum_low)*100
  df$cell1_gcct_newinf[i] <- (quantile(mn$cell1_gcct_newinf, probs=0.50, na.rm = TRUE, names = FALSE) / cell_gcct_newinf_sum)*100
  df$cell1_gcct_newinf.high[i] <- (quantile(mn$cell1_gcct_newinf, probs=qnt.high, na.rm = TRUE, names = FALSE) / cell_gcct_newinf_sum_high)*100
  
  df$cell2_gcct_newinf.low[i] <- (quantile(mn$cell2_gcct_newinf, probs=qnt.low, na.rm = TRUE, names = FALSE) / cell_gcct_newinf_sum_low)*100
  df$cell2_gcct_newinf[i] <- (quantile(mn$cell2_gcct_newinf, probs=0.50, na.rm = TRUE, names = FALSE) / cell_gcct_newinf_sum)*100
  df$cell2_gcct_newinf.high[i] <- (quantile(mn$cell2_gcct_newinf, probs=qnt.high, na.rm = TRUE, names = FALSE) / cell_gcct_newinf_sum_high)*100
  
  df$cell3_gcct_newinf.low[i] <- (quantile(mn$cell3_gcct_newinf, probs=qnt.low, na.rm = TRUE, names = FALSE) / cell_gcct_newinf_sum_low)*100
  df$cell3_gcct_newinf[i] <- (quantile(mn$cell3_gcct_newinf, probs=0.50, na.rm = TRUE, names = FALSE) / cell_gcct_newinf_sum)*100
  df$cell3_gcct_newinf.high[i] <- (quantile(mn$cell3_gcct_newinf, probs=qnt.high, na.rm = TRUE, names = FALSE) / cell_gcct_newinf_sum_high)*100
  
  df$cell4_gcct_newinf.low[i] <- (quantile(mn$cell4_gcct_newinf, probs=qnt.low, na.rm = TRUE, names = FALSE) / cell_gcct_newinf_sum_low)*100
  df$cell4_gcct_newinf[i] <- (quantile(mn$cell4_gcct_newinf, probs=0.50, na.rm = TRUE, names = FALSE) / cell_gcct_newinf_sum)*100
  df$cell4_gcct_newinf.high[i] <- (quantile(mn$cell4_gcct_newinf, probs=qnt.high, na.rm = TRUE, names = FALSE) / cell_gcct_newinf_sum_high)*100
  
  #PAF
  df$paf.low[i] <- (quantile(mn$ir100, probs = qnt.low, na.rm = TRUE, names = FALSE) - df$ir100[1]) / quantile(mn$ir100, probs = qnt.low, na.rm = TRUE, names = FALSE)
  df$paf[i] <- (quantile(mn$ir100, probs = 0.50, na.rm = TRUE, names = FALSE) - df$ir100[1]) / quantile(mn$ir100, probs = 0.50, na.rm = TRUE, names = FALSE)
  df$paf.high[i] <- (quantile(mn$ir100, probs = qnt.high, na.rm = TRUE, names = FALSE) - df$ir100[1]) / quantile(mn$ir100, probs = qnt.high, na.rm = TRUE, names = FALSE)
  
  cat("*")
  
}

df
write.csv(df, "C:/Users/jsjone2/Desktop/Table 1 - Final Year - Highest Prob.csv")

#Final year of follow-up - multiplicative
#Obtain incidence rate and simulation interval for each of the 30 combinations in Table 1
load("data/followup/sim.n6001.rda")
sim.base <- sim
epi_stats(sim.base, at = 520, qnt.low = 0.25, qnt.high = 0.75)

haz.base <- as.numeric(colMeans(tail(sim.base$epi$ir100, 52)))
ir.base <- unname(colMeans(sim.base$epi$ir100))
ir.base
incid.base <- unname(colSums(sim.base$epi$incid))

sims <- c(6011,
          6001,
          6012,
          6013,
          6014,
          6015,
          6016,
          6017,
          6018,
          6019)

qnt.low <- 0.25
qnt.high <- 0.75

scenario_num <- rep(NA, length(sims))

hiv.incid <- rep(NA, length(sims))
hiv.incid.low <- rep(NA, length(sims))
hiv.incid.high <- rep(NA, length(sims))

ir100 <- rep(NA, length(sims))
ir100.low <- rep(NA, length(sims))
ir100.high <- rep(NA, length(sims))
ir100.mean <- rep(NA, length(sims))

cell1_rectureth.low <- rep(NA, length(sims))
cell1_rectureth <- rep(NA, length(sims))
cell1_rectureth.high <- rep(NA, length(sims))

cell2_rectureth.low <- rep(NA, length(sims))
cell2_rectureth <- rep(NA, length(sims))
cell2_rectureth.high <- rep(NA, length(sims))

cell3_rectureth.low <- rep(NA, length(sims))
cell3_rectureth <- rep(NA, length(sims))
cell3_rectureth.high <- rep(NA, length(sims))

cell4_rectureth.low <- rep(NA, length(sims))
cell4_rectureth <- rep(NA, length(sims))
cell4_rectureth.high <- rep(NA, length(sims))

cell1_gcct_newinf.low <- rep(NA, length(sims))
cell1_gcct_newinf <- rep(NA, length(sims))
cell1_gcct_newinf.high <- rep(NA, length(sims))

cell2_gcct_newinf.low <- rep(NA, length(sims))
cell2_gcct_newinf <- rep(NA, length(sims))
cell2_gcct_newinf.high <- rep(NA, length(sims))

cell3_gcct_newinf.low <- rep(NA, length(sims))
cell3_gcct_newinf <- rep(NA, length(sims))
cell3_gcct_newinf.high <- rep(NA, length(sims))

cell4_gcct_newinf.low <- rep(NA, length(sims))
cell4_gcct_newinf <- rep(NA, length(sims))
cell4_gcct_newinf.high <- rep(NA, length(sims))

paf.low <- rep(NA, length(sims))
paf <- rep(NA, length(sims))
paf.high <- rep(NA, length(sims))

#add sims to dataframe as object
df <- data.frame(scenario_num, 
                 hiv.incid.low, hiv.incid, hiv.incid.high, 
                 ir100.low, ir100, ir100.high, ir100.mean,
                 cell1_rectureth.low, cell1_rectureth, cell1_rectureth.high, 
                 cell2_rectureth.low, cell2_rectureth, cell2_rectureth.high,
                 cell3_rectureth.low, cell3_rectureth, cell3_rectureth.high, 
                 cell4_rectureth.low, cell4_rectureth, cell4_rectureth.high,
                 cell1_gcct_newinf.low, cell1_gcct_newinf, cell1_gcct_newinf.high, 
                 cell2_gcct_newinf.low, cell2_gcct_newinf, cell2_gcct_newinf.high,
                 cell3_gcct_newinf.low, cell3_gcct_newinf, cell3_gcct_newinf.high, 
                 cell4_gcct_newinf.low, cell4_gcct_newinf, cell4_gcct_newinf.high,
                 paf.low, paf, paf.high)
df
for (i in seq_along(sims)) {
  
  fn <- list.files("data/followup/", pattern = as.character(sims[i]), full.names = TRUE)
  load(fn)
  
  #sim <- truncate_sim(sim, at = 2600)
  mn <- as.data.frame(sim)
  
  df$scenario_num[i] <- sims[i]
  
  #HIV diagnosis counts
  df$hiv.incid.low[i] <- quantile(mn$incid, probs = qnt.low, na.rm = TRUE, names = FALSE)
  df$hiv.incid[i] <- quantile(mn$incid, probs = 0.50, na.rm = TRUE, names = FALSE)
  df$hiv.incid.high[i] <- quantile(mn$incid, probs = qnt.high, na.rm = TRUE, names = FALSE)
  
  #Incidence Rate
  df$ir100.low[i] <- quantile(mn$ir100, probs=qnt.low, na.rm = TRUE, names = FALSE)
  df$ir100[i] <- quantile(mn$ir100, probs=0.50, na.rm = TRUE, names = FALSE)
  df$ir100.high[i] <- quantile(mn$ir100, probs=qnt.high, na.rm = TRUE, names = FALSE)
  df$ir100.mean[i] <- mean(mn$ir100)
  
  #sum cell1 - cell4 variables for proportion
  #rectureth
  cell_rectureth_sum_low <- quantile(mn$cell1_rectureth, probs = qnt.low, na.rm = TRUE, names = FALSE) +
    quantile(mn$cell2_rectureth, probs = qnt.low, na.rm = TRUE, names = FALSE) +
    quantile(mn$cell3_rectureth, probs = qnt.low, na.rm = TRUE, names = FALSE) +
    quantile(mn$cell4_rectureth, probs = qnt.low, na.rm = TRUE, names = FALSE) 
  
  cell_rectureth_sum <- quantile(mn$cell1_rectureth, probs = 0.50, na.rm = TRUE, names = FALSE) +
    quantile(mn$cell2_rectureth, probs = 0.50, na.rm = TRUE, names = FALSE) +
    quantile(mn$cell3_rectureth, probs = 0.50, na.rm = TRUE, names = FALSE) +
    quantile(mn$cell4_rectureth, probs = 0.50, na.rm = TRUE, names = FALSE) 
  
  cell_rectureth_sum_high <- quantile(mn$cell1_rectureth, probs = qnt.high, na.rm = TRUE, names = FALSE) +
    quantile(mn$cell2_rectureth, probs = qnt.high, na.rm = TRUE, names = FALSE) +
    quantile(mn$cell3_rectureth, probs = qnt.high, na.rm = TRUE, names = FALSE) +
    quantile(mn$cell4_rectureth, probs = qnt.high, na.rm = TRUE, names = FALSE)
  
  #gcct_newinf
  cell_gcct_newinf_sum_low <- quantile(mn$cell1_gcct_newinf, probs = qnt.low, na.rm = TRUE, names = FALSE) +
    quantile(mn$cell2_gcct_newinf, probs = qnt.low, na.rm = TRUE, names = FALSE) +
    quantile(mn$cell3_gcct_newinf, probs = qnt.low, na.rm = TRUE, names = FALSE) +
    quantile(mn$cell4_gcct_newinf, probs = qnt.low, na.rm = TRUE, names = FALSE) 
  
  cell_gcct_newinf_sum <- quantile(mn$cell1_gcct_newinf, probs = 0.50, na.rm = TRUE, names = FALSE) +
    quantile(mn$cell2_gcct_newinf, probs = 0.50, na.rm = TRUE, names = FALSE) +
    quantile(mn$cell3_gcct_newinf, probs = 0.50, na.rm = TRUE, names = FALSE) +
    quantile(mn$cell4_gcct_newinf, probs = 0.50, na.rm = TRUE, names = FALSE) 
  
  cell_gcct_newinf_sum_high <- quantile(mn$cell1_gcct_newinf, probs = qnt.high, na.rm = TRUE, names = FALSE) +
    quantile(mn$cell2_gcct_newinf, probs = qnt.high, na.rm = TRUE, names = FALSE) +
    quantile(mn$cell3_gcct_newinf, probs = qnt.high, na.rm = TRUE, names = FALSE) +
    quantile(mn$cell4_gcct_newinf, probs = qnt.high, na.rm = TRUE, names = FALSE)  
  
  #Cell1 - Cell4
  #rectureth
  df$cell1_rectureth.low[i] <- (quantile(mn$cell1_rectureth, probs=qnt.low, na.rm = TRUE, names = FALSE) / cell_rectureth_sum_low)*100
  df$cell1_rectureth[i] <- (quantile(mn$cell1_rectureth, probs=0.50, na.rm = TRUE, names = FALSE) / cell_rectureth_sum)*100
  df$cell1_rectureth.high[i] <- (quantile(mn$cell1_rectureth, probs=qnt.high, na.rm = TRUE, names = FALSE) / cell_rectureth_sum_high)*100
  
  df$cell2_rectureth.low[i] <- (quantile(mn$cell2_rectureth, probs=qnt.low, na.rm = TRUE, names = FALSE) / cell_rectureth_sum_low)*100
  df$cell2_rectureth[i] <- (quantile(mn$cell2_rectureth, probs=0.50, na.rm = TRUE, names = FALSE) / cell_rectureth_sum)*100
  df$cell2_rectureth.high[i] <- (quantile(mn$cell2_rectureth, probs=qnt.high, na.rm = TRUE, names = FALSE) / cell_rectureth_sum_high)*100
  
  df$cell3_rectureth.low[i] <- (quantile(mn$cell3_rectureth, probs=qnt.low, na.rm = TRUE, names = FALSE) / cell_rectureth_sum_low)*100
  df$cell3_rectureth[i] <- (quantile(mn$cell3_rectureth, probs=0.50, na.rm = TRUE, names = FALSE) / cell_rectureth_sum)*100
  df$cell3_rectureth.high[i] <- (quantile(mn$cell3_rectureth, probs=qnt.high, na.rm = TRUE, names = FALSE) / cell_rectureth_sum_high)*100
  
  df$cell4_rectureth.low[i] <- (quantile(mn$cell4_rectureth, probs=qnt.low, na.rm = TRUE, names = FALSE) / cell_rectureth_sum_low)*100
  df$cell4_rectureth[i] <- (quantile(mn$cell4_rectureth, probs=0.50, na.rm = TRUE, names = FALSE) / cell_rectureth_sum)*100
  df$cell4_rectureth.high[i] <- (quantile(mn$cell4_rectureth, probs=qnt.high, na.rm = TRUE, names = FALSE) / cell_rectureth_sum_high)*100
  
  #gcct_newinf
  #rectureth
  df$cell1_gcct_newinf.low[i] <- (quantile(mn$cell1_gcct_newinf, probs=qnt.low, na.rm = TRUE, names = FALSE) / cell_gcct_newinf_sum_low)*100
  df$cell1_gcct_newinf[i] <- (quantile(mn$cell1_gcct_newinf, probs=0.50, na.rm = TRUE, names = FALSE) / cell_gcct_newinf_sum)*100
  df$cell1_gcct_newinf.high[i] <- (quantile(mn$cell1_gcct_newinf, probs=qnt.high, na.rm = TRUE, names = FALSE) / cell_gcct_newinf_sum_high)*100
  
  df$cell2_gcct_newinf.low[i] <- (quantile(mn$cell2_gcct_newinf, probs=qnt.low, na.rm = TRUE, names = FALSE) / cell_gcct_newinf_sum_low)*100
  df$cell2_gcct_newinf[i] <- (quantile(mn$cell2_gcct_newinf, probs=0.50, na.rm = TRUE, names = FALSE) / cell_gcct_newinf_sum)*100
  df$cell2_gcct_newinf.high[i] <- (quantile(mn$cell2_gcct_newinf, probs=qnt.high, na.rm = TRUE, names = FALSE) / cell_gcct_newinf_sum_high)*100
  
  df$cell3_gcct_newinf.low[i] <- (quantile(mn$cell3_gcct_newinf, probs=qnt.low, na.rm = TRUE, names = FALSE) / cell_gcct_newinf_sum_low)*100
  df$cell3_gcct_newinf[i] <- (quantile(mn$cell3_gcct_newinf, probs=0.50, na.rm = TRUE, names = FALSE) / cell_gcct_newinf_sum)*100
  df$cell3_gcct_newinf.high[i] <- (quantile(mn$cell3_gcct_newinf, probs=qnt.high, na.rm = TRUE, names = FALSE) / cell_gcct_newinf_sum_high)*100
  
  df$cell4_gcct_newinf.low[i] <- (quantile(mn$cell4_gcct_newinf, probs=qnt.low, na.rm = TRUE, names = FALSE) / cell_gcct_newinf_sum_low)*100
  df$cell4_gcct_newinf[i] <- (quantile(mn$cell4_gcct_newinf, probs=0.50, na.rm = TRUE, names = FALSE) / cell_gcct_newinf_sum)*100
  df$cell4_gcct_newinf.high[i] <- (quantile(mn$cell4_gcct_newinf, probs=qnt.high, na.rm = TRUE, names = FALSE) / cell_gcct_newinf_sum_high)*100
  
  #PAF
  df$paf.low[i] <- (quantile(mn$ir100, probs = qnt.low, na.rm = TRUE, names = FALSE) - df$ir100[1]) / quantile(mn$ir100, probs = qnt.low, na.rm = TRUE, names = FALSE)
  df$paf[i] <- (quantile(mn$ir100, probs = 0.50, na.rm = TRUE, names = FALSE) - df$ir100[1]) / quantile(mn$ir100, probs = 0.50, na.rm = TRUE, names = FALSE)
  df$paf.high[i] <- (quantile(mn$ir100, probs = qnt.high, na.rm = TRUE, names = FALSE) - df$ir100[1]) / quantile(mn$ir100, probs = qnt.high, na.rm = TRUE, names = FALSE)
  
  cat("*")
  
}

df
write.csv(df, "C:/Users/jsjone2/Desktop/Table 1 - Final Year - Multiplicative.csv")



###OLD CODE WITH SYPHILIS SIMS#####

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
ir.base <- unname(colMeans(sim.base$epi$ir100))
ir.base
incid.base <- unname(colSums(sim.base$epi$incid))

sims <- c(5001,
          5004,
          5007,
          5010,
          5013,
          5002,
          5005,
          5008,
          5011,
          5014,
          5003,
          5006,
          5009,
          5012,
          5015)

qnt.low <- 0.25
qnt.high <- 0.75

scenario_num <- rep(NA, length(sims))

hiv.incid <- rep(NA, length(sims))
hiv.incid.low <- rep(NA, length(sims))
hiv.incid.high <- rep(NA, length(sims))

ir100 <- rep(NA, length(sims))
ir100.low <- rep(NA, length(sims))
ir100.high <- rep(NA, length(sims))
ir100.mean <- rep(NA, length(sims))

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
                 ir100.low, ir100, ir100.high, ir100.mean,
                 cell1.low, cell1, cell1.high, 
                 cell2.low, cell2, cell2.high,
                 cell3.low, cell3, cell3.high, 
                 cell4.low, cell4, cell4.high,
                 paf.low, paf, paf.high)
df
for (i in seq_along(sims)) {
  
  fn <- list.files("data/followup/", pattern = as.character(sims[i]), full.names = TRUE)
  load(fn)
  
  #sim <- truncate_sim(sim, at = 2600)
  mn <- as.data.frame(sim)
  
  df$scenario_num[i] <- sims[i]
  
  #HIV diagnosis counts
  df$hiv.incid.low[i] <- quantile(mn$incid, probs = qnt.low, na.rm = TRUE, names = FALSE)
  df$hiv.incid[i] <- quantile(mn$incid, probs = 0.50, na.rm = TRUE, names = FALSE)
  df$hiv.incid.high[i] <- quantile(mn$incid, probs = qnt.high, na.rm = TRUE, names = FALSE)
  
  #Incidence Rate
  df$ir100.low[i] <- quantile(mn$ir100, probs=qnt.low, na.rm = TRUE, names = FALSE)
  df$ir100[i] <- quantile(mn$ir100, probs=0.50, na.rm = TRUE, names = FALSE)
  df$ir100.high[i] <- quantile(mn$ir100, probs=qnt.high, na.rm = TRUE, names = FALSE)
  df$ir100.mean[i] <- mean(mn$ir100)
  
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
  df$paf.low[i] <- (quantile(mn$ir100, probs = qnt.low, na.rm = TRUE, names = FALSE) - df$ir100[1]) / quantile(mn$ir100, probs = qnt.low, na.rm = TRUE, names = FALSE)
  df$paf[i] <- (quantile(mn$ir100, probs = 0.50, na.rm = TRUE, names = FALSE) - df$ir100[1]) / quantile(mn$ir100, probs = 0.50, na.rm = TRUE, names = FALSE)
  df$paf.high[i] <- (quantile(mn$ir100, probs = qnt.high, na.rm = TRUE, names = FALSE) - df$ir100[1]) / quantile(mn$ir100, probs = qnt.high, na.rm = TRUE, names = FALSE)
  
  cat("*")
  
}

df
write.csv(df, "C:/Users/jsjone2/Desktop/Table 1 - Final Year - Highest Prob.csv")

#Final year of follow-up - multiplicative
#Obtain incidence rate and simulation interval for each of the 30 combinations in Table 1
load("data/followup/sim.n5016.rda")
sim.base <- sim
epi_stats(sim.base, at = 520, qnt.low = 0.25, qnt.high = 0.75)

haz.base <- as.numeric(colMeans(tail(sim.base$epi$ir100, 52)))
ir.base <- unname(colMeans(sim.base$epi$ir100)) * 1000
incid.base <- unname(colSums(sim.base$epi$incid))

sims <- c(5016,
          5019,
          5022,
          5025,
          5028,
          5017,
          5020,
          5023,
          5026,
          5029,
          5018,
          5021,
          5024,
          5027,
          5030)

qnt.low <- 0.25
qnt.high <- 0.75

scenario_num <- rep(NA, length(sims))

hiv.incid <- rep(NA, length(sims))
hiv.incid.low <- rep(NA, length(sims))
hiv.incid.high <- rep(NA, length(sims))

ir100 <- rep(NA, length(sims))
ir100.low <- rep(NA, length(sims))
ir100.high <- rep(NA, length(sims))
ir100.mean <- rep(NA, length(sims))

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
                 ir100.low, ir100, ir100.high, ir100.mean,
                 cell1.low, cell1, cell1.high, 
                 cell2.low, cell2, cell2.high,
                 cell3.low, cell3, cell3.high, 
                 cell4.low, cell4, cell4.high,
                 paf.low, paf, paf.high)
df
for (i in seq_along(sims)) {
  
  fn <- list.files("data/followup/", pattern = as.character(sims[i]), full.names = TRUE)
  load(fn)
  
  #sim <- truncate_sim(sim, at = 2600)
  mn <- as.data.frame(sim)
  
  df$scenario_num[i] <- sims[i]
  
  #HIV diagnosis counts
  df$hiv.incid.low[i] <- quantile(mn$incid, probs = qnt.low, na.rm = TRUE, names = FALSE)
  df$hiv.incid[i] <- quantile(mn$incid, probs = 0.50, na.rm = TRUE, names = FALSE)
  df$hiv.incid.high[i] <- quantile(mn$incid, probs = qnt.high, na.rm = TRUE, names = FALSE)
  
  #Incidence Rate
  df$ir100.low[i] <- quantile(mn$ir100, probs=qnt.low, na.rm = TRUE, names = FALSE)
  df$ir100[i] <- quantile(mn$ir100, probs=0.50, na.rm = TRUE, names = FALSE)
  df$ir100.high[i] <- quantile(mn$ir100, probs=qnt.high, na.rm = TRUE, names = FALSE)
  df$ir100.mean[i] <- mean(mn$ir100)
  
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
  df$paf.low[i] <- (quantile(mn$ir100, probs = qnt.low, na.rm = TRUE, names = FALSE) - df$ir100[1]) / quantile(mn$ir100, probs = qnt.low, na.rm = TRUE, names = FALSE)
  df$paf[i] <- (quantile(mn$ir100, probs = 0.50, na.rm = TRUE, names = FALSE) - df$ir100[1]) / quantile(mn$ir100, probs = 0.50, na.rm = TRUE, names = FALSE)
  df$paf.high[i] <- (quantile(mn$ir100, probs = qnt.high, na.rm = TRUE, names = FALSE) - df$ir100[1]) / quantile(mn$ir100, probs = qnt.high, na.rm = TRUE, names = FALSE)
  
  cat("*")
  
}

df

write.csv(df, "C:/Users/jsjone2/Desktop/Table 1 - Final Year - Multiplicative.csv")
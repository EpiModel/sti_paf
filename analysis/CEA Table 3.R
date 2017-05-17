## STI PrEP CEA Table 3
# Varying coverage of annual and high-risk testing

rm(list = ls())
library("EpiModelHIV")
library("EpiModelHPC")
library("dplyr")
source("analysis/fx.R")

## Base STI lower-risk testing interval (364 days): n3054
## Varying STI lower-risk testing interval: n3131-n3141
## Base STI higher-risk testing interval: n3014
## Varying STI higher-risk testing interval: n3153 - n3173

# Lower and Higher-Risk
sims <- c(3000, 3131:3141, 3054, 3153, 3156, 3158, 3160, 3162, 3164, 3166, 3168, 3170, 3172, 3014)

qnt.low <- 0.25
qnt.high <- 0.75

annint <- rep(NA, length(sims))
hrint <- rep(NA, length(sims))
anncov <- rep(NA, length(sims))
hrcov <- rep(NA, length(sims))
QALY <- rep(NA, length(sims))
hivtestcosts <- rep(NA, length(sims))
hivhealthcosts <- rep(NA, length(sims))
stisympttestcosts <- rep(NA, length(sims))
gc.asympttestcosts <- rep(NA, length(sims))
ct.asympttestcosts <- rep(NA, length(sims))
syph.asympttestcosts <- rep(NA, length(sims))
rect.asympttestcosts <- rep(NA, length(sims))
ureth.asympttestcosts <- rep(NA, length(sims))
stiasympttestcosts <- rep(NA, length(sims))
stiasympttests <- rep(NA, length(sims))
stitxcosts <- rep(NA, length(sims))
overallcosts <- rep(NA, length(sims))

incremcost <- rep(NA, length(sims))
incremover <- rep(NA, length(sims))
#incrembasescreen <- rep(NA, length(sims))

incremQALY <- rep(NA, length(sims))
incremQALYover <- rep(NA, length(sims))
#incremQALYbasescreen <- rep(NA, length(sims))

costpqaly <- rep(NA, length(sims))
costpincremqaly <- rep(NA, length(sims))
costpoverqaly <- rep(NA, length(sims))
#costpoverqalybasescreen <- rep(NA, length(sims))

# add sims to data frame as an object?
df <- data.frame(annint, hrint, anncov, hrcov, QALY, stiasympttests, hivtestcosts, hivhealthcosts, stisympttestcosts,
                 gc.asympttestcosts, ct.asympttestcosts, syph.asympttestcosts,
                 rect.asympttestcosts, ureth.asympttestcosts, stitxcosts,
                 stiasympttestcosts, overallcosts, incremcost, incremover, 
                 #incrembasescreen,
                 incremQALY, incremQALYover, 
                 #incremQALYbasescreen, 
                 costpqaly, costpincremqaly, costpoverqaly)
                #, costpoverqalybasescreen)

for (i in seq_along(sims)) {
    
    fn <- list.files("data/followup/", pattern = as.character(sims[i]), full.names = TRUE)
    load(fn)
    
    #sim <- truncate_sim(sim, at = 2600)
    
    df$annint[i] <- sim$param$stitest.active.int
    df$hrint[i] <- sim$param$sti.highrisktest.int
    df$anncov[i] <- sim$param$stianntest.coverage
    df$hrcov[i] <- sim$param$stihighrisktest.coverage
    
    # Time in HIV stages
    time.hivneg <- quantile(unname(colSums(sim$epi$time.hivneg) / 52), probs = 0.5)
    #time.on.prep <- quantile(unname(colSums(sim$epi$time.on.prep) / 52), probs = 0.5)
    #time.off.prep <- quantile(unname(colSums(sim$epi$time.off.prep) / 52), probs = 0.5)
    stage.time.ar.ndx <- quantile(unname(colSums(sim$epi$stage.time.ar.ndx) / 52), probs = 0.5)
    stage.time.af.ndx <- quantile(unname(colSums(sim$epi$stage.time.af.ndx) / 52), probs = 0.5)
    stage.time.chronic.ndx <- quantile(unname(colSums(sim$epi$stage.time.chronic.ndx) / 52), probs = 0.5)
    stage.time.aids.ndx <- quantile(unname(colSums(sim$epi$stage.time.aids.ndx) / 52), probs = 0.5)
    stage.time.ar.dx <- quantile(unname(colSums(sim$epi$stage.time.ar.dx) / 52), probs = 0.5)
    stage.time.af.dx <- quantile(unname(colSums(sim$epi$stage.time.af.dx) / 52), probs = 0.5)
    stage.time.chronic.dx <- quantile(unname(colSums(sim$epi$stage.time.chronic.dx) / 52), probs = 0.5)
    stage.time.aids.dx <- quantile(unname(colSums(sim$epi$stage.time.aids.dx) / 52), probs = 0.5)
    stage.time.ar.art <- quantile(unname(colSums(sim$epi$stage.time.ar.art) / 52), probs = 0.5)
    stage.time.af.art <- quantile(unname(colSums(sim$epi$stage.time.af.art) / 52), probs = 0.5)
    stage.time.chronic.art <- quantile(unname(colSums(sim$epi$stage.time.chronic.art) / 52), probs = 0.5)
    stage.time.aids.art <- quantile(unname(colSums(sim$epi$stage.time.aids.art) / 52), probs = 0.5)
    
    # Total tests
    totalhivtests.pos <- totalrGCasympttests.pos <- totaluGCasympttests.pos <- totalGCasympttests.pos <- totalrCTasympttests.pos <- totaluCTasympttests.pos <- totalCTasympttests.pos <- totalsyphasympttests.pos <- 0
    totalhivtests <- quantile(colMeans(unname(sim$epi$totalhivtests[521, ])), probs = 0.5)
    totalrGCsympttests <- quantile(colMeans(unname(sim$epi$totalrGCsympttests[521, ])), probs = 0.5)
    totaluGCsympttests <- quantile(colMeans(unname(sim$epi$totaluGCsympttests[521, ])), probs = 0.5)
    totalGCsympttests <- quantile(colMeans(unname(sim$epi$totalGCsympttests[521, ])), probs = 0.5)
    totalrCTsympttests <- quantile(colMeans(unname(sim$epi$totalrCTsympttests[521, ])), probs = 0.5)
    totaluCTsympttests <- quantile(colMeans(unname(sim$epi$totaluCTsympttests[521, ])), probs = 0.5)
    totalCTsympttests <- quantile(colMeans(unname(sim$epi$totalCTsympttests[521, ])), probs = 0.5)
    totalsyphsympttests <- quantile(colMeans(unname(sim$epi$totalsyphsympttests[521, ])), probs = 0.5)
    totalrGCasympttests <- quantile(colMeans(unname(sim$epi$totalrGCasympttests[521, ])), probs = 0.5)
    totaluGCasympttests <- quantile(colMeans(unname(sim$epi$totaluGCasympttests[521, ])), probs = 0.5)
    totalGCasympttests <- quantile(colMeans(unname(sim$epi$totalGCasympttests[521, ])), probs = 0.5)
    totalrCTasympttests <- quantile(colMeans(unname(sim$epi$totalrCTasympttests[521, ])), probs = 0.5)
    totaluCTasympttests <- quantile(colMeans(unname(sim$epi$totaluCTasympttests[521, ])), probs = 0.5)
    totalCTasympttests <- quantile(colMeans(unname(sim$epi$totalCTasympttests[521, ])), probs = 0.5)
    totalsyphasympttests <- quantile(colMeans(unname(sim$epi$totalsyphasympttests[521, ])), probs = 0.5)
    # totalhivtests.pos <- quantile(colMeans(unname(sim$epi$totalhivtests.pos[521, ])), probs = 0.5)
    # totalrGCasympttests.pos <- quantile(colMeans(unname(sim$epi$totalrGCasympttests.pos[521, ])), probs = 0.5)
    # totaluGCasympttests.pos <- quantile(colMeans(unname(sim$epi$totaluGCasympttests.pos[521, ])), probs = 0.5)
    # totalGCasympttests.pos <- quantile(colMeans(unname(sim$epi$totalGCasympttests.pos[521, ])), probs = 0.5)
    # totalrCTasympttests.pos <- quantile(colMeans(unname(sim$epi$totalrCTasympttests.pos[521, ])), probs = 0.5)
    # totaluCTasympttests.pos <- quantile(colMeans(unname(sim$epi$totaluCTasympttests.pos[521, ])), probs = 0.5)
    # totalCTasympttests.pos <- quantile(colMeans(unname(sim$epi$totalCTasympttests.pos[521, ])), probs = 0.5)
    # totalsyphasympttests.pos <- quantile(colMeans(unname(sim$epi$totalsyphasympttests.pos[521, ])), probs = 0.5)
    
    # STI Tx
    txCT <- quantile(unname(colSums(sim$epi$txCT)), probs = 0.5)
    txGC <- quantile(unname(colSums(sim$epi$txGC)), probs = 0.5)
    txsyph <- quantile(unname(colSums(sim$epi$txsyph)), probs = 0.5)
    # txearlysyph <- quantile(unname(colSums(sim$epi$txearlysyph)), probs = 0.5)
    # txlatesyph <- quantile(unname(colSums(sim$epi$txlatesyph)), probs = 0.5)
    
    # QALY
    df$QALY[i] <- sum((time.hivneg * 1), 
                      (stage.time.ar.ndx * 0.865), (stage.time.af.ndx * 0.865),
                      (stage.time.chronic.ndx * 0.865), (stage.time.aids.ndx * 0.72),
                      (stage.time.ar.dx * 0.795), (stage.time.af.dx * 0.795),
                      (stage.time.chronic.dx * 0.795), (stage.time.aids.dx * 0.72),
                      (stage.time.ar.art * 0.83), (stage.time.af.art * 0.83),
                      (stage.time.chronic.art * 0.83), (stage.time.aids.art * 0.82))

    df$stiasympttests[i] <- sum((totalrGCasympttests), (totaluGCasympttests),
                                (totalrCTasympttests), (totaluCTasympttests),
                                (totalsyphasympttests))
    df$hivtestcosts[i] <- sum(((totalhivtests - totalhivtests.pos) * 64.75), (totalhivtests.pos * 614.75))
    df$hivhealthcosts[i] <- sum((time.hivneg * 4469.81), 
                             (stage.time.ar.ndx * 4502.83), (stage.time.ar.dx * 4502.83), (stage.time.af.art * 17071.8),
                             (stage.time.af.ndx * 4502.83),  (stage.time.af.dx * 4502.83), (stage.time.af.art * 17071.8),
                             (stage.time.chronic.ndx * 10558.70), (stage.time.chronic.dx * 10558.70), (stage.time.chronic.art * 23842),
                             (stage.time.aids.ndx * 28533.69), (stage.time.aids.dx * 28533.69), (stage.time.aids.art * 27990.41))
    df$stisympttestcosts[i] <- sum((totalrGCsympttests * 45.62), (totaluGCsympttests * 45.62),
                                   (totalrCTsympttests * 45.62), (totaluCTsympttests * 45.62),
                                   (totalsyphsympttests * 33.96))
    df$gc.asympttestcosts[i] <- sum((totalGCasympttests * 45.62))
    df$ct.asympttestcosts[i] <- sum((totalCTasympttests * 45.62))
    df$syph.asympttestcosts[i] <- sum((totalsyphasympttests * 33.96))
    df$rect.asympttestcosts[i] <- sum((totalrCTasympttests * 45.62), (totalrGCasympttests * 45.62))
    df$ureth.asympttestcosts[i] <- sum((totaluCTasympttests * 45.62), (totaluGCasympttests * 45.62))
    df$stiasympttestcosts[i] <- sum(sum((totalrGCasympttests * 45.62), (totaluGCasympttests * 45.62),
                                        (totalrCTasympttests * 45.62), (totaluCTasympttests * 45.62),
                                        (totalsyphasympttests * 33.96)))
    df$stitxcosts[i] <- sum((txsyph * 99.35), (txGC * 53.35), (txCT * 53.35))
    # df$stitxcosts[i] <- sum((txlatesyph * 99.35 * 3), (txearlysyph * 99.35), (txGC * 53.35), (txCT * 53.35))
    
    df$overallcosts[i] <- sum(df$hivtestcosts[i], df$hivhealthcosts[i], df$stisympttestcosts[i], df$stiasympttestcosts[i], df$stitxcosts[i])
    
    df$incremcost[1] <- df$overallcosts[i]
    df$incremover[1] <- df$overallcosts[i]
    #df$incrembasescreen[1] <- 0
    
    df$incremQALY[1] <- df$QALY[i]
    df$incremQALYover[1] <- df$QALY[i]
    #df$incremQALYbasescreen[1] <- 0
    
    df$costpqaly[1] <- df$overallcosts[i] / df$QALY[i]
    df$costpincremqaly[1] <- 0
    df$costpoverqaly[1] <- 0
    #df$costpoverqalybasescreen[1:2] <- 0
    
    
    if (i >= 2) {
        
    df$incremcost[i] <- df$stiasympttestcosts[i] - df$stiasympttestcosts[i - 1]
    df$incremover[i] <- df$stiasympttestcosts[i] - df$stiasympttestcosts[1]
    #df$incrembasescreen[i] <- df$stiasympttestcosts[i] - df$stiasympttestcosts[2]
    
    df$incremQALY[i] <- df$QALY[i] - df$QALY[i - 1]
    df$incremQALYover[i] <- df$QALY[i] - df$QALY[1]
    #df$incremQALYbasescreen[i] <- df$QALY[i] - df$QALY[2]
    
    df$costpqaly[i] <- df$overallcosts[i] / df$QALY[i]
    df$costpincremqaly[i] <- (df$incremcost[i] / df$incremQALY[i])
    df$costpoverqaly[i] <- (df$incremover[i] / df$incremQALYover[i])
    #df$costpoverqalybasescreen[i] <- (df$incrembasescreen[i] / df$incremQALYbasescreen[i])
    
    }
    
    cat("*")
    
}
View(df)
write.csv(df, "C:/Users/kweiss2/Documents/GitHub/stitestguidelines/analysis/CEA Table 3.csv") 

# On laptop
#write.csv(df, "/Users/kvnweiss/stitestguidelines/analysis/Table 3.csv")


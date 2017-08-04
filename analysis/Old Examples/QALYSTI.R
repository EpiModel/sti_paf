library("EpiModelHIV")
library("EpiModelHPC")
library("dplyr")
source("analysis/fx.R")


sim <- truncate_sim(sim, at = 2600)

# Time in HIV stages
time.hivneg <- as.numeric(sim$epi$time.hivneg[520, ]) / 52
round(quantile(time.hivneg, probs = c(0.025, 0.5, 0.975)), 3)

time.on.prep <- as.numeric(sim$epi$time.on.prep[520, ]) / 52
round(quantile(time.on.prep, probs = c(0.025, 0.5, 0.975)), 3)

time.off.prep <- as.numeric(sim$epi$time.off.prep[520, ]) / 52
round(quantile(time.off.prep, probs = c(0.025, 0.5, 0.975)), 3)

stage.time.ar.ndx <- as.numeric(sim$epi$stage.time.ar.ndx[520, ]) / 52
round(quantile(stage.time.ar.ndx, probs = c(0.025, 0.5, 0.975)), 3)

stage.time.af.ndx <- as.numeric(sim$epi$stage.time.af.ndx[520, ]) / 52
round(quantile(stage.time.af.ndx, probs = c(0.025, 0.5, 0.975)), 3)

stage.time.chronic.ndx <- as.numeric(sim$epi$stage.time.chronic.ndx[520, ]) / 52
round(quantile(stage.time.chronic.ndx, probs = c(0.025, 0.5, 0.975)), 3)

stage.time.aids.ndx <- as.numeric(sim$epi$stage.time.aids.ndx[520,]) / 52
round(quantile(stage.time.aids.ndx, probs = c(0.025, 0.5, 0.975)), 3)

stage.time.ar.dx <- as.numeric(sim$epi$stage.time.ar.dx[520, ]) / 52
round(quantile(stage.time.ar.dx, probs = c(0.025, 0.5, 0.975)), 3)

stage.time.af.dx <- as.numeric(sim$epi$stage.time.af.dx[520, ]) / 52
round(quantile(stage.time.af.dx, probs = c(0.025, 0.5, 0.975)), 3)

stage.time.chronic.dx <- as.numeric(sim$epi$stage.time.chronic.dx[520, ]) / 52
round(quantile(stage.time.chronic.dx, probs = c(0.025, 0.5, 0.975)), 3)

stage.time.aids.dx <- as.numeric(sim$epi$stage.time.aids.dx[520,]) / 52
round(quantile(stage.time.aids.dx, probs = c(0.025, 0.5, 0.975)), 3)

stage.time.ar.art <- as.numeric(sim$epi$stage.time.ar.art[520, ]) / 52
round(quantile(stage.time.ar.art, probs = c(0.025, 0.5, 0.975)), 3)

stage.time.af.art <- as.numeric(sim$epi$stage.time.af.art[520, ]) / 52
round(quantile(stage.time.af.art, probs = c(0.025, 0.5, 0.975)), 3)

stage.time.chronic.art <- as.numeric(sim$epi$stage.time.chronic.art[520, ]) / 52
round(quantile(stage.time.chronic.art, probs = c(0.025, 0.5, 0.975)), 3)

stage.time.aids.art <- as.numeric(sim$epi$stage.time.aids.art[520,]) / 52
round(quantile(stage.time.aids.art, probs = c(0.025, 0.5, 0.975)), 3)


# Total HIV Tests (PrEP or not)
totalhivtests <- as.numeric(sim$epi$totalhivtests[520, ])
round(quantile(totalhivtests, probs = c(0.025, 0.5, 0.975)), 3)


# Total Symptomatic STI Tests
totalrGCsympttests <- as.numeric(sim$epi$totalrGCsympttests[520, ])
round(quantile(totalrGCsympttests, probs = c(0.025, 0.5, 0.975)), 3)

totaluGCsympttests <- as.numeric(sim$epi$totaluGCsympttests[520, ])
round(quantile(totaluGCsympttests, probs = c(0.025, 0.5, 0.975)), 3)

totalGCsympttests <- as.numeric(sim$epi$totalGCsympttests[520, ])
round(quantile(totalGCsympttests, probs = c(0.025, 0.5, 0.975)), 3)

totalrCTsympttests <- as.numeric(sim$epi$totalrCTsympttests[520, ])
round(quantile(totalrCTsympttests, probs = c(0.025, 0.5, 0.975)), 3)

totaluCTsympttests <- as.numeric(sim$epi$totaluCTsympttests[520, ])
round(quantile(totaluCTsympttests, probs = c(0.025, 0.5, 0.975)), 3)

totalCTsympttests <- as.numeric(sim$epi$totalCTsympttests[520, ])
round(quantile(totalCTsympttests, probs = c(0.025, 0.5, 0.975)), 3)

totalsyphsympttests <- as.numeric(sim$epi$totalsyphsympttests[520, ])
round(quantile(totalsyphsympttests, probs = c(0.025, 0.5, 0.975)), 3)


# Total Asymptomatic STI Tests
totalrGCasympttests <- as.numeric(sim$epi$totalrGCasympttests[520, ])
round(quantile(totalrGCasympttests, probs = c(0.025, 0.5, 0.975)), 3)

totaluGCasympttests <- as.numeric(sim$epi$totaluGCasympttests[520, ])
round(quantile(totaluGCasympttests, probs = c(0.025, 0.5, 0.975)), 3)

totalGCasympttests <- as.numeric(sim$epi$totalGCasympttests[520, ])
round(quantile(totalGCasympttests, probs = c(0.025, 0.5, 0.975)), 3)

totalrCTasympttests <- as.numeric(sim$epi$totalrCTasympttests[520, ])
round(quantile(totalrCTasympttests, probs = c(0.025, 0.5, 0.975)), 3)

totaluCTasympttests <- as.numeric(sim$epi$totaluCTasympttests[520, ])
round(quantile(totaluCTasympttests, probs = c(0.025, 0.5, 0.975)), 3)

totalCTasympttests <- as.numeric(sim$epi$totalCTasympttests[520, ])
round(quantile(totalCTasympttests, probs = c(0.025, 0.5, 0.975)), 3)

totalsyphasympttests <- as.numeric(sim$epi$totalsyphasympttests[520, ])
round(quantile(totalsyphasympttests, probs = c(0.025, 0.5, 0.975)), 3)


# Total Positive HIV Tests
totalhivtests.pos <- as.numeric(sim$epi$totalhivtests.pos[520, ])
round(quantile(totalhivtests, probs = c(0.025, 0.5, 0.975)), 3)


# Total Positive Asymptomatic STI Tests
totalrGCasympttests.pos <- as.numeric(sim$epi$totalrGCasympttests.pos[520, ])
round(quantile(totalrGCasympttests, probs = c(0.025, 0.5, 0.975)), 3)

totaluGCasympttests.pos <- as.numeric(sim$epi$totaluGCasympttests.pos[520, ])
round(quantile(totaluGCasympttests, probs = c(0.025, 0.5, 0.975)), 3)

totalGCasympttests.pos <- as.numeric(sim$epi$totalGCasympttests.pos[520, ])
round(quantile(totalGCasympttests, probs = c(0.025, 0.5, 0.975)), 3)

totalrCTasympttests.pos <- as.numeric(sim$epi$totalrCTasympttests.pos[520, ])
round(quantile(totalrCTasympttests, probs = c(0.025, 0.5, 0.975)), 3)

totaluCTasympttests.pos <- as.numeric(sim$epi$totaluCTasympttests.pos[520, ])
round(quantile(totaluCTasympttests, probs = c(0.025, 0.5, 0.975)), 3)

totalCTasympttests.pos <- as.numeric(sim$epi$totalCTasympttests.pos[520, ])
round(quantile(totalCTasympttests, probs = c(0.025, 0.5, 0.975)), 3)

totalsyphasympttests.pos <- as.numeric(sim$epi$totalsyphasympttests.pos[520, ])
round(quantile(totalsyphasympttests, probs = c(0.025, 0.5, 0.975)), 3)


## Create QALY Data Frame
df <- as.data.frame(cbind(
    rbind(quantile(time.hivneg, probs = 0.25), # quantile(time.on.prep, probs = 0.25), 
          quantile(stage.time.ar.ndx, probs = 0.25),
          quantile(stage.time.af.ndx, probs = 0.25), quantile(stage.time.chronic.ndx, probs = 0.25),
          quantile(stage.time.aids.ndx, probs = 0.25),
          quantile(stage.time.ar.dx, probs = 0.25),
          quantile(stage.time.af.dx, probs = 0.25), quantile(stage.time.chronic.dx, probs = 0.25),
          quantile(stage.time.aids.dx, probs = 0.25),
          quantile(stage.time.ar.art, probs = 0.25),
          quantile(stage.time.af.art, probs = 0.25), quantile(stage.time.chronic.art, probs = 0.25),
          quantile(stage.time.aids.art, probs = 0.25)),
    rbind(quantile(time.hivneg, probs = 0.5), # quantile(time.on.prep, probs = 0.5), 
          quantile(stage.time.ar.ndx, probs = 0.5),
          quantile(stage.time.af.ndx, probs = 0.5), quantile(stage.time.chronic.ndx, probs = 0.5), 
          quantile(stage.time.aids.ndx, probs = 0.5),
          quantile(stage.time.ar.dx, probs = 0.5),
          quantile(stage.time.af.dx, probs = 0.5), quantile(stage.time.chronic.dx, probs = 0.5), 
          quantile(stage.time.aids.dx, probs = 0.5),
          quantile(stage.time.ar.art, probs = 0.5),
          quantile(stage.time.af.art, probs = 0.5), quantile(stage.time.chronic.art, probs = 0.5), 
          quantile(stage.time.aids.art, probs = 0.5)),
    rbind(quantile(time.hivneg, probs = 0.75), # quantile(time.on.prep, probs = 0.75), 
          quantile(stage.time.ar.ndx, probs = 0.75),
          quantile(stage.time.af.ndx, probs = 0.75), quantile(stage.time.chronic.ndx, probs = 0.75), 
          quantile(stage.time.aids.ndx, probs = 0.75),
          quantile(stage.time.ar.dx, probs = 0.75),
          quantile(stage.time.af.dx, probs = 0.75), quantile(stage.time.chronic.dx, probs = 0.75), 
          quantile(stage.time.aids.dx, probs = 0.75),
          quantile(stage.time.ar.art, probs = 0.75),
          quantile(stage.time.af.art, probs = 0.75), quantile(stage.time.chronic.art, probs = 0.75), 
          quantile(stage.time.aids.art, probs = 0.75))))

rownames(df) <- c("time.hivneg", #"time.on.prep", #"time.off.prep", 
                  "stage.time.ar.ndx", "stage.time.af.ndx", "stage.time.chronic.ndx", "stage.time.aids.ndx",
                  "stage.time.ar.dx", "stage.time.af.dx", "stage.time.chronic.dx", "stage.time.aids.dx",
                  "stage.time.ar.art", "stage.time.af.art", "stage.time.chronic.art", "stage.time.aids.art")
colnames(df) <- c("Twentyfive", "Fifty", "Seventyfive")

##########################
df$weight <- c(1, 
               0.865, 0.865, 0.865, 0.72, 
               0.795, 0.795, 0.795, 0.72,
               0.83, 0.83, 0.83, 0.82)
df$Twentyfivevalue <- df$weight*df$Twentyfive
df$Fiftyvalue <- df$weight*df$Fifty
df$Seventyfivevalue <- df$weight*df$Seventyfive


## Create QALY/Cost Data Frame for First
QALY <- cbind(colSums(df[1:3] / 52),
              quantile(totalhivtests, probs = 0.25), quantile(totalhivtests, probs = 0.5), quantile(totalhivtests, probs = 0.75),
              quantile(totalhivtests.pos, probs = 0.25), quantile(totalhivtests.pos, probs = 0.5), quantile(totalhivtests.pos, probs = 0.75),
              quantile(stage.time.ar.art, probs = 0.25), quantile(stage.time.ar.art, probs = 0.5), quantile(stage.time.ar.art, probs = 0.75),
              quantile(stage.time.af.art, probs = 0.25), quantile(stage.time.af.art, probs = 0.5), quantile(stage.time.af.art, probs = 0.75),
              quantile(stage.time.chronic.art, probs = 0.25), quantile(stage.time.chronic.art, probs = 0.5), quantile(stage.time.chronic.art, probs = 0.75),
              quantile(stage.time.aids.art, probs = 0.25), quantile(stage.time.aids.art, probs = 0.5), quantile(stage.time.aids.art, probs = 0.75),
              quantile(totalrGCasympttests, probs = 0.25), quantile(totalrGCasympttests, probs = 0.5), quantile(totalrGCasympttests, probs = 0.75),
              quantile(totaluGCasympttests, probs = 0.25), quantile(totaluGCasympttests, probs = 0.5), quantile(totaluGCasympttests, probs = 0.75),
              quantile(totalGCasympttests, probs = 0.25), quantile(totalGCasympttests, probs = 0.5), quantile(totalGCasympttests, probs = 0.75),
              quantile(totalrCTasympttests, probs = 0.25), quantile(totalrGCasympttests, probs = 0.5), quantile(totalrGCasympttests, probs = 0.75),
              quantile(totaluCTasympttests, probs = 0.25), quantile(totaluCTasympttests, probs = 0.5), quantile(totaluCTasympttests, probs = 0.75),
              quantile(totalCTasympttests, probs = 0.25), quantile(totalCTasympttests, probs = 0.5), quantile(totalCTasympttests, probs = 0.75),
              quantile(totalsyphasympttests, probs = 0.25), quantile(totalsyphasympttests, probs = 0.5), quantile(totalsyphasympttests, probs = 0.75),
              quantile(totalrGCasympttests.pos, probs = 0.25), quantile(totalrGCasympttests.pos, probs = 0.5), quantile(totalrGCasympttests.pos, probs = 0.75),
              quantile(totaluGCasympttests.pos, probs = 0.25), quantile(totaluGCasympttests.pos, probs = 0.5), quantile(totaluGCasympttests.pos, probs = 0.75),
              quantile(totalGCasympttests.pos, probs = 0.25), quantile(totalGCasympttests.pos, probs = 0.5), quantile(totalGCasympttests.pos, probs = 0.75),
              quantile(totalrCTasympttests.pos, probs = 0.25), quantile(totalrGCasympttests.pos, probs = 0.5), quantile(totalrGCasympttests.pos, probs = 0.75),
              quantile(totaluCTasympttests.pos, probs = 0.25), quantile(totaluCTasympttests.pos, probs = 0.5), quantile(totaluCTasympttests.pos, probs = 0.75),
              quantile(totalCTasympttests.pos, probs = 0.25), quantile(totalCTasympttests.pos, probs = 0.5), quantile(totalCTasympttests.pos, probs = 0.75),
              quantile(totalsyphasympttests.pos, probs = 0.25), quantile(totalsyphasympttests.pos, probs = 0.5), quantile(totalsyphasympttests.pos, probs = 0.75))

QALY2 <- QALY
QALY2



# For others
sim <- truncate_sim(sim, at = 2600)

# Time in HIV stages
time.hivneg <- as.numeric(sim$epi$time.hivneg[520, ])
round(quantile(time.hivneg, probs = c(0.025, 0.5, 0.975)), 3)

time.on.prep <- as.numeric(sim$epi$time.on.prep[520, ])
round(quantile(time.on.prep, probs = c(0.025, 0.5, 0.975)), 3)

time.off.prep <- as.numeric(sim$epi$time.off.prep[520, ])
round(quantile(time.off.prep, probs = c(0.025, 0.5, 0.975)), 3)

stage.time.ar.ndx <- as.numeric(sim$epi$stage.time.ar.ndx[520, ]) / 52
round(quantile(stage.time.ar.ndx, probs = c(0.025, 0.5, 0.975)), 3)

stage.time.af.ndx <- as.numeric(sim$epi$stage.time.af.ndx[520, ]) / 52
round(quantile(stage.time.af.ndx, probs = c(0.025, 0.5, 0.975)), 3)

stage.time.chronic.ndx <- as.numeric(sim$epi$stage.time.chronic.ndx[520, ]) / 52
round(quantile(stage.time.chronic.ndx, probs = c(0.025, 0.5, 0.975)), 3)

stage.time.aids.ndx <- as.numeric(sim$epi$stage.time.aids.ndx[520,]) / 52
round(quantile(stage.time.aids.ndx, probs = c(0.025, 0.5, 0.975)), 3)

stage.time.ar.dx <- as.numeric(sim$epi$stage.time.ar.dx[520, ]) / 52
round(quantile(stage.time.ar.dx, probs = c(0.025, 0.5, 0.975)), 3)

stage.time.af.dx <- as.numeric(sim$epi$stage.time.af.dx[520, ]) / 52
round(quantile(stage.time.af.dx, probs = c(0.025, 0.5, 0.975)), 3)

stage.time.chronic.dx <- as.numeric(sim$epi$stage.time.chronic.dx[520, ]) / 52
round(quantile(stage.time.chronic.dx, probs = c(0.025, 0.5, 0.975)), 3)

stage.time.aids.dx <- as.numeric(sim$epi$stage.time.aids.dx[520,]) / 52
round(quantile(stage.time.aids.dx, probs = c(0.025, 0.5, 0.975)), 3)

stage.time.ar.art <- as.numeric(sim$epi$stage.time.ar.art[520, ]) / 52
round(quantile(stage.time.ar.art, probs = c(0.025, 0.5, 0.975)), 3)

stage.time.af.art <- as.numeric(sim$epi$stage.time.af.art[520, ]) / 52
round(quantile(stage.time.af.art, probs = c(0.025, 0.5, 0.975)), 3)

stage.time.chronic.art <- as.numeric(sim$epi$stage.time.chronic.art[520, ]) / 52
round(quantile(stage.time.chronic.art, probs = c(0.025, 0.5, 0.975)), 3)

stage.time.aids.art <- as.numeric(sim$epi$stage.time.aids.art[520,]) / 52
round(quantile(stage.time.aids.art, probs = c(0.025, 0.5, 0.975)), 3)


# Total HIV Tests (PrEP or not)
totalhivtests <- as.numeric(sim$epi$totalhivtests[520, ])
round(quantile(totalhivtests, probs = c(0.025, 0.5, 0.975)), 3)


# Total Symptomatic STI Tests
totalrGCsympttests <- as.numeric(sim$epi$totalrGCsympttests[520, ])
round(quantile(totalrGCsympttests, probs = c(0.025, 0.5, 0.975)), 3)

totaluGCsympttests <- as.numeric(sim$epi$totaluGCsympttests[520, ])
round(quantile(totaluGCsympttests, probs = c(0.025, 0.5, 0.975)), 3)

totalGCsympttests <- as.numeric(sim$epi$totalGCsympttests[520, ])
round(quantile(totalGCsympttests, probs = c(0.025, 0.5, 0.975)), 3)

totalrCTsympttests <- as.numeric(sim$epi$totalrCTsympttests[520, ])
round(quantile(totalrCTsympttests, probs = c(0.025, 0.5, 0.975)), 3)

totaluCTsympttests <- as.numeric(sim$epi$totaluCTsympttests[520, ])
round(quantile(totaluCTsympttests, probs = c(0.025, 0.5, 0.975)), 3)

totalCTsympttests <- as.numeric(sim$epi$totalCTsympttests[520, ])
round(quantile(totalCTsympttests, probs = c(0.025, 0.5, 0.975)), 3)

totalsyphsympttests <- as.numeric(sim$epi$totalsyphsympttests[520, ])
round(quantile(totalsyphsympttests, probs = c(0.025, 0.5, 0.975)), 3)


# Total Asymptomatic STI Tests
totalrGCasympttests <- as.numeric(sim$epi$totalrGCasympttests[520, ])
round(quantile(totalrGCasympttests, probs = c(0.025, 0.5, 0.975)), 3)

totaluGCasympttests <- as.numeric(sim$epi$totaluGCasympttests[520, ])
round(quantile(totaluGCasympttests, probs = c(0.025, 0.5, 0.975)), 3)

totalGCasympttests <- as.numeric(sim$epi$totalGCasympttests[520, ])
round(quantile(totalGCasympttests, probs = c(0.025, 0.5, 0.975)), 3)

totalrCTasympttests <- as.numeric(sim$epi$totalrCTasympttests[520, ])
round(quantile(totalrCTasympttests, probs = c(0.025, 0.5, 0.975)), 3)

totaluCTasympttests <- as.numeric(sim$epi$totaluCTasympttests[520, ])
round(quantile(totaluCTasympttests, probs = c(0.025, 0.5, 0.975)), 3)

totalCTasympttests <- as.numeric(sim$epi$totalCTasympttests[520, ])
round(quantile(totalCTasympttests, probs = c(0.025, 0.5, 0.975)), 3)

totalsyphasympttests <- as.numeric(sim$epi$totalsyphasympttests[520, ])
round(quantile(totalsyphasympttests, probs = c(0.025, 0.5, 0.975)), 3)


# Total Positive HIV Tests
totalhivtests.pos <- as.numeric(sim$epi$totalhivtests.pos[520, ])
round(quantile(totalhivtests, probs = c(0.025, 0.5, 0.975)), 3)


# Total Positive Asymptomatic STI Tests
totalrGCasympttests.pos <- as.numeric(sim$epi$totalrGCasympttests.pos[520, ])
round(quantile(totalrGCasympttests, probs = c(0.025, 0.5, 0.975)), 3)

totaluGCasympttests.pos <- as.numeric(sim$epi$totaluGCasympttests.pos[520, ])
round(quantile(totaluGCasympttests, probs = c(0.025, 0.5, 0.975)), 3)

totalGCasympttests.pos <- as.numeric(sim$epi$totalGCasympttests.pos[520, ])
round(quantile(totalGCasympttests, probs = c(0.025, 0.5, 0.975)), 3)

totalrCTasympttests.pos <- as.numeric(sim$epi$totalrCTasympttests.pos[520, ])
round(quantile(totalrCTasympttests, probs = c(0.025, 0.5, 0.975)), 3)

totaluCTasympttests.pos <- as.numeric(sim$epi$totaluCTasympttests.pos[520, ])
round(quantile(totaluCTasympttests, probs = c(0.025, 0.5, 0.975)), 3)

totalCTasympttests.pos <- as.numeric(sim$epi$totalCTasympttests.pos[520, ])
round(quantile(totalCTasympttests, probs = c(0.025, 0.5, 0.975)), 3)

totalsyphasympttests.pos <- as.numeric(sim$epi$totalsyphasympttests.pos[520, ])
round(quantile(totalsyphasympttests, probs = c(0.025, 0.5, 0.975)), 3)


## Create QALY Data Frame
df <- as.data.frame(cbind(rbind(quantile(time.hivneg, probs = 0.25), # quantile(time.on.prep, probs = 0.25), 
                                quantile(stage.time.ar.ndx, probs = 0.25),
                                quantile(stage.time.af.ndx, probs = 0.25), quantile(stage.time.chronic.ndx, probs = 0.25),
                                quantile(stage.time.aids.ndx, probs = 0.25),
                                quantile(stage.time.ar.dx, probs = 0.25),
                                quantile(stage.time.af.dx, probs = 0.25), quantile(stage.time.chronic.dx, probs = 0.25),
                                quantile(stage.time.aids.dx, probs = 0.25)),
                          rbind(quantile(time.hivneg, probs = 0.5), # quantile(time.on.prep, probs = 0.5), 
                                quantile(stage.time.ar.ndx, probs = 0.5),
                                quantile(stage.time.af.ndx, probs = 0.5), quantile(stage.time.chronic.ndx, probs = 0.5), 
                                quantile(stage.time.aids.ndx, probs = 0.5),
                                quantile(stage.time.ar.dx, probs = 0.5),
                                quantile(stage.time.af.dx, probs = 0.5), quantile(stage.time.chronic.dx, probs = 0.5), 
                                quantile(stage.time.aids.dx, probs = 0.5)),
                          rbind(quantile(time.hivneg, probs = 0.75), # quantile(time.on.prep, probs = 0.75), 
                                quantile(stage.time.ar.ndx, probs = 0.75),
                                quantile(stage.time.af.ndx, probs = 0.75), quantile(stage.time.chronic.ndx, probs = 0.75), 
                                quantile(stage.time.aids.ndx, probs = 0.75),
                                quantile(stage.time.ar.dx, probs = 0.75),
                                quantile(stage.time.af.dx, probs = 0.75), quantile(stage.time.chronic.dx, probs = 0.75), 
                                quantile(stage.time.aids.dx, probs = 0.75))))
rownames(df) <- c("time.hivneg", # "time.on.prep", #"time.off.prep", 
                  "stage.time.ar.ndx", "stage.time.af.ndx", "stage.time.chronic.ndx", "stage.time.aids.ndx",
                  "stage.time.ar.dx", "stage.time.af.dx", "stage.time.chronic.dx", "stage.time.aids.dx",
                  "stage.time.ar.art", "stage.time.af.art", "stage.time.chronic.art", "stage.time.aids.art")
colnames(df) <- c("Twentyfive", "Fifty", "Seventyfive")

######
df$weight <- c(1, 
               0.865, 0.865, 0.865, 0.72, 
               0.795, 0.795, 0.795, 0.72,
               0.83, 0.83, 0.83, 0.82)
df$Twentyfivevalue <- df$weight*df$Twentyfive
df$Fiftyvalue <- df$weight*df$Fifty
df$Seventyfivevalue <- df$weight*df$Seventyfive


## Create QALY/Cost Data Frame
QALY <- cbind(colSums(df[1:3] / 52),
              quantile(totalhivtests, probs = 0.25), quantile(totalhivtests, probs = 0.5), quantile(totalhivtests, probs = 0.75),
              quantile(totalhivtests.pos, probs = 0.25), quantile(totalhivtests.pos, probs = 0.5), quantile(totalhivtests.pos, probs = 0.75),
              quantile(stage.time.ar.art, probs = 0.25), quantile(stage.time.ar.art, probs = 0.5), quantile(stage.time.ar.art, probs = 0.75),
              quantile(stage.time.af.art, probs = 0.25), quantile(stage.time.af.art, probs = 0.5), quantile(stage.time.af.art, probs = 0.75),
              quantile(stage.time.chronic.art, probs = 0.25), quantile(stage.time.chronic.art, probs = 0.5), quantile(stage.time.chronic.art, probs = 0.75),
              quantile(stage.time.aids.art, probs = 0.25), quantile(stage.time.aids.art, probs = 0.5), quantile(stage.time.aids.art, probs = 0.75),
              quantile(totalrGCasympttests, probs = 0.25), quantile(totalrGCasympttests, probs = 0.5), quantile(totalrGCasympttests, probs = 0.75),
              quantile(totaluGCasympttests, probs = 0.25), quantile(totaluGCasympttests, probs = 0.5), quantile(totaluGCasympttests, probs = 0.75),
              quantile(totalGCasympttests, probs = 0.25), quantile(totalGCasympttests, probs = 0.5), quantile(totalGCasympttests, probs = 0.75),
              quantile(totalrCTasympttests, probs = 0.25), quantile(totalrGCasympttests, probs = 0.5), quantile(totalrGCasympttests, probs = 0.75),
              quantile(totaluCTasympttests, probs = 0.25), quantile(totaluCTasympttests, probs = 0.5), quantile(totaluCTasympttests, probs = 0.75),
              quantile(totalCTasympttests, probs = 0.25), quantile(totalCTasympttests, probs = 0.5), quantile(totalCTasympttests, probs = 0.75),
              quantile(totalsyphasympttests, probs = 0.25), quantile(totalsyphasympttests, probs = 0.5), quantile(totalsyphasympttests, probs = 0.75),
              quantile(totalrGCasympttests.pos, probs = 0.25), quantile(totalrGCasympttests.pos, probs = 0.5), quantile(totalrGCasympttests.pos, probs = 0.75),
              quantile(totaluGCasympttests.pos, probs = 0.25), quantile(totaluGCasympttests.pos, probs = 0.5), quantile(totaluGCasympttests.pos, probs = 0.75),
              quantile(totalGCasympttests.pos, probs = 0.25), quantile(totalGCasympttests.pos, probs = 0.5), quantile(totalGCasympttests.pos, probs = 0.75),
              quantile(totalrCTasympttests.pos, probs = 0.25), quantile(totalrGCasympttests.pos, probs = 0.5), quantile(totalrGCasympttests.pos, probs = 0.75),
              quantile(totaluCTasympttests.pos, probs = 0.25), quantile(totaluCTasympttests.pos, probs = 0.5), quantile(totaluCTasympttests.pos, probs = 0.75),
              quantile(totalCTasympttests.pos, probs = 0.25), quantile(totalCTasympttests.pos, probs = 0.5), quantile(totalCTasympttests.pos, probs = 0.75),
              quantile(totalsyphasympttests.pos, probs = 0.25), quantile(totalsyphasympttests.pos, probs = 0.5), quantile(totalsyphasympttests.pos, probs = 0.75))
QALY2 <- rbind(QALY2, QALY)


# Create Table for Output
#rownames(QALY2) <- c("Cov = 0%", "Cov = 10%", "Cov = 20%", "Cov = 30%", "Cov = 40%", "Cov = 50%",
#                     "Cov = 60%", "Cov = 70%", "Cov = 80%", "Cov = 90%", "Cov = 100%")
colnames(QALY2) <- c("25% QALY", "50% QALY", "75% QALY", 
                     "25% HIV Tests", "50% HIV Tests", "75% HIV Tests", 
                     "25% HIV Positive Tests", "50% HIV Positive Tests", "75% HIV Positive Tests", 
                     "25% Years on PrEP", "50% Years on PrEP", "75% Years on PrEP",
                     "25% AR Time on ART", "50% AR Time on ART", "75% AR Time on ART",
                     "25% AF Time on ART", "50% AF Time on ART", "75% AF Time on ART",
                     "25% Chronic Time on ART", "50% Chronic Time on ART", "75% Chronic Time on ART",
                     "25% AIDS Time on ART", "50% AIDS Time on ART", "75% AIDS Time on ART",
                     "25% RGC Tests", "50% RGC Tests", "75% RGC Tests",
                     "25% UGC Tests", "50% UGC Tests", "75% UGC Tests",
                     "25% GC Tests", "50% GC Tests", "75% GC Tests",
                     "25% RCT Tests", "50% RCT Tests", "75% RCT Tests",
                     "25% UCT Tests", "50% UCT Tests", "75% UCT Tests",
                     "25% CT Tests", "50% CT Tests", "75% CT Tests",
                     "25% Syph Tests", "50% Syph Tests", "75% Syph Tests",
                     "25% RGC Positive  Tests", "50% RGC Positive  Tests", "75% RGC Positive Tests",
                     "25% UGC Positive Tests", "50% UGC Positive Tests", "75% UGC Positive Tests",
                     "25% GC Positive Tests", "50% GC Positive Tests", "75% GC Positive Tests",
                     "25% RCT Positive Tests", "50% RCT Positive Tests", "75% RCT Positive Tests",
                     "25% UCT Positive Tests", "50% UCT Positive Tests", "75% UCT Positive Tests",
                     "25% CT Positive Tests", "50% CT Positive Tests", "75% CT Positive Tests",
                     "25% Syph Positive Tests", "50% Syph Positive Tests", "75% Syph Positive Tests")

QALY2

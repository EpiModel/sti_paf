
## Process burn-in
library("EpiModelHPC")
library("EpiModelHIV")

rm(list = ls())
# Merge sim files
sim <- merge_simfiles(simno = 331, indir = "data/", ftype = "max")

# Create function for selecting sim closest to target
mean_sim <- function(sim, targets) {

  nsims <- sim$control$nsims

  # Initialize distance vector
  dist <- rep(NA, nsims)
  #dist2 <- rep(NA, nsims)
  dist3 <- rep(NA, nsims)
  prev.hiv <- rep(NA, nsims)
  prev.pssyph <- rep(NA, nsims)
  prev.gc <- rep(NA, nsims)
  prev.ct <- rep(NA, nsims)

  # Obtain statistics and perform multivariable Euclidean distance calculation
  for (i in 1:nsims) {

      # Create data frame to draw statistics from
      df <- as.data.frame(x = sim, out = "vals", sim = i)

      # Create a vector of statistics
      calib <- c(mean(tail(df$ir100.gc, 10)),
                 mean(tail(df$ir100.ct, 10)),
                 mean(tail(df$ir100.syph, 10)),
                 mean(tail(df$i.prev, 10)),
                 mean(tail(df$prev.syph, 10)),
                 mean(tail(df$prev.primsecosyph, 10)),
                 mean(df$ir100.syph[5200] - df$ir100.syph[5190]),
                 mean(df$ir100.gc[5200] - df$ir100.gc[5190]),
                 mean(df$ir100.ct[5200] - df$ir100.ct[5190]),
                 mean(df$ir100.syph[5190] - df$ir100.syph[5180]),
                 mean(df$ir100.gc[5190] - df$ir100.gc[5180]),
                 mean(df$ir100.ct[5190] - df$ir100.ct[5180]))

      #Syph
      wts <- c(2, 2, 2, 2, 2, 2,
               1, 1, 1, 1, 1, 1)

      # No Syph
      # wts <- c(2, 2, 2,
      #          1, 1, 1, 1)

      # Iteratively calculate distance
      dist[i] <- sqrt(sum(((calib - targets)*wts)^2))
      #dist2[i] <- mean(abs(wts * (calib - targets)/targets))
      dist3[i] <- mean(abs((calib - targets)/targets))
      prev.hiv[i] <- df$i.prev[5200]
      prev.pssyph[i] <- df$prev.primsecosyph[5200]
      prev.gc[i] <- df$prev.gc[5200]
      prev.ct[i] <- df$prev.ct[5200]
  }

  # Which sims minimizes distance
  meansims <<- which(dist < 4 & prev.hiv != 0 & prev.pssyph != 0 &
                       prev.gc != 0 & prev.ct != 0)
  meansims2 <<- which(dist3 < 0.2 & prev.hiv != 0 & prev.pssyph != 0 &
                        prev.gc != 0 & prev.ct != 0)

  return(meansims2)
}

# Run function
mean_sim(sim, targets = c(3.5, 5.6, 2.6,
                          0.15, 0.012, 0.006,
                          0, 0, 0, 0))

# Save burn-in file for FU sims
sim2 <- get_sims(sim, sims = meansims2)


# Syphilis
par(mfrow = c(2,2), oma = c(0,0,2,0))
plot(sim, y = "ir100.syph")
abline(h = 2.6, col = "red", lty = 2)
title("Syph IR")
plot(sim, y = "prev.primsecosyph")
abline(h = 0.006, lty = c(2), col = 'red')
title("P&S Syphilis Prevalence")
plot(sim, y = "prev.syph")
abline(h = 0.012, col = "red", lty = 2)
title("Syphilis (All Stages) Prevalence")
title(paste0("All Sims ", sim$control$simno), outer = TRUE)

#NG/CT
par(mfrow = c(2,2), oma = c(0,0,2,0))
plot(sim, y = "ir100.ct")
abline(h = 5.6, col = "red", lty = 2)
title("CT IR")
plot(sim, y = "prev.ct")
title("CT Prev")
plot(sim, y = "ir100.gc")
abline(h = 3.5, col = "red", lty = 2)
title("NG IR")
plot(sim, y = "prev.gc")
title("NG Prev")
title(paste0("All Sims ", sim$control$simno), outer = TRUE)

#HIV
par(mfrow = c(1,2), oma = c(0,0,2,0))
plot(sim, y = "ir100")
title("HIV IR")
abline(h = 2.0, lty = 2, col = "red")
plot(sim, y = "i.prev")
title("HIV Prev")
abline(h = 0.15, lty = 2, col = "red")
title(paste0("All Sims ", sim$control$simno), outer = TRUE)

# Overall
par(mfrow = c(2, 2), oma = c(0,0,2,0))
plot(sim, y = "i.prev", qnts = 0.90)
abline(h = 0.15, col = "red", lty = 2)
title("HIV Prevalence")
plot(sim, y = "ir100.gc", ylim = c(0, 15))
abline(h = 3.5, col = "red", lty = 2)
title("NG Incidence")
plot(sim, y = "ir100.ct", ylim = c(0, 15))
abline(h = 5.6, col = "red", lty = 2)
title("CT Incidence")
plot(sim, y = "ir100.syph", ylim = c(0, 15))
abline(h = 2.6, col = "red", lty = 2)
#abline(h = 1.5, col = "red", lty = 2)
title("Syph Incidence")
title(paste0("All Sims ", sim$control$simno), outer = TRUE)


# Best-fitting ------------------
par(mfrow = c(2,2), oma = c(0,0,2,0))
# plot(sim2, y = "ir100")
# abline(h = 3.8, col = "red", lty = 2)
plot(sim2, y = "i.prev")
abline(h = 0.15, col = "red", lty = 2)
title("HIV Prevalence")
plot(sim2, y = "ir100.gc")
abline(h = 3.5, col = "red", lty = 2)
title("GC Incidence")
plot(sim2, y = "ir100.ct")
abline(h = 5.6, col = "red", lty = 2)
title("CT Incidence")
plot(sim2, y = "ir100.syph")
abline(h = 2.6, col = "red", lty = 2)
title("Syph Incidence")
title(paste0("Best-fitting Sim ", sim2$control$simno), outer = TRUE)

# Other statistics--------------------------------------------------------------

# Syphilis prevalence and ratios
par(mfrow = c(1,2), oma = c(0,0,2,0))
plot(sim2, y = "prev.primsecosyph", qnts = 0.90)
abline(h = 0.01, lty = c(2), col = 'red')
abline(h = 0.006, lty = c(2), col = 'red')
title("P&S Syphilis Prevalence")
plot(sim2, y = "prev.syph", qnts = 0.90)
abline(h = 0.02, col = "red", lty = 2)
abline(h = 0.012, col = "red", lty = 2)
title("Syphilis (All Stages) Prevalence")

par(mfrow = c(1,2), oma = c(0,0,2,0))
plot(sim, y = "prev.gc", qnts = 0.90)
title("NG Prevalence")
plot(sim, y = "prev.ct", qnts = 0.90)
title("CT Prevalence")

plot(sim, y = "early.late.syphratio", ylim = c(0, 1.0))
title("Ratio of Early to Late \n Syphilis Cases")
#abline(h = 0.2, lty = c(2), col = 'red')
plot(sim, y = "early.late.diagsyphratio", ylim = c(0, 10))
title("Ratio of Diagnosed Early to Late \n Syphilis Cases")
#abline(h = 0.5, lty = c(2), col = 'red')
title("Syphilis Prevalence Measures", outer = TRUE)



par(mfrow = c(1,2), oma = c(0,0,2,0))
plot(sim2, y = "prev.primsecosyph", qnts = 0.90)
abline(h = 0.01, lty = c(2), col = 'red')
title("P&S Syphilis Prevalence")
plot(sim2, y = "prev.syph", qnts = 0.90)
title("Syphilis Prevalence")
abline(h = 0.02, col = "red", lty = 2)
abline(h = 0.012, col = "red", lty = 2)

par(mfrow = c(1,2), oma = c(0,0,2,0))
plot(sim2, y = "prev.gc", qnts = 0.90)
title("NG Prevalence")
plot(sim2, y = "prev.ct", qnts = 0.90)
title("CT Prevalence")

## Tested in Last 12 months by serostatus
par(mfrow = c(2, 2))
plot(sim, y = 'test.gc.12mo', ylab = "Proportion")
plot(sim, y = 'test.gc.12mo.hivpos', add = TRUE, mean.col = "red", qnts.col = "red")
plot(sim, y = 'test.gc.12mo.hivneg', add = TRUE, mean.col = "green", qnts.col = "green")
abline(h = c(0.641, 0.462), lty = c(2, 2), col = c("red", "green"))
legend("topright", legend = c("All", "HIV+", "HIV-"),
       lty = c(1, 1, 1), col = c("blue", "red", "green"))
title("Tested for NG in last 12 months")

plot(sim, y = 'test.ct.12mo', ylab = "Proportion")
plot(sim, y = 'test.ct.12mo.hivpos', add = TRUE, mean.col = "red", qnts.col = "red")
plot(sim, y = 'test.ct.12mo.hivneg', add = TRUE, mean.col = "green", qnts.col = "green")
abline(h = c(0.628, 0.458), lty = c(2, 2), col = c("red", "green"))
legend("topright", legend = c("All", "HIV+", "HIV-"),
       lty = c(1, 1, 1), col = c("blue", "red", "green"))
title("Tested for CT in last 12 months")

plot(sim, y = 'test.syph.12mo', ylab = "Proportion")
plot(sim, y = 'test.syph.12mo.hivpos', add = TRUE, mean.col = "red", qnts.col = "red")
plot(sim, y = 'test.syph.12mo.hivneg', add = TRUE, mean.col = "green", qnts.col = "green")
abline(h = c(0.68, 0.45), lty = c(2, 2), col = c("red", "green"))
legend("topright", legend = c("All", "HIV+", "HIV"),
       lty = c(1, 1, 1), col = c("blue", "red", "green"))
title("Tested for Syph in last 12 months")

# NG: 46.2% (Hoots unpublished NHBS)
# CT: 45.8 % (Hoots unpublished NHBS)
# Syphilis: 45% (An 2017 - NHBS)

#   HIV-positive MSM:
# NG: 43% (Mattson 2017 - MMP), 47.2% (Patel 2017 – MMP)
# CT: 43% (Mattson 2017 - MMP), 47.2% (Patel 2017 – MMP)
# Syphilis: 69% (Mattson 2017 - MMP)

# Ratio approach
# NG among HIV-negative MSM: (47.2) x (46.2/64.1) = 34.0%
# CT among HIV-negative MSM:(47.2) x (45.8 /62.8) = 34.4%
## Tested in Last 12 months by diagnosis status
# par(mfrow = c(2, 2))
# plot(sim, y = 'test.gc.12mo', ylab = "Proportion")
# plot(sim, y = 'test.gc.12mo.hivdiag', add = TRUE, mean.col = "red", qnts.col = "red")
# plot(sim, y = 'test.gc.12mo.nonhivdiag', add = TRUE, mean.col = "green", qnts.col = "green")
# abline(h = c(0.472, 0.340), lty = c(2, 2), col = c("red", "green"))
# legend("topright", legend = c("All", "HIV-diag", "Non HIV-diag"),
#        lty = c(1, 1, 1), col = c("blue", "red", "green"), ncol = 2)
# title("Tested for NG in last 12 months")
#
# plot(sim, y = 'test.ct.12mo', ylab = "Proportion")
# plot(sim, y = 'test.ct.12mo.hivdiag', add = TRUE, mean.col = "red", qnts.col = "red")
# plot(sim, y = 'test.ct.12mo.nonhivdiag', add = TRUE, mean.col = "green", qnts.col = "green")
# abline(h = c(0.472, 0.344), lty = c(2, 2), col = c("red", "green"))
# legend("topright", legend = c("All", "HIV-diag", "Non HIV-diag"),
#        lty = c(1, 1, 1), col = c("blue", "red", "green"), ncol = 2)
# title("Tested for CT in last 12 months")
#
# plot(sim, y = 'test.syph.12mo', ylab = "Proportion")
# plot(sim, y = 'test.syph.12mo.hivdiag', add = TRUE, mean.col = "red", qnts.col = "red")
# plot(sim, y = 'test.syph.12mo.nonhivdiag', add = TRUE, mean.col = "green", qnts.col = "green")
# abline(h = c(0.69, 0.45), lty = c(2, 2), col = c("red", "green"))
# legend("topright", legend = c("All", "HIV-diag", "Non HIV-diag"),
#        lty = c(1, 1, 1), col = c("blue", "red", "green"), ncol = 2)
# title("Tested for Syph in last 12 months")

## STI asymptomatic Testing - all sims
par(mfrow = c(2,2), oma = c(0,0,2,0))
plot(sim, y = "GCasympttests", mean.col = "black", qnts.col = "black", qnts = 0.5, ylab = "Tests/week")
plot(sim, y = "GCasympttests.hivneg", mean.col = "red", qnts.col = "red", qnts = 0.5, add = TRUE)
plot(sim, y = "GCasympttests.hivpos", mean.col = "blue", qnts.col = "blue", qnts = 0.5, add = TRUE)
plot(sim, y = "GCasympttests.pos.hivneg", mean.col = "purple", qnts.col = "purple", qnts = 0.5, add = TRUE)
plot(sim, y = "GCasympttests.pos.hivpos", mean.col = "yellow", qnts.col = "yellow", qnts = 0.5, add = TRUE)
legend("topleft", lty = c(1, 1, 1, 1, 1), col = c("black", "red", "blue", "purple", "yellow"),
       legend = c("All", "NG-non-HIV dx", "NG-HIV dx", "NG-pos non-HIV", "NG-pos HIV dx"),
       ncol = 2)

plot(sim, y = "CTasympttests", mean.col = "black", qnts.col = "black", qnts = 0.5, ylab = "Tests/week")
plot(sim, y = "CTasympttests.hivneg", mean.col = "green", qnts.col = "green", qnts = 0.5, add = TRUE)
plot(sim, y = "CTasympttests.hivpos", mean.col = "orange", qnts.col = "orange", qnts = 0.5, add = TRUE)
plot(sim, y = "CTasympttests.pos.hivneg", mean.col = "purple", qnts.col = "purple", qnts = 0.5, add = TRUE)
plot(sim, y = "CTasympttests.pos.hivpos", mean.col = "yellow", qnts.col = "yellow", qnts = 0.5, add = TRUE)
legend("topleft", lty = c(1, 1, 1, 1, 1), col = c("black", "green", "orange", "purple", "yellow"),
       legend = c("All", "CT-non-HIV dx", "CT-HIV dx", "CT-pos non-HIV", "CT-pos HIV dx"),
       ncol = 2)

plot(sim, y = "syphasympttests", mean.col = "black", qnts.col = "black", qnts = 0.5, ylab = "Tests/week")
plot(sim, y = "syphasympttests.hivneg", mean.col = "purple", qnts.col = "purple", qnts = 0.5, add = TRUE)
plot(sim, y = "syphasympttests.hivpos", mean.col = "brown", qnts.col = "brown", qnts = 0.5, add = TRUE)
plot(sim, y = "syphasympttests.pos.hivneg", mean.col = "green", qnts.col = "green", qnts = 0.5, add = TRUE)
plot(sim, y = "syphasympttests.pos.hivpos", mean.col = "orange", qnts.col = "orange", qnts = 0.5, add = TRUE)
legend("topleft", lty = c(1, 1, 1, 1, 1), col = c("black", "purple", "brown", "green", "orange"),
       legend = c("All", "Syph-non-HIV dx", "Syph-HIV dx", "Syph-pos non-HIV", "Syph-pos HIV dx"),
       ncol = 2)

plot(sim, y = "stiasympttests", mean.col = "black", qnts.col = "black", qnts = 0.5, ylab = "Tests/week", ylim = c(0, 600))
plot(sim, y = "stiasympttests.hivneg", mean.col = "purple", qnts.col = "purple", qnts = 0.5, add = TRUE)
plot(sim, y = "stiasympttests.hivpos", mean.col = "brown", qnts.col = "brown", qnts = 0.5, add = TRUE)
plot(sim, y = "stiasympttests.pos.hivneg", mean.col = "green", qnts.col = "green", qnts = 0.5, add = TRUE)
plot(sim, y = "stiasympttests.pos.hivpos", mean.col = "orange", qnts.col = "orange", qnts = 0.5, add = TRUE)
legend("topleft", lty = c(1, 1, 1, 1, 1), col = c("black", "purple", "brown", "green", "orange"),
       legend = c("All", "STI non-HIV dx", "STI HIV dx", "STI-pos non-HIV", "STI-pos HIV-dx"), ncol = 2)
title("STI Testing - Serostatus-Specific", outer = TRUE)

# STI Symptomatic testing
par(mfrow = c(1, 1))
plot(sim, y = "GCsympttests", mean.col = "black", qnts.col = "black", qnts = 0.5, ylab = "Tests/week")
plot(sim, y = "CTsympttests", mean.col = "purple", qnts.col = "purple", qnts = 0.5, add = TRUE)
plot(sim, y = "stisympttests", mean.col = "brown", qnts.col = "brown", qnts = 0.5, add = TRUE)
legend("topleft", lty = c(1, 1, 1), col = c("black", "purple", "brown"),
       legend = c("NG", "CT", "STI"))
title("STI Symptomatic Testing", outer = TRUE)

## STI Testing - mean sim
par(mfrow = c(2,2), oma = c(0,0,2,0))
plot(sim2, y = "GCasympttests", mean.col = "black", ylab = "Tests/week")
plot(sim2, y = "GCasympttests.hivneg", mean.col = "red", add = TRUE)
plot(sim2, y = "GCasympttests.hivpos", mean.col = "blue", add = TRUE)
legend("topleft", lty = c(1, 1, 1), col = c("black", "red", "blue"),
       legend = c("All", "NG-non-HIV dx", "NG-HIV dx"))
plot(sim2, y = "GCasympttests", mean.col = "black", ylab = "Tests/week")
plot(sim2, y = "CTasympttests.hivneg", mean.col = "green", add = TRUE)
plot(sim2, y = "CTasympttests.hivpos", mean.col = "orange", add = TRUE)
legend("topleft", lty = c(1, 1, 1), col = c("black","green", "orange"),
       legend = c("All", "CT-non-HIV dx", "CT-HIV dx"))
plot(sim2, y = "syphasympttests", mean.col = "black", ylab = "Tests/week")
plot(sim2, y = "syphasympttests.hivneg", mean.col = "purple", add = TRUE)
plot(sim2, y = "syphasympttests.hivpos", mean.col = "brown", add = TRUE)
legend("topleft", lty = c(1, 1, 1), col = c("black", "purple", "brown"),
       legend = c("All", "Syph-non-HIV dx", "Syph-HIV dx"))
#title("STI Testing - 10% Coverage \n Serostatus-specific")
plot(sim2, y = "stiasympttests", mean.col = "black", ylab = "Tests/week")
plot(sim2, y = "stiasympttests.hivneg", mean.col = "purple", add = TRUE)
plot(sim2, y = "stiasympttests.hivpos", mean.col = "brown", add = TRUE)
legend("topleft", lty = c(1, 1, 1), col = c("black", "purple", "brown"),
       legend = c("All", "Syph-non-HIV dx", "Syph-HIV dx"))
title("STI Testing - Serostatus-Specific", outer = TRUE)

# Asymptomatic testing
par(mfrow = c(2,2), oma = c(0,0,2,0))
plot(sim2, y = "GCasympttests", mean.col = "black", qnts.col = "black", qnts = 0.5, ylab = "Tests/week")
plot(sim2, y = "GCasympttests.hivneg", mean.col = "red", qnts.col = "red", qnts = 0.5, add = TRUE)
plot(sim2, y = "GCasympttests.hivpos", mean.col = "blue", qnts.col = "blue", qnts = 0.5, add = TRUE)
plot(sim2, y = "GCasympttests.pos.hivneg", mean.col = "purple", qnts.col = "purple", qnts = 0.5, add = TRUE)
plot(sim2, y = "GCasympttests.pos.hivpos", mean.col = "yellow", qnts.col = "yellow", qnts = 0.5, add = TRUE)
legend("topleft", lty = c(1, 1, 1, 1, 1), col = c("black", "red", "blue", "purple", "yellow"),
       legend = c("All", "NG-non-HIV dx", "NG-HIV dx", "NG-pos non-HIV", "NG-pos HIV dx"),
       ncol = 2)

plot(sim2, y = "CTasympttests", mean.col = "black", qnts.col = "black", qnts = 0.5, ylab = "Tests/week")
plot(sim2, y = "CTasympttests.hivneg", mean.col = "green", qnts.col = "green", qnts = 0.5, add = TRUE)
plot(sim2, y = "CTasympttests.hivpos", mean.col = "orange", qnts.col = "orange", qnts = 0.5, add = TRUE)
plot(sim2, y = "CTasympttests.pos.hivneg", mean.col = "purple", qnts.col = "purple", qnts = 0.5, add = TRUE)
plot(sim2, y = "CTasympttests.pos.hivpos", mean.col = "yellow", qnts.col = "yellow", qnts = 0.5, add = TRUE)
legend("topleft", lty = c(1, 1, 1, 1, 1), col = c("black", "green", "orange", "purple", "yellow"),
       legend = c("All", "CT-non-HIV dx", "CT-HIV dx", "CT-pos non-HIV", "CT-pos HIV dx"),
       ncol = 2)

plot(sim2, y = "syphasympttests", mean.col = "black", qnts.col = "black", qnts = 0.5, ylab = "Tests/week")
plot(sim2, y = "syphasympttests.hivneg", mean.col = "purple", qnts.col = "purple", qnts = 0.5, add = TRUE)
plot(sim2, y = "syphasympttests.hivpos", mean.col = "brown", qnts.col = "brown", qnts = 0.5, add = TRUE)
plot(sim2, y = "syphasympttests.pos.hivneg", mean.col = "green", qnts.col = "green", qnts = 0.5, add = TRUE)
plot(sim2, y = "syphasympttests.pos.hivpos", mean.col = "orange", qnts.col = "orange", qnts = 0.5, add = TRUE)
legend("topleft", lty = c(1, 1, 1, 1, 1), col = c("black", "purple", "brown", "green", "orange"),
       legend = c("All", "Syph-non-HIV dx", "Syph-HIV dx", "Syph-pos non-HIV", "Syph-pos HIV dx"),
       ncol = 2)

plot(sim2, y = "stiasympttests", mean.col = "black", qnts.col = "black", qnts = 0.5, ylab = "Tests/week")
plot(sim2, y = "stiasympttests.hivneg", mean.col = "purple", qnts.col = "purple", qnts = 0.5, add = TRUE)
plot(sim2, y = "stiasympttests.hivpos", mean.col = "brown", qnts.col = "brown", qnts = 0.5, add = TRUE)
plot(sim2, y = "stiasympttests.pos.hivneg", mean.col = "green", qnts.col = "green", qnts = 0.5, add = TRUE)
plot(sim2, y = "stiasympttests.pos.hivpos", mean.col = "orange", qnts.col = "orange", qnts = 0.5, add = TRUE)
legend("topleft", lty = c(1, 1, 1, 1, 1), col = c("black", "purple", "brown", "green", "orange"),
       legend = c("All", "STI non-HIV dx", "STI HIV dx", "STI-pos non-HIV", "STI-pos HIV-dx"), ncol = 2)
title("STI Testing - Serostatus-Specific", outer = TRUE)

# STI Symptomatic testing
par(mfrow = c(1, 1))
plot(sim2, y = "GCsympttests", mean.col = "black", qnts.col = "black", qnts = 0.5, ylab = "Tests/week")
plot(sim2, y = "CTsympttests", mean.col = "purple", qnts.col = "purple", qnts = 0.5, add = TRUE)
plot(sim2, y = "stisympttests", mean.col = "brown", qnts.col = "brown", qnts = 0.5, add = TRUE)
legend("topleft", lty = c(1, 1, 1), col = c("black", "purple", "brown"),
       legend = c("NG", "CT", "STI"))
title("STI Symptomatic Testing", outer = TRUE)

# Syphilis stage-specific prevalence (conditional on infection)
par(mfrow = c(1, 2))
plot(sim, y = "prev.stage.incub", mean.col = "green", qnts.col = "green",
     qnts = 0.5, qnts.alpha = 0.5, ylim = c(0, 0.4), ylab = "Proportion of all prevalent infections")
plot(sim, y = "prev.stage.prim", mean.col = "blue", qnts.col = "blue",
     qnts = 0.5, qnts.alpha = 0.5, add = TRUE)
plot(sim, y = "prev.stage.seco", mean.col = "red", qnts.col = "red",
     qnts = 0.5, qnts.alpha = 0.5, add = TRUE)
plot(sim, y = "prev.stage.earlat", mean.col = "orange", qnts.col = "orange",
     qnts = 0.5, qnts.alpha = 0.5, add = TRUE)
legend("topleft", lty = c(1,1,1,1), col = c("green", "blue", "red", "orange"),
       legend = c("Incub", "Prim", "Seco", "Early Latent"))

plot(sim, y = "prev.stage.latelat", mean.col = "purple", qnts.col = "purple",
     qnts = 0.5, ylab = "Proportion of all prevalent infections")
plot(sim, y = "prev.stage.tert", mean.col = "black", qnts.col = "black",
     qnts = 0.5, add = TRUE)
legend("topleft", lty = c(1,1), col = c("purple", "black"),
       legend = c("Late Latent", "Tertiary"))
title("Prevalence of Stage-Specific Syphilis", outer = TRUE)

# Check values over last time steps
df <- as.data.frame(x = sim, out = "vals")
sum(tail(df$num.newearlydiagsyph, 52))
sum(tail(df$num.newlatediagsyph, 52))

mean(tail(as.data.frame(sim2)$ir100.gc, 10))
mean(tail(as.data.frame(sim2)$ir100.ct, 10))
mean(tail(as.data.frame(sim2)$ir100.syph, 10))
mean(tail(as.data.frame(sim2)$i.prev, 10))
mean(tail(as.data.frame(sim2)$ir100.gc, 5))
mean(tail(as.data.frame(sim2)$ir100.ct, 5))
mean(tail(as.data.frame(sim2)$ir100.syph, 5))
mean(tail(as.data.frame(sim2)$i.prev, 5))

# Save as best-fitting ----------------------------------------------------
sim <- sim2
save(sim, file = "est/stimod.burnin.rda")
#save(sim, file = "est/stidiagmodchange.burnin.rda")
system("scp est/stimod.burnin.rda hyak:/gscratch/csde/sjenness/sti/est/")

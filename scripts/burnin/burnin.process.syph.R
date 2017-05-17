
## Process burn-in
library("EpiModelHPC")
library("EpiModelHIV")

# system("scp scripts/burnin/*.burn.[Rs]* hyak:/gscratch/csde/sjenness/sti")
# 
# # Examine output
# system("scp hyak:/gscratch/csde/sjenness/sti/data/*.rda data/")
# 
# list.files("data/")

load("data/sim.n100.rda")

df <- as.data.frame(sim)

par(mar = c(3,3,1,1), mgp = c(2,1,0))
plot(sim, y = "i.prev", ylim = c(0.1, 0.3), qnts = 0.5, mean.lwd = 1)
abline(h = 0.26, lty = 2)
text(x = 0, y = 0.28, round(mean(tail(df$i.prev, 100)), 3))

plot(sim, y = "ir100.gc", mean.smooth = FALSE, mean.lwd = 1, ylim = c(0, 10))
abline(h = 4.2, lty = 2)
text(0, 1, round(mean(tail(df$ir100.gc, 520)), 2))

plot(sim, y = "ir100.ct", mean.smooth = FALSE, mean.lwd = 1, ylim = c(0, 10))
abline(h = 6.6, lty = 2)
text(0, 1, round(mean(tail(df$ir100.ct, 520)), 2))


sim <- truncate_sim(sim, at = 2600)
par(mar = c(3,3,1,1), mgp = c(2,1,0))
plot(sim, y = "i.prev", ylim = c(0.1, 0.3), qnts = 0.5, mean.lwd = 1)
abline(h = 0.26, lty = 2)
text(x = 0, y = 0.28, round(mean(tail(df$i.prev, 100)), 3))

plot(sim, y = "ir100.gc", mean.smooth = TRUE, mean.lwd = 1, ylim = c(0, 10))
abline(h = 4.2, lty = 2)
text(0, 1, round(mean(tail(df$ir100.gc, 52)), 2))

plot(sim, y = "ir100.ct", mean.smooth = TRUE, mean.lwd = 1, ylim = c(0, 10))
abline(h = 6.6, lty = 2)
text(0, 1, round(mean(tail(df$ir100.ct, 52)), 2))

plot(sim, y = "ir100", mean.smooth = TRUE, mean.lwd = 1, ylim = c(0, 10))
abline(h = 3.8, lty = 2)

save(sim, file = "data/sim.n100.rda")

# Other Calibration ---------------------------------------------------

# Merge sim files
sim <- merge_simfiles(simno = 151, indir = "data/", ftype = "max")

# Create function for selecting sim closest to target
mean_sim <- function(sim, targets) {
    

  nsims <- sim$control$nsims

  # Initialize distance vector
  dist <- rep(NA, nsims)

  # Obtain statistics and perform multivariable Euclidean distance calculation
  for (i in 1:nsims) {

      # Create data frame to draw statistics from
      df <- as.data.frame(x = sim, out = "vals", sim = i)

      # Create a vector of statistics
      calib <- c(mean(tail(df$ir100.gc, 52)),
                 mean(tail(df$ir100.ct, 52)),
                 mean(tail(df$ir100, 52)),
                 mean(tail(df$i.prev, 10)),
                 mean(tail(df$prev.primsecosyph.hivpos, 10)),
                 mean(tail(df$prev.primsecosyph.hivneg, 10)),
                 mean(tail(df$prev.primsecosyph, 10)),
                 mean(tail(df$prev.hiv.primsecosyphpos, 10)),
                 mean(df$ir100.gc[2600] - df$ir100.gc[2595]),
                 mean(df$ir100.ct[2600] - df$ir100.ct[2595]),
                 mean(df$ir100[2600] - df$ir100[2595]),
                 mean(df$i.prev[2600] - df$i.prev[2595]),
                 mean(df$prev.primsecosyph[2600] - df$prev.primsecosyph[2595])
                )#,
                 #mean(tail(df$prev.earlysyph, 1)),
                 #mean(tail(df$prev.latesyph, 1)))

      wts <- c(3, 3, 3, 3, 2, 1, 3, 1, 1, 1, 1, 1, 1)

      # Iteratively calculate distance
      dist[i] <- sqrt(sum(((calib - targets)*wts)^2))
  }

  # Which sim minimizes distance
  meansim <- which.min(dist)
  return(meansim)
}

# Run function
mean_sim(sim, targets = c(4.2, 6.6, 3.8, 0.26, 0.103, 0.026, 0.046, 0.498, 0, 0, 0, 0, 0))

# Save burn-in file for FU sims
sim2 <- get_sims(sim, sims = 43)

# Check means
mean(tail(as.data.frame(sim2)$ir100.gc, 52))
mean(tail(as.data.frame(sim2)$ir100.ct, 52))
mean(tail(as.data.frame(sim2)$i.prev, 52))
mean(tail(as.data.frame(sim2)$ir100, 52))
mean(tail(as.data.frame(sim2)$prev.primsecosyph.hivpos, 52))
mean(tail(as.data.frame(sim2)$prev.primsecosyph.hivneg, 52))
mean(tail(as.data.frame(sim2)$prev.hiv.primsecosyphpos, 52))
mean(tail(as.data.frame(sim2)$prev.primsecosyph, 52))

mean(tail(as.data.frame(sim2)$ir100.syph, 52))
mean(tail(as.data.frame(sim2)$prev.syph.hivpos, 52))
mean(tail(as.data.frame(sim2)$prev.syph.hivneg, 52))
mean(tail(as.data.frame(sim2)$prev.earlysyph, 52))
mean(tail(as.data.frame(sim2)$prev.latesyph, 52))

mean(as.data.frame(sim2)$ir100.gc[2600] - as.data.frame(sim2)$ir100.gc[2595])
mean(as.data.frame(sim2)$ir100.ct[2600] - as.data.frame(sim2)$ir100.ct[2595])
mean(as.data.frame(sim2)$ir100[2600] - as.data.frame(sim2)$ir100[2595])
mean(as.data.frame(sim2)$i.prev[2600] - as.data.frame(sim2)$i.prev[2595])
mean(as.data.frame(sim2)$prev.primsecosyph[2600] - as.data.frame(sim2)$prev.primsecosyph[2595])

par(mfrow = c(2,2), oma = c(0,0,2,0))
plot(sim2, y = "ir100")
abline(h = 3.8, col = "red", lty = 2)
title("HIV Incidence")
plot(sim2, y = "ir100.gc")
abline(h = 4.2, col = "red", lty = 2)
title("GC Incidence")
plot(sim2, y = "ir100.ct")
abline(h = 6.6, col = "red", lty = 2)
title("CT Incidence")
plot(sim2, y = "ir100.syph")
abline(h = 0.9, col = "red", lty = 2)
title("Syph Incidence")

plot(sim2, y = "ir100", xlim = c(2200, 2600), ylim = c(3.0, 5.0))
abline(h = 3.8, col = "red", lty = 2)
title("HIV Incidence")
plot(sim2, y = "ir100.gc", xlim = c(2200, 2600), ylim = c(3.0, 5.0))
abline(h = 4.2, col = "red", lty = 2)
title("GC Incidence")
plot(sim2, y = "ir100.ct", xlim = c(2200, 2600), ylim = c(5.5, 7.5))
abline(h = 6.6, col = "red", lty = 2)
title("CT Incidence")


par(mfrow = c(2,2), oma = c(0,0,2,0))
plot(sim2, y = "prev.primsecosyph", ylab = "Prevalence")
title("P and S Syphilis Prevalence")
abline(h = 0.046, col = "red", lty = 2)
plot(sim2, y = "prev.ct", ylab = "Prevalence")
title("CT Prevalence")
plot(sim2, y = "prev.gc", ylab = "Prevalence")
title("GC Prevalence")
plot(sim2, y = "i.prev", ylim = c(0, 0.3), ylab = "Prevalence")
abline(h = 0.26, col = "red",  lty = 2)
title("HIV Prevalence")


sim <- sim2

save(sim, file = "est/stimod.burnin.rda")
# system("scp est/stimod.burnin.rda hyak:/gscratch/csde/sjenness/sti/est/")

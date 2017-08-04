## STI Test Guidelines Figure 4

rm(list = ls())
library("EpiModelHIV")
library("EpiModelHPC")
library("dplyr")
source("analysis/fx.R")

## Base STI treatment completion: n3054
## Varying Treatment completion for 40% lower risk and 0% higher-risk: n3175-n3184
tiff(filename = "analysis/Fig4a.tiff", height = 6, width = 11, units = "in", res = 250)
par(mfrow = c(1,2), mar = c(3,3,2,1.2), mgp = c(2,1,0))
sims <- c(3175:3184, 3054, 3000)
pal <- viridis::viridis(n = length(sims), option = "D")

for (i in seq_along(sims)) {
    fn <- list.files("data/followup/", pattern = as.character(sims[i]), full.names = TRUE)
    load(fn)
    plot(sim, y = "ir100.sti", add = i > 1, ylim = c(4, 8),
         mean.col = pal[i], qnts.col = pal[i], qnts.alpha = 0.3, qnts = 0,
         main = "STI Incidence by Treatment Probability\n 40% Ann / 0% HR Cov",
         xlab = "Week", ylab = "IR per 100 PYAR")
    abline(h = seq(4, 8, 0.2), lty = 2, col = "gray")
}
legend("bottomleft", legend = c("0%", "10%", "20%", "30%", "40%", "50%", "60%", "70%", "80%", "90%", "100%", "No testing"),
       col = pal, lwd = 3, cex = 0.85, bty = "n")

## Base STI treatment completion: n3014
## Varying Treatment completion for 40% higher-risk and 0% lower-risk: n3185-n3194
sims <- c(3185:3194, 3014, 3000)
pal <- viridis::viridis(n = length(sims), option = "D")
for (i in seq_along(sims)) {
    fn <- list.files("data/followup/", pattern = as.character(sims[i]), full.names = TRUE)
    load(fn)
    plot(sim, y = "ir100.sti", add = i > 1, ylim = c(4, 8),
         mean.col = pal[i], qnts.col = pal[i], qnts.alpha = 0.3, qnts = 0,
         main = "STI Incidence by Treatment Probability \n 40% HR / 0% Ann Cov",
         xlab = "Week", ylab = "IR per 100 PYAR")
    abline(h = seq(4, 8, 0.2), lty = 2, col = "gray")
}
legend("bottomleft", legend = c("0%", "10%", "20%", "30%", "40%", "50%", "60%", "70%", "80%", "90%", "100%", "No testing"),
       col = pal, lwd = 3, cex = 0.85, bty = "n")
dev.off()


##### HIV

## Base STI treatment completion: n3054
## Varying Treatment completion for 40% lower risk and 0% higher-risk: n3175-n3184
tiff(filename = "analysis/Fig4b.tiff", height = 6, width = 11, units = "in", res = 250)
par(mfrow = c(1,2), mar = c(3,3,2,1.2), mgp = c(2,1,0))
sims <- c(3175:3184, 3054, 3000)
pal <- viridis::viridis(n = length(sims), option = "D")

for (i in seq_along(sims)) {
    fn <- list.files("data/followup/", pattern = as.character(sims[i]), full.names = TRUE)
    load(fn)
    plot(sim, y = "ir100", add = i > 1, ylim = c(3, 4),
         mean.col = pal[i], qnts.col = pal[i], qnts.alpha = 0.3, qnts = 0,
         main = "HIV Incidence by STI Treatment Probability \n 40% Ann / 0% HR Cov",
         xlab = "Week", ylab = "IR per 100 PYAR")
    abline(h = seq(3, 4, 0.1), lty = 2, col = "gray")
}
legend("bottomleft", legend = c("0%", "10%", "20%", "30%", "40%", "50%", "60%", "70%", "80%", "90%", "100%", "No testing"), 
       col = pal, lwd = 3, cex = 0.85, bty = "n")

## Base STI treatment completion: n3014
## Varying Treatment completion for 40% higher-risk and 0% lower-risk: n3185-n3194
sims <- c(3185:3194, 3014, 3000)
pal <- viridis::viridis(n = length(sims), option = "D")
for (i in seq_along(sims)) {
    fn <- list.files("data/followup/", pattern = as.character(sims[i]), full.names = TRUE)
    load(fn)
    plot(sim, y = "ir100", add = i > 1, ylim = c(3, 4),
         mean.col = pal[i], qnts.col = pal[i], qnts.alpha = 0.3, qnts = 0,
         main = "HIV Incidence by STI Treatment Probability \n 40% HR / 0% Ann Cov",
         xlab = "Week", ylab = "IR per 100 PYAR")
    abline(h = seq(3, 4, 0.1), lty = 2, col = "gray")
}
legend("bottomleft", legend = c("0%", "10%", "20%", "30%", "40%", "50%", "60%", "70%", "80%", "90%", "100%", "No testing"), 
       col = pal, lwd = 3, cex = 0.85, bty = "n")
dev.off()

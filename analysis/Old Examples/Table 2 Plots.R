## STI Test Guidelines Plots for Table 2

rm(list = ls())
library("EpiModelHIV")
library("EpiModelHPC")
library("dplyr")
source("analysis/fx.R")

par(mfrow = c(1,2), mar = c(3,3,2,1.2), mgp = c(2,1,0))
sims <- c(3000, 3014, 3025, 3036, 3047, 3058, 3069, 3080, 3091, 3102, 3113, 3124)
pal <- viridis::viridis(n = length(sims), option = "D")

for (i in seq_along(sims)) {
    fn <- list.files("data/followup/", pattern = as.character(sims[i]), full.names = TRUE)
    load(fn)
    plot(sim, y = "ir100.sti", add = i > 1, ylim = c(0, 8),
         mean.col = pal[i], qnts.col = pal[i], qnts.alpha = 0.3, qnts = 0,
         main = "STI Incidence by Lower-Risk Coverage \n (40% Higher-risk coverage)",
         xlab = "Week", ylab = "IR per 100 PYAR")
    abline(h = seq(0, 8, 0.25), lty = 2, col = "gray")
}
legend("bottomleft", legend = c("0%", "10%", "20%", "30%", "40%", "50%", "60%", "70%", "80%", "90%", "100%"), 
       col = pal, lwd = 3, cex = 0.85, bty = "n")


sims <- c(3000, 3054:3064)
for (i in seq_along(sims)) {
    fn <- list.files("data/followup/", pattern = as.character(sims[i]), full.names = TRUE)
    load(fn)
    plot(sim, y = "ir100.sti", add = i > 1, ylim = c(0, 8),
         mean.col = pal[i], qnts.col = pal[i], qnts.alpha = 0.3, qnts = 0,
         main = "STI Incidence by Higher-Risk Coverage \n (40% Lower-risk coverage)",
         xlab = "Week", ylab = "IR per 100 PYAR")
    abline(h = seq(0, 8, 0.25), lty = 2, col = "gray")
}
legend("bottomleft", legend = c("0%", "10%", "20%", "30%", "40%", "50%", "60%", "70%", "80%", "90%", "100%"),
       col = pal, lwd = 3, cex = 0.85, bty = "n")
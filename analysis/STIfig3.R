## STI Test Guidelines Figure 3

rm(list = ls())
library("EpiModelHIV")
library("EpiModelHPC")
library("dplyr")
source("analysis/fx.R")

## Base STI lower-risk testing interval (364 days): n3054
## Varying STI lower-risk testing interval: n3131-n3141

tiff(filename = "analysis/Fig3a.tiff", height = 6, width = 11, units = "in", res = 250)
par(mfrow = c(1,2), mar = c(3,3,2,1.2), mgp = c(2,1,0))
sims <- c(3131:3141, 3054, 3000)
pal <- viridis::viridis(n = length(sims), option = "D")

for (i in seq_along(sims)) {
    fn <- list.files("data/followup/", pattern = as.character(sims[i]), full.names = TRUE)
    load(fn)
    plot(sim, y = "ir100.sti", add = i > 1, ylim = c(0, 8),
         mean.col = pal[i], qnts.col = pal[i], qnts.alpha = 0.3, qnts = 0.2,
         main = "STI Incidence by Lower-Risk STI Screening Interval \n 40% Ann / 0% HR Coverage",
         xlab = "Week", ylab = "IR per 100 PYAR")
    abline(h = seq(0, 8, 0.2), lty = 2, col = "gray")
}
legend("bottomleft", legend = c("28 days", "63 days", "91 days","119 days", "147 days", "182 days", 
                                "210 days", "238 days", "273 days", "301 days", "329 days", "364 days", "No testing"),
       col = pal, lwd = 3, cex = 0.85, bty = "n")


## Base STI higher-risk testing interval: n3014
## Varying STI higher-risk testing interval: n3153 - n3173
sims <- c(3153:3174, 3014, 3000)
pal <- viridis::viridis(n = length(sims), option = "D")
for (i in seq_along(sims)) {
    fn <- list.files("data/followup/", pattern = as.character(sims[i]), full.names = TRUE)
    load(fn)
    plot(sim, y = "ir100.sti", add = i > 1, ylim = c(0, 8),
         mean.col = pal[i], qnts.col = pal[i], qnts.alpha = 0.3, qnts = 0.2,
         main = "STI Incidence by Higher-Risk STI Screening Interval \n 40% HR / 0% Ann Cov",
         xlab = "Week", ylab = "IR per 100 PYAR")
    abline(h = seq(0, 8, 0.2), lty = 2, col = "gray")
}
legend("bottomleft", legend = c("28 days", "42 days", "49 days", "56 days", "63 days",
                                "70 days", "77 days", "84 days", "91 days","119 days", 
                                "126 days", "133 days", "140 days", "147 days", "154 days",
                                "161 days", "168 days", "175 days", "182 days", "No testing"), 
       col = pal, lwd = 3, cex = 0.85, bty = "n")
dev.off()


### HIV
tiff(filename = "analysis/Fig3b.tiff", height = 6, width = 11, units = "in", res = 250)
par(mfrow = c(1,2), mar = c(3,3,2,1.2), mgp = c(2,1,0))
sims <- c(3131:3141, 3054, 3000)
pal <- viridis::viridis(n = length(sims), option = "D")

for (i in seq_along(sims)) {
    fn <- list.files("data/followup/", pattern = as.character(sims[i]), full.names = TRUE)
    load(fn)
    plot(sim, y = "ir100", add = i > 1, ylim = c(0, 4),
         mean.col = pal[i], qnts.col = pal[i], qnts.alpha = 0.3, qnts = 0,
         main = "HIV Incidence by Lower-Risk STI Screening Interval \n 40% Ann / 0% HR Cov",
         xlab = "Week", ylab = "IR per 100 PYAR")
    abline(h = seq(0, 4, 0.2), lty = 2, col = "gray")
}
legend("bottomleft", legend = c("28 days", "63 days", "91 days","119 days", "147 days", "182 days", 
                                "210 days", "238 days", "273 days", "301 days", "329 days", "364 days", "No testing"),
       col = pal, lwd = 3, cex = 0.85, bty = "n")


## Base STI higher-risk testing interval: n3014
## Varying STI higher-risk testing interval: n3153 - n3173
sims <- c(3153:3174, 3014, 3000)
pal <- viridis::viridis(n = length(sims), option = "D")
for (i in seq_along(sims)) {
    fn <- list.files("data/followup/", pattern = as.character(sims[i]), full.names = TRUE)
    load(fn)
    plot(sim, y = "ir100", add = i > 1, ylim = c(0, 4),
         mean.col = pal[i], qnts.col = pal[i], qnts.alpha = 0.3, qnts = 0,
         main = "HIV Incidence by Higher-Risk STI Screening Interval \n 40% HR / 0% Ann Cov",
         xlab = "Week", ylab = "IR per 100 PYAR")
         abline(h = seq(0, 4, 0.2), lty = 2, col = "gray")
}

legend("bottomleft", legend = c("28 days", "42 days", "49 days", "56 days", "63 days",
                                "70 days", "77 days", "84 days", "91 days","119 days", 
                                "126 days", "133 days", "140 days", "147 days", "154 days",
                                "161 days", "168 days", "175 days", "182 days", "No testing"), 
       col = pal, lwd = 3, cex = 0.85, bty = "n")
dev.off()

## STI Testing Guidelines Figure 2 - Boxplots
## Partner Cutoff --------------------------------------------------------------

rm(list = ls())
suppressMessages(library("EpiModelHIV"))
library("EpiModelHPC")
library("dplyr")
source("analysis/fx.R")

# Box Plots by Indications
library('wesanderson')
par(mfrow = c(1, 2), mar = c(3,3,2.5,1), mgp = c(2,1,0))
pal <- wesanderson::wes_palette("Moonrise1", n = 9, type = "continuous")

# Baseline
load("data/followup/sim.n3003.rda")
sim.base <- sim
mn.base <- as.data.frame(sim.base)
ir.base <- (sum(mn.base$incid)/sum((1 - mn.base$i.prev) * mn.base$num)) * 52 * 1e5
ir.gc.base <- (sum(mn.base$incid.gc)/sum((1 - mn.base$prev.gc) * mn.base$num)) * 52 * 1e5
ir.ct.base <- (sum(mn.base$incid.ct)/sum((1 - mn.base$prev.ct) * mn.base$num)) * 52 * 1e5
ir.syph.base <- (sum(mn.base$incid.syph)/sum((1 - mn.base$prev.syph) * mn.base$num)) * 52 * 1e5
ir.sti.base <- (sum(mn.base$incid.sti)/sum((1 - mn.base$prev.sti) * mn.base$num)) * 52 * 1e5
incid.base <- sum(mn.base$incid)
incid.gc.base <- sum(mn.base$incid.gc)
incid.ct.base <- sum(mn.base$incid.ct)
incid.syph.base <- sum(mn.base$incid.syph)
incid.sti.base <- sum(mn.base$incid.sti)

# Sims 3142-3152 for Number of partners - 20% higher-risk coverage
sims <- c(3003, 3490:3498)
df.hiv.pia <- data.frame(rep(NA, 256))
df.hivonly.nnt <- data.frame(rep(NA, 256))
df.hiv.nnt <- data.frame(rep(NA, 256))
df.gc.pia <- data.frame(rep(NA, 256))
df.gc.nnt <- data.frame(rep(NA, 256))
df.ct.pia <- data.frame(rep(NA, 256))
df.ct.nnt <- data.frame(rep(NA, 256))
df.syph.pia <- data.frame(rep(NA, 256))
df.syph.nnt <- data.frame(rep(NA, 256))
df.sti.pia <- data.frame(rep(NA, 256))
df.sti.nnt <- data.frame(rep(NA, 256))

for (i in seq_along(sims)) {

  load(list.files("data/followup/", pattern = as.character(sims[i]), full.names = TRUE))
  mn <- as.data.frame(sim)
  ir <- (colSums(sim$epi$incid, na.rm = TRUE)) /
    sum((1 - mn$i.prev)  * mn$num) * 52 * 1e5
  ir.gc <- (colSums(sim$epi$incid.gc, na.rm = TRUE)) /
    sum((1 - mn$prev.gc)  * mn$num) * 52 * 1e5
  ir.ct <- (colSums(sim$epi$incid.ct, na.rm = TRUE)) /
    sum((1 - mn$prev.ct)  * mn$num) * 52 * 1e5
  ir.syph <- (colSums(sim$epi$incid.syph, na.rm = TRUE)) /
    sum((1 - mn$prev.syph)  * mn$num) * 52 * 1e5
  ir.sti <- (colSums(sim$epi$incid.sti, na.rm = TRUE)) /
    sum((1 - mn$prev.sti)  * mn$num) * 52 * 1e5

  vec.hiv.nia <- round(ir.base - unname(ir), 1)
  df.hiv.pia[, i] <- vec.hiv.nia / ir.base

  vec.gc.nia <- round(ir.gc.base - unname(ir.gc), 1)
  df.gc.pia[, i] <- vec.gc.nia / ir.gc.base

  vec.ct.nia <- round(ir.ct.base - unname(ir.ct), 1)
  df.ct.pia[, i] <- vec.ct.nia / ir.ct.base

  vec.syph.nia <- round(ir.syph.base - unname(ir.syph), 1)
  df.syph.pia[, i] <- vec.syph.nia / ir.syph.base

  vec.sti.nia <- round(ir.sti.base - unname(ir.sti), 1)
  df.sti.pia[, i] <- vec.sti.nia / ir.sti.base

  hiv.tests <- unname(colSums(tail(sim$epi$hivtests.nprep)))
  gc.asympt.tests <- unname(colSums(tail(sim$epi$GCasympttests)))
  ct.asympt.tests <- unname(colSums(tail(sim$epi$CTasympttests)))
  syph.asympt.tests <- unname(colSums(tail(sim$epi$syphasympttests)))

  #HIV could be HIV tests or total STI tests
  df.hivonly.nnt[, i] <- (hiv.tests) / (incid.base - unname(colSums(sim$epi$incid)))
  df.hiv.nnt[, i] <- (gc.asympt.tests + gc.asympt.tests + syph.asympt.tests) / (incid.base - unname(colSums(sim$epi$incid)))
  df.gc.nnt[, i] <- gc.asympt.tests / (incid.gc.base - unname(colSums(sim$epi$incid.gc)))
  df.ct.nnt[, i] <- ct.asympt.tests / (incid.ct.base - unname(colSums(sim$epi$incid.ct)))
  df.syph.nnt[, i] <- syph.asympt.tests / (incid.syph.base - unname(colSums(sim$epi$incid.syph)))
  df.sti.nnt[, i] <- syph.asympt.tests / (incid.syph.base - unname(colSums(sim$epi$incid.syph)))

}
names(df.hiv.pia) <- names(df.gc.pia) <- names(df.gc.nnt) <- names(df.ct.pia) <- names(df.sti.pia) <-
  names(df.hivonly.nnt) <- names(df.hiv.nnt) <- names(df.ct.nnt) <- names(df.syph.pia) <- names(df.syph.nnt) <- names(df.sti.nnt) <-
  c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10")

head(df.hiv.pia)
head(df.gc.pia)
head(df.ct.pia)
head(df.syph.pia)
head(df.sti.pia)
head(df.hiv.nnt)
head(df.gc.nnt)
head(df.ct.nnt)
head(df.syph.nnt)
head(df.sti.nnt)

pal <- wes_palette("Zissou")[c(1, 5)]
# HIV
# tiff(filename = "analysis/Fig2a.tiff", height = 4, width = 10, units = "in", res = 250)
# par(mfrow = c(2, 1), mar = c(4,4,2.5,1), oma = c(0, 0, 2, 0), mgp = c(3, 0.75, 0))
# boxplot(df.hiv.pia, outline = FALSE, medlwd = 1.1,
#         col = c(rep(pal[1], 8), rep(pal[2], 4)),  ylim = c(0, 1),
#         main = "PIA by Partner Cutoff", las = 2,
#         xlab = "Partner Cutoff", ylab = "Proportion HIV Infections Averted",
#         cex.axis = 0.7)
#
# boxplot(df.hiv.nnt, outline = FALSE, medlwd = 1.1,
#         col = c(rep(pal[1], 8), rep(pal[2], 4)),
#         main = "NNT by Partner Cutoff", las = 2,
#         xlab = "Partner Cutoff", ylab = "HIV Number Needed to Treat",
#         cex.axis = 0.7)
# title("20% HR Cov (6 months), 0% Ann Cov", outer = TRUE)
# dev.off()

# GC
# tiff(filename = "analysis/Fig2b.tiff", height = 4, width = 8, units = "in", res = 250)
par(mfrow = c(2, 2), mar = c(4,4,2.5,1), oma = c(0, 0, 2, 0), mgp = c(3, 0.75, 0))
boxplot(df.gc.pia, outline = FALSE, medlwd = 1.1,
        col = c(rep(pal[1], 8), rep(pal[2], 4)), ylim = c(0, 1),
        main = "PIA by Partner Cutoff", las = 2,
        xlab = "Partner Cutoff", ylab = "Proportion NG Infections Averted",
        cex.axis = 0.7)
# boxplot(df.gc.nnt, outline = FALSE, medlwd = 1.1,
#         col = c(rep(pal[1], 8), rep(pal[2], 4)),
#         main = "NNT by Partner Cutoff", las = 2,
#         xlab = "Partner Cutoff", ylab = "NG Number Needed to Treat",
#         cex.axis = 0.7)
# title("20% HR Cov (6 months), 10% Ann Cov", outer = TRUE)
# dev.off()

# CT
# tiff(filename = "analysis/Fig2c.tiff", height = 4, width = 8, units = "in", res = 250)
# par(mfrow = c(1, 2), mar = c(4,4,2.5,1), oma = c(0, 0, 2, 0), mgp = c(3, 0.75, 0))
boxplot(df.ct.pia, outline = FALSE, medlwd = 1.1,
        col = c(rep(pal[1], 8), rep(pal[2], 4)), ylim = c(0, 1),
        main = "PIA by Partner Cutoff", las = 2,
        xlab = "Partner Cutoff", ylab = "Proportion CT Infections Averted",
        cex.axis = 0.7)
# boxplot(df.ct.nnt, outline = FALSE, medlwd = 1.1,
#         col = c(rep(pal[1], 8), rep(pal[2], 4)),
#         main = "NNT by Partner Cutoff",las = 2,
#         xlab = "Partner Cutoff", ylab = "CT Number Needed to Treat",
#         cex.axis = 0.7)
# title("20% HR Cov (6 months), 10% Ann Cov", outer = TRUE)
# dev.off()

# Syph
# tiff(filename = "analysis/Fig2d.tiff", height = 4, width = 8, units = "in", res = 250)
# par(mfrow = c(1, 2), mar = c(4,4,2.5,1), oma = c(0, 0, 2, 0), mgp = c(3, 0.75, 0))
boxplot(df.syph.pia, outline = FALSE, medlwd = 1.1,
        col = c(rep(pal[1], 8), rep(pal[2], 4)), ylim = c(0, 1),
        main = "PIA by Partner Cutoff", las = 2,
        xlab = "Partner Cutoff", ylab = "Proportion Syph Infections Averted",
        cex.axis = 0.7)
# boxplot(df.syph.nnt, outline = FALSE, medlwd = 1.1,
#         col = c(rep(pal[1], 8), rep(pal[2], 4)),
#         main = "NNT by Partner Cutoff",las = 2,
#         xlab = "Partner Cutoff", ylab = "Syph Number Needed to Treat",
#         cex.axis = 0.7)
# title("20% HR Cov (6 months), 10% Ann Cov", outer = TRUE)
# dev.off()

# Combined STI
boxplot(df.syph.pia, outline = FALSE, medlwd = 1.1,
col = c(rep(pal[1], 8), rep(pal[2], 4)), ylim = c(0, 1),
main = "PIA by Partner Cutoff", las = 2,
xlab = "Partner Cutoff", ylab = "Proportion STI Infections Averted",
cex.axis = 0.7)
title("20% HR Cov (6 months), 10% Ann Cov", outer = TRUE)
dev.off()


## Screening Intervals - Lower Risk --------------------------------------------
## STI Testing Guidelines Figure 2 - Version 2 - Boxplots

rm(list = ls())
suppressMessages(library("EpiModelHIV"))
library("EpiModelHPC")
library("dplyr")
source("analysis/fx.R")

# Box Plots by Indications
library('wesanderson')
par(mfrow = c(1, 2), mar = c(3,3,2.5,1), mgp = c(2,1,0))
pal <- wesanderson::wes_palette("Moonrise1", n = 9, type = "continuous")

# Baseline
load("data/followup/sim.n3003.rda")
sim.base <- sim
mn.base <- as.data.frame(sim.base)
ir.base <- (sum(mn.base$incid)/sum((1 - mn.base$i.prev) * mn.base$num)) * 52 * 1e5
ir.gc.base <- (sum(mn.base$incid.gc)/sum((1 - mn.base$prev.gc) * mn.base$num)) * 52 * 1e5
ir.ct.base <- (sum(mn.base$incid.ct)/sum((1 - mn.base$prev.ct) * mn.base$num)) * 52 * 1e5
ir.syph.base <- (sum(mn.base$incid.syph)/sum((1 - mn.base$prev.syph) * mn.base$num)) * 52 * 1e5
ir.sti.base <- (sum(mn.base$incid.sti)/sum((1 - mn.base$prev.sti) * mn.base$num)) * 52 * 1e5
incid.base <- sum(mn.base$incid)
incid.gc.base <- sum(mn.base$incid.gc)
incid.ct.base <- sum(mn.base$incid.ct)
incid.syph.base <- sum(mn.base$incid.syph)
incid.sti.base <- sum(mn.base$incid.sti)

sims <- c(3442, 3443, 3003, 3444, 3445)
df.hiv.pia <- data.frame(rep(NA, 256))
df.hivonly.nnt <- data.frame(rep(NA, 256))
df.hiv.nnt <- data.frame(rep(NA, 256))
df.gc.pia <- data.frame(rep(NA, 256))
df.gc.nnt <- data.frame(rep(NA, 256))
df.ct.pia <- data.frame(rep(NA, 256))
df.ct.nnt <- data.frame(rep(NA, 256))
df.syph.pia <- data.frame(rep(NA, 256))
df.syph.nnt <- data.frame(rep(NA, 256))
df.sti.pia <- data.frame(rep(NA, 256))
df.sti.nnt <- data.frame(rep(NA, 256))

for (i in seq_along(sims)) {

  load(list.files("data/followup/", pattern = as.character(sims[i]), full.names = TRUE))
  mn <- as.data.frame(sim)
  ir <- (colSums(sim$epi$incid, na.rm = TRUE)) /
    sum((1 - mn$i.prev)  * mn$num) * 52 * 1e5
  ir.gc <- (colSums(sim$epi$incid.gc, na.rm = TRUE)) /
    sum((1 - mn$prev.gc)  * mn$num) * 52 * 1e5
  ir.ct <- (colSums(sim$epi$incid.ct, na.rm = TRUE)) /
    sum((1 - mn$prev.ct)  * mn$num) * 52 * 1e5
  ir.syph <- (colSums(sim$epi$incid.syph, na.rm = TRUE)) /
    sum((1 - mn$prev.syph)  * mn$num) * 52 * 1e5
  ir.sti <- (colSums(sim$epi$incid.sti, na.rm = TRUE)) /
    sum((1 - mn$prev.sti)  * mn$num) * 52 * 1e5

  vec.hiv.nia <- round(ir.base - unname(ir), 1)
  df.hiv.pia[, i] <- vec.hiv.nia / ir.base

  vec.gc.nia <- round(ir.gc.base - unname(ir.gc), 1)
  df.gc.pia[, i] <- vec.gc.nia / ir.gc.base

  vec.ct.nia <- round(ir.ct.base - unname(ir.ct), 1)
  df.ct.pia[, i] <- vec.ct.nia / ir.ct.base

  vec.syph.nia <- round(ir.syph.base - unname(ir.syph), 1)
  df.syph.pia[, i] <- vec.syph.nia / ir.syph.base

  vec.sti.nia <- round(ir.sti.base - unname(ir.sti), 1)
  df.sti.pia[, i] <- vec.sti.nia / ir.sti.base

  hiv.tests <- unname(colSums(tail(sim$epi$hivtests.nprep)))
  gc.asympt.tests <- unname(colSums(tail(sim$epi$GCasympttests)))
  ct.asympt.tests <- unname(colSums(tail(sim$epi$CTasympttests)))
  syph.asympt.tests <- unname(colSums(tail(sim$epi$syphasympttests)))

  #HIV could be HIV tests or total STI tests
  df.hivonly.nnt[, i] <- (hiv.tests) / (incid.base - unname(colSums(sim$epi$incid)))
  df.hiv.nnt[, i] <- (gc.asympt.tests + gc.asympt.tests + syph.asympt.tests) / (incid.base - unname(colSums(sim$epi$incid)))
  df.gc.nnt[, i] <- gc.asympt.tests / (incid.gc.base - unname(colSums(sim$epi$incid.gc)))
  df.ct.nnt[, i] <- ct.asympt.tests / (incid.ct.base - unname(colSums(sim$epi$incid.ct)))
  df.syph.nnt[, i] <- syph.asympt.tests / (incid.syph.base - unname(colSums(sim$epi$incid.syph)))
  df.sti.nnt[, i] <- syph.asympt.tests / (incid.syph.base - unname(colSums(sim$epi$incid.syph)))

}
names(df.hiv.pia) <- names(df.gc.pia) <- names(df.gc.nnt) <- names(df.ct.pia) <- names(df.sti.pia) <-
  names(df.hivonly.nnt) <- names(df.hiv.nnt) <- names(df.ct.nnt) <- names(df.syph.pia) <- names(df.syph.nnt) <- names(df.sti.nnt) <-
  c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10")

head(df.hiv.pia)
head(df.gc.pia)
head(df.ct.pia)
head(df.syph.pia)
head(df.sti.pia)
head(df.hiv.nnt)
head(df.gc.nnt)
head(df.ct.nnt)
head(df.syph.nnt)
head(df.sti.nnt)

pal <- wes_palette("Zissou")[c(1, 5)]
# HIV
# tiff(filename = "analysis/Fig2a.tiff", height = 4, width = 10, units = "in", res = 250)
# par(mfrow = c(2, 1), mar = c(4,4,2.5,1), oma = c(0, 0, 2, 0), mgp = c(3, 0.75, 0))
# boxplot(df.hiv.pia, outline = FALSE, medlwd = 1.1,
#         col = c(rep(pal[1], 8), rep(pal[2], 4)),  ylim = c(0, 1),
#         main = "PIA by Partner Cutoff", las = 2,
#         xlab = "Partner Cutoff", ylab = "Proportion HIV Infections Averted",
#         cex.axis = 0.7)
#
# boxplot(df.hiv.nnt, outline = FALSE, medlwd = 1.1,
#         col = c(rep(pal[1], 8), rep(pal[2], 4)),
#         main = "NNT by Partner Cutoff", las = 2,
#         xlab = "Partner Cutoff", ylab = "HIV Number Needed to Treat",
#         cex.axis = 0.7)
# title("20% HR Cov (6 months), 10% Ann Cov", outer = TRUE)
# dev.off()

# GC
# tiff(filename = "analysis/Fig2b.tiff", height = 4, width = 8, units = "in", res = 250)
par(mfrow = c(2, 2), mar = c(4,4,2.5,1), oma = c(0, 0, 2, 0), mgp = c(3, 0.75, 0))
boxplot(df.gc.pia, outline = FALSE, medlwd = 1.1,
        col = c(rep(pal[1], 8), rep(pal[2], 4)), ylim = c(0, 1),
        main = "PIA by Partner Cutoff", las = 2,
        xlab = "Partner Cutoff", ylab = "Proportion NG Infections Averted",
        cex.axis = 0.7)
# boxplot(df.gc.nnt, outline = FALSE, medlwd = 1.1,
#         col = c(rep(pal[1], 8), rep(pal[2], 4)),
#         main = "NNT by Partner Cutoff", las = 2,
#         xlab = "Partner Cutoff", ylab = "NG Number Needed to Treat",
#         cex.axis = 0.7)
# title("20% HR Cov (6 months), 10% Ann Cov", outer = TRUE)
# dev.off()

# CT
# tiff(filename = "analysis/Fig2c.tiff", height = 4, width = 8, units = "in", res = 250)
# par(mfrow = c(1, 2), mar = c(4,4,2.5,1), oma = c(0, 0, 2, 0), mgp = c(3, 0.75, 0))
boxplot(df.ct.pia, outline = FALSE, medlwd = 1.1,
        col = c(rep(pal[1], 8), rep(pal[2], 4)), ylim = c(0, 1),
        main = "PIA by Partner Cutoff", las = 2,
        xlab = "Partner Cutoff", ylab = "Proportion CT Infections Averted",
        cex.axis = 0.7)
# boxplot(df.ct.nnt, outline = FALSE, medlwd = 1.1,
#         col = c(rep(pal[1], 8), rep(pal[2], 4)),
#         main = "NNT by Partner Cutoff",las = 2,
#         xlab = "Partner Cutoff", ylab = "CT Number Needed to Treat",
#         cex.axis = 0.7)
# title("20% HR Cov (6 months), 10% Ann Cov", outer = TRUE)
# dev.off()

# Syph
# tiff(filename = "analysis/Fig2d.tiff", height = 4, width = 8, units = "in", res = 250)
# par(mfrow = c(1, 2), mar = c(4,4,2.5,1), oma = c(0, 0, 2, 0), mgp = c(3, 0.75, 0))
boxplot(df.syph.pia, outline = FALSE, medlwd = 1.1,
        col = c(rep(pal[1], 8), rep(pal[2], 4)), ylim = c(0, 1),
        main = "PIA by Partner Cutoff", las = 2,
        xlab = "Partner Cutoff", ylab = "Proportion Syph Infections Averted",
        cex.axis = 0.7)
# boxplot(df.syph.nnt, outline = FALSE, medlwd = 1.1,
#         col = c(rep(pal[1], 8), rep(pal[2], 4)),
#         main = "NNT by Partner Cutoff",las = 2,
#         xlab = "Partner Cutoff", ylab = "Syph Number Needed to Treat",
#         cex.axis = 0.7)
# title("20% HR Cov (6 months), 10% Ann Cov", outer = TRUE)
# dev.off()

# Combined STI
boxplot(df.syph.pia, outline = FALSE, medlwd = 1.1,
        col = c(rep(pal[1], 8), rep(pal[2], 4)), ylim = c(0, 1),
        main = "PIA by Partner Cutoff", las = 2,
        xlab = "Partner Cutoff", ylab = "Proportion STI Infections Averted",
        cex.axis = 0.7)
title("20% HR Cov (6 months), 10% Ann Cov", outer = TRUE)
dev.off()



## Screening Intervals - Higher-Risk -------------------------------------------
## STI Testing Guidelines Figure 2 - Version 2 - Boxplots

rm(list = ls())
suppressMessages(library("EpiModelHIV"))
library("EpiModelHPC")
library("dplyr")
source("analysis/fx.R")

# Box Plots by Indications
library('wesanderson')
par(mfrow = c(1, 2), mar = c(3,3,2.5,1), mgp = c(2,1,0))
pal <- wesanderson::wes_palette("Moonrise1", n = 9, type = "continuous")

# Baseline
load("data/followup/sim.n3003.rda")
sim.base <- sim
mn.base <- as.data.frame(sim.base)
ir.base <- (sum(mn.base$incid)/sum((1 - mn.base$i.prev) * mn.base$num)) * 52 * 1e5
ir.gc.base <- (sum(mn.base$incid.gc)/sum((1 - mn.base$prev.gc) * mn.base$num)) * 52 * 1e5
ir.ct.base <- (sum(mn.base$incid.ct)/sum((1 - mn.base$prev.ct) * mn.base$num)) * 52 * 1e5
ir.syph.base <- (sum(mn.base$incid.syph)/sum((1 - mn.base$prev.syph) * mn.base$num)) * 52 * 1e5
ir.sti.base <- (sum(mn.base$incid.sti)/sum((1 - mn.base$prev.sti) * mn.base$num)) * 52 * 1e5
incid.base <- sum(mn.base$incid)
incid.gc.base <- sum(mn.base$incid.gc)
incid.ct.base <- sum(mn.base$incid.ct)
incid.syph.base <- sum(mn.base$incid.syph)
incid.sti.base <- sum(mn.base$incid.sti)

sims <- c(3446, 3447, 3003, 3448, 3449)
df.hiv.pia <- data.frame(rep(NA, 256))
df.hivonly.nnt <- data.frame(rep(NA, 256))
df.hiv.nnt <- data.frame(rep(NA, 256))
df.gc.pia <- data.frame(rep(NA, 256))
df.gc.nnt <- data.frame(rep(NA, 256))
df.ct.pia <- data.frame(rep(NA, 256))
df.ct.nnt <- data.frame(rep(NA, 256))
df.syph.pia <- data.frame(rep(NA, 256))
df.syph.nnt <- data.frame(rep(NA, 256))
df.sti.pia <- data.frame(rep(NA, 256))
df.sti.nnt <- data.frame(rep(NA, 256))

for (i in seq_along(sims)) {

  load(list.files("data/followup/", pattern = as.character(sims[i]), full.names = TRUE))
  mn <- as.data.frame(sim)
  ir <- (colSums(sim$epi$incid, na.rm = TRUE)) /
    sum((1 - mn$i.prev)  * mn$num) * 52 * 1e5
  ir.gc <- (colSums(sim$epi$incid.gc, na.rm = TRUE)) /
    sum((1 - mn$prev.gc)  * mn$num) * 52 * 1e5
  ir.ct <- (colSums(sim$epi$incid.ct, na.rm = TRUE)) /
    sum((1 - mn$prev.ct)  * mn$num) * 52 * 1e5
  ir.syph <- (colSums(sim$epi$incid.syph, na.rm = TRUE)) /
    sum((1 - mn$prev.syph)  * mn$num) * 52 * 1e5
  ir.sti <- (colSums(sim$epi$incid.sti, na.rm = TRUE)) /
    sum((1 - mn$prev.sti)  * mn$num) * 52 * 1e5

  vec.hiv.nia <- round(ir.base - unname(ir), 1)
  df.hiv.pia[, i] <- vec.hiv.nia / ir.base

  vec.gc.nia <- round(ir.gc.base - unname(ir.gc), 1)
  df.gc.pia[, i] <- vec.gc.nia / ir.gc.base

  vec.ct.nia <- round(ir.ct.base - unname(ir.ct), 1)
  df.ct.pia[, i] <- vec.ct.nia / ir.ct.base

  vec.syph.nia <- round(ir.syph.base - unname(ir.syph), 1)
  df.syph.pia[, i] <- vec.syph.nia / ir.syph.base

  vec.sti.nia <- round(ir.sti.base - unname(ir.sti), 1)
  df.sti.pia[, i] <- vec.sti.nia / ir.sti.base

  hiv.tests <- unname(colSums(tail(sim$epi$hivtests.nprep)))
  gc.asympt.tests <- unname(colSums(tail(sim$epi$GCasympttests)))
  ct.asympt.tests <- unname(colSums(tail(sim$epi$CTasympttests)))
  syph.asympt.tests <- unname(colSums(tail(sim$epi$syphasympttests)))

  #HIV could be HIV tests or total STI tests
  df.hivonly.nnt[, i] <- (hiv.tests) / (incid.base - unname(colSums(sim$epi$incid)))
  df.hiv.nnt[, i] <- (gc.asympt.tests + gc.asympt.tests + syph.asympt.tests) / (incid.base - unname(colSums(sim$epi$incid)))
  df.gc.nnt[, i] <- gc.asympt.tests / (incid.gc.base - unname(colSums(sim$epi$incid.gc)))
  df.ct.nnt[, i] <- ct.asympt.tests / (incid.ct.base - unname(colSums(sim$epi$incid.ct)))
  df.syph.nnt[, i] <- syph.asympt.tests / (incid.syph.base - unname(colSums(sim$epi$incid.syph)))
  df.sti.nnt[, i] <- syph.asympt.tests / (incid.syph.base - unname(colSums(sim$epi$incid.syph)))

}
names(df.hiv.pia) <- names(df.gc.pia) <- names(df.gc.nnt) <- names(df.ct.pia) <- names(df.sti.pia) <-
  names(df.hivonly.nnt) <- names(df.hiv.nnt) <- names(df.ct.nnt) <- names(df.syph.pia) <- names(df.syph.nnt) <- names(df.sti.nnt) <-
  c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10")

head(df.hiv.pia)
head(df.gc.pia)
head(df.ct.pia)
head(df.syph.pia)
head(df.sti.pia)
head(df.hiv.nnt)
head(df.gc.nnt)
head(df.ct.nnt)
head(df.syph.nnt)
head(df.sti.nnt)

pal <- wes_palette("Zissou")[c(1, 5)]
# HIV
# tiff(filename = "analysis/Fig2a.tiff", height = 4, width = 10, units = "in", res = 250)
# par(mfrow = c(2, 1), mar = c(4,4,2.5,1), oma = c(0, 0, 2, 0), mgp = c(3, 0.75, 0))
# boxplot(df.hiv.pia, outline = FALSE, medlwd = 1.1,
#         col = c(rep(pal[1], 8), rep(pal[2], 4)),  ylim = c(0, 1),
#         main = "PIA by Partner Cutoff", las = 2,
#         xlab = "Partner Cutoff", ylab = "Proportion HIV Infections Averted",
#         cex.axis = 0.7)
#
# boxplot(df.hiv.nnt, outline = FALSE, medlwd = 1.1,
#         col = c(rep(pal[1], 8), rep(pal[2], 4)),
#         main = "NNT by Partner Cutoff", las = 2,
#         xlab = "Partner Cutoff", ylab = "HIV Number Needed to Treat",
#         cex.axis = 0.7)
# title("20% HR Cov (6 months), 0% Ann Cov", outer = TRUE)
# dev.off()

# GC
# tiff(filename = "analysis/Fig2b.tiff", height = 4, width = 8, units = "in", res = 250)
par(mfrow = c(2, 2), mar = c(4,4,2.5,1), oma = c(0, 0, 2, 0), mgp = c(3, 0.75, 0))
boxplot(df.gc.pia, outline = FALSE, medlwd = 1.1,
        col = c(rep(pal[1], 8), rep(pal[2], 4)), ylim = c(0, 1),
        main = "PIA by Partner Cutoff", las = 2,
        xlab = "Partner Cutoff", ylab = "Proportion NG Infections Averted",
        cex.axis = 0.7)
# boxplot(df.gc.nnt, outline = FALSE, medlwd = 1.1,
#         col = c(rep(pal[1], 8), rep(pal[2], 4)),
#         main = "NNT by Partner Cutoff", las = 2,
#         xlab = "Partner Cutoff", ylab = "NG Number Needed to Treat",
#         cex.axis = 0.7)
# title("20% HR Cov (6 months), 10% Ann Cov", outer = TRUE)
# dev.off()

# CT
# tiff(filename = "analysis/Fig2c.tiff", height = 4, width = 8, units = "in", res = 250)
# par(mfrow = c(1, 2), mar = c(4,4,2.5,1), oma = c(0, 0, 2, 0), mgp = c(3, 0.75, 0))
boxplot(df.ct.pia, outline = FALSE, medlwd = 1.1,
        col = c(rep(pal[1], 8), rep(pal[2], 4)), ylim = c(0, 1),
        main = "PIA by Partner Cutoff", las = 2,
        xlab = "Partner Cutoff", ylab = "Proportion CT Infections Averted",
        cex.axis = 0.7)
# boxplot(df.ct.nnt, outline = FALSE, medlwd = 1.1,
#         col = c(rep(pal[1], 8), rep(pal[2], 4)),
#         main = "NNT by Partner Cutoff",las = 2,
#         xlab = "Partner Cutoff", ylab = "CT Number Needed to Treat",
#         cex.axis = 0.7)
# title("20% HR Cov (6 months), 10% Ann Cov", outer = TRUE)
# dev.off()

# Syph
# tiff(filename = "analysis/Fig2d.tiff", height = 4, width = 8, units = "in", res = 250)
# par(mfrow = c(1, 2), mar = c(4,4,2.5,1), oma = c(0, 0, 2, 0), mgp = c(3, 0.75, 0))
boxplot(df.syph.pia, outline = FALSE, medlwd = 1.1,
        col = c(rep(pal[1], 8), rep(pal[2], 4)), ylim = c(0, 1),
        main = "PIA by Partner Cutoff", las = 2,
        xlab = "Partner Cutoff", ylab = "Proportion Syph Infections Averted",
        cex.axis = 0.7)
# boxplot(df.syph.nnt, outline = FALSE, medlwd = 1.1,
#         col = c(rep(pal[1], 8), rep(pal[2], 4)),
#         main = "NNT by Partner Cutoff",las = 2,
#         xlab = "Partner Cutoff", ylab = "Syph Number Needed to Treat",
#         cex.axis = 0.7)
# title("20% HR Cov (6 months), 10% Ann Cov", outer = TRUE)
# dev.off()

# Combined STI
boxplot(df.syph.pia, outline = FALSE, medlwd = 1.1,
        col = c(rep(pal[1], 8), rep(pal[2], 4)), ylim = c(0, 1),
        main = "PIA by Partner Cutoff", las = 2,
        xlab = "Partner Cutoff", ylab = "Proportion STI Infections Averted",
        cex.axis = 0.7)
title("20% HR Cov (6 months), 10% Ann Cov", outer = TRUE)
dev.off()

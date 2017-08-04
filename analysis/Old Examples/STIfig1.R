# Figure 1 ----------------------------------------------------------------
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
load("data/followup/sim.n3000.rda")
sim.base <- sim
mn.base <- as.data.frame(sim.base)
ir.base <- (sum(mn.base$incid)/sum((1 - mn.base$i.prev) * mn.base$num)) * 52 * 1e5
ir.gc.base <- (sum(mn.base$incid.gc)/sum((1 - mn.base$prev.gc) * mn.base$num)) * 52 * 1e5
ir.ct.base <- (sum(mn.base$incid.ct)/sum((1 - mn.base$prev.ct) * mn.base$num)) * 52 * 1e5
ir.syph.base <- (sum(mn.base$incid.syph)/sum((1 - mn.base$prev.syph) * mn.base$num)) * 52 * 1e5
incid.base <- sum(mn.base$incid)
incid.gc.base <- sum(mn.base$incid.gc)
incid.ct.base <- sum(mn.base$incid.ct)
incid.syph.base <- sum(mn.base$incid.syph)

# Sims 3142-3152 for comparisons to all - 6 month interval
# 3142 - STI
# 3143 - recent partners
# 3144 - new partners
# 3145 - partner who has multiple partners
# 3146 - partner with a STI
# 3147 - any CAI in a non-main
# 3148 - any CAI
# 3149 - recent or new partners
# 3150 - sti, recent, or new partners
# 3151 - CAI in non-main or any CAI
# 3152 - partner with multiple partners or with a STI

sims <- c(3014, 3142:3152)
df.hiv.pia <- data.frame(rep(NA, 256))
df.hiv.nnt <- data.frame(rep(NA, 256))
df.gc.pia <- data.frame(rep(NA, 256))
df.gc.nnt <- data.frame(rep(NA, 256))
df.ct.pia <- data.frame(rep(NA, 256))
df.ct.nnt <- data.frame(rep(NA, 256))
df.syph.pia <- data.frame(rep(NA, 256))
df.syph.nnt <- data.frame(rep(NA, 256))

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
    
    vec.hiv.nia <- round(ir.base - unname(ir), 1)
    df.hiv.pia[, i] <- vec.hiv.nia / ir.base
    
    vec.gc.nia <- round(ir.gc.base - unname(ir.gc), 1)
    df.gc.pia[, i] <- vec.gc.nia / ir.gc.base
    
    vec.ct.nia <- round(ir.ct.base - unname(ir.ct), 1)
    df.ct.pia[, i] <- vec.ct.nia / ir.ct.base
    
    vec.syph.nia <- round(ir.syph.base - unname(ir.syph), 1)
    df.syph.pia[, i] <- vec.syph.nia / ir.syph.base
    
    total.hiv.tests <- unname(colMeans(tail(sim$epi$totalhivtests, 1)))
    gc.asympt.tests <- unname(colMeans(tail(sim$epi$totalGCasympttests, 1)))
    ct.asympt.tests <- unname(colMeans(tail(sim$epi$totalCTasympttests, 1)))
    syph.asympt.tests <- unname(colMeans(tail(sim$epi$totalsyphasympttests, 1)))
    total.asympt.tests <- unname(colMeans(tail(sim$epi$totalstiasympttests, 1)))
    
    #HIV could be HIV tests or total STI tests
    df.hiv.nnt[, i] <- total.asympt.tests / (incid.base - unname(colSums(sim$epi$incid)))
    df.gc.nnt[, i] <- gc.asympt.tests / (incid.gc.base - unname(colSums(sim$epi$incid.gc)))
    df.ct.nnt[, i] <- ct.asympt.tests / (incid.ct.base - unname(colSums(sim$epi$incid.ct)))
    df.syph.nnt[, i] <- syph.asympt.tests / (incid.syph.base - unname(colSums(sim$epi$incid.syph)))
    
}
names(df.hiv.pia) <- names(df.gc.pia) <- names(df.gc.nnt) <- names(df.ct.pia) <- 
    names(df.hiv.nnt) <- names(df.ct.nnt) <- names(df.syph.pia) <- names(df.syph.nnt) <- 
    c("All", "C1", "C2", "C3", "C4", "C5", "C6", "C7", "J1", "J2", "J3", "J4")

head(df.hiv.pia)
head(df.gc.pia)
head(df.ct.pia)
head(df.syph.pia)
head(df.hiv.nnt)
head(df.gc.nnt)
head(df.ct.nnt)
head(df.syph.nnt)

pal <- wes_palette("Zissou")[c(1, 5)]
# pdf(file = "analysis/P1Fig2.pdf", height = 6, width = 12, pointsize = 16)

# HIV
tiff(filename = "analysis/Fig1a.tiff", height = 4, width = 10, units = "in", res = 250)
par(mfrow = c(1, 2), mar = c(4,4,2.5,1), oma = c(0, 0, 2, 0), mgp = c(3, 0.75, 0))
boxplot(df.hiv.pia, outline = FALSE, medlwd = 1.1,
        col = c(rep(pal[1], 8), rep(pal[2], 4)),  ylim = c(0, 1),
        main = "PIA by Behavioral Indication", las = 2,
        xlab = "Behavioral Indication", ylab = "Proportion HIV Infections Averted")

boxplot(df.hiv.nnt, outline = FALSE, medlwd = 1.1,
        col = c(rep(pal[1], 8), rep(pal[2], 4)),
        main = "NNT by Behavioral Indication", las = 2,
        xlab = "Behavioral Indication", ylab = "HIV Number Needed to Treat")
title("40% HR Cov (6 months), 0% Ann Cov", outer = TRUE)
dev.off()

# GC
tiff(filename = "analysis/Fig1b.tiff", height = 4, width = 8, units = "in", res = 250)
par(mfrow = c(1, 2), mar = c(4,4,2.5,1), oma = c(0, 0, 2, 0), mgp = c(3, 0.75, 0))
boxplot(df.gc.pia, outline = FALSE, medlwd = 1.1,
        col = c(rep(pal[1], 8), rep(pal[2], 4)), ylim = c(0, 1),
        main = "PIA by Behavioral Indication", las = 2,
        xlab = "Behavioral Indication", ylab = "Proportion NG Infections Averted")
boxplot(df.gc.nnt, outline = FALSE, medlwd = 1.1,
        col = c(rep(pal[1], 8), rep(pal[2], 4)),
        main = "NNT by Behavioral Indication", las = 2,
        xlab = "Behavioral Indication", ylab = "NG Number Needed to Treat")
title("40% HR Cov (6 months), 0% Ann Cov", outer = TRUE)
dev.off()

# CT
tiff(filename = "analysis/Fig1c.tiff", height = 4, width = 8, units = "in", res = 250)
par(mfrow = c(1, 2), mar = c(4,4,2.5,1), oma = c(0, 0, 2, 0), mgp = c(3, 0.75, 0))
boxplot(df.ct.pia, outline = FALSE, medlwd = 1.1,
        col = c(rep(pal[1], 8), rep(pal[2], 4)), ylim = c(0, 1),
        main = "PIA by Behavioral Indication", las = 2,
        xlab = "Behavioral Indication", ylab = "Proportion CT Infections Averted")
boxplot(df.ct.nnt, outline = FALSE, medlwd = 1.1,
        col = c(rep(pal[1], 8), rep(pal[2], 4)),
        main = "NNT by Behavioral Indication",las = 2,
        xlab = "Behavioral Indication", ylab = "CT Number Needed to Treat")
title("40% HR Cov (6 months), 0% Ann Cov", outer = TRUE)
dev.off()

# Syph
tiff(filename = "analysis/Fig1d.tiff", height = 4, width = 8, units = "in", res = 250)
par(mfrow = c(1, 2), mar = c(4,4,2.5,1), oma = c(0, 0, 2, 0), mgp = c(3, 0.75, 0))
boxplot(df.syph.pia, outline = FALSE, medlwd = 1.1,
        col = c(rep(pal[1], 8), rep(pal[2], 4)), ylim = c(0, 1),
        main = "PIA by Behavioral Indication", las = 2,
        xlab = "Behavioral Indication", ylab = "Proportion Syph Infections Averted")
boxplot(df.syph.nnt, outline = FALSE, medlwd = 1.1,
        col = c(rep(pal[1], 8), rep(pal[2], 4)),
        main = "NNT by Behavioral Indication",las = 2,
        xlab = "Behavioral Indication", ylab = "Syph Number Needed to Treat")
title("40% HR Cov (6 months), 0% Ann Cov", outer = TRUE)
dev.off()

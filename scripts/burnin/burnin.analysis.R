library("EpiModelHPC")
library("EpiModelHIV")

## Incidence

ir100.gc <- as.numeric(sim$epi$ir100.gc[2600, ])
round(quantile(ir100.gc, probs = c(0.025, 0.5, 0.975)), 3)

ir100.ct <- as.numeric(sim$epi$ir100.ct[2600, ])
round(quantile(ir100.ct, probs = c(0.025, 0.5, 0.975)), 3)

ir100 <- as.numeric(sim$epi$ir100[2600, ])
round(quantile(ir100, probs = c(0.025, 0.5, 0.975)), 3)

ir100.syph <- as.numeric(sim$epi$ir100.syph[2600, ])
round(quantile(ir100.syph, probs = c(0.025, 0.5, 0.975)), 3)

i.prev <- as.numeric(sim$epi$i.prev[2600, ])
round(quantile(i.prev, probs = c(0.025, 0.5, 0.975)), 3)

# Syphilis stage prevalence and HIV-syphilis overlap

prev.primsecosyph.hivpos <- as.numeric(sim$epi$prev.primsecosyph.hivpos[2600, ])
round(quantile(prev.primsecosyph.hivpos, probs = c(0.025, 0.5, 0.975)), 3)

prev.primsecosyph.hivneg <- as.numeric(sim$epi$prev.primsecosyph.hivneg[2600, ])
round(quantile(prev.primsecosyph.hivneg, probs = c(0.025, 0.5, 0.975)), 3)

prev.syph.hivpos <- as.numeric(sim$epi$prev.syph.hivpos[2600, ])
round(quantile(prev.syph.hivpos, probs = c(0.025, 0.5, 0.975)), 3)

prev.syph.hivneg <- as.numeric(sim$epi$prev.syph.hivneg[2600, ])
round(quantile(prev.syph.hivneg, probs = c(0.025, 0.5, 0.975)), 3)

prev.syph <- as.numeric(sim$epi$prev.syph[2600, ])
round(quantile(prev.syph, probs = c(0.025, 0.5, 0.975)), 3)

prev.primsecosyph <- as.numeric(sim$epi$prev.primsecosyph[2600, ])
round(quantile(prev.primsecosyph, probs = c(0.025, 0.5, 0.975)), 3)

prev.hiv.syphpos <- as.numeric(sim$epi$prev.hiv.syphpos[2600, ])
round(quantile(prev.hiv.syphpos, probs = c(0.025, 0.5, 0.975)), 3)

prev.hiv.primsecosyphpos <- as.numeric(sim$epi$prev.hiv.primsecosyphpos[2600, ])
round(quantile(prev.hiv.primsecosyphpos, probs = c(0.025, 0.5, 0.975)), 3)

prev.earlysyph <- as.numeric(sim$epi$prev.earlysyph[2600, ])
round(quantile(prev.earlysyph, probs = c(0.025, 0.5, 0.975)), 3)

prev.latesyph <- as.numeric(sim$epi$prev.latesyph[2600, ])
round(quantile(prev.latesyph, probs = c(0.025, 0.5, 0.975)), 3)

prev.stage.incubprim <- as.numeric(sim$epi$prev.stage.incubprim[2600, ])
round(quantile(prev.stage.incubprim, probs = c(0.025, 0.5, 0.975)), 3)

prev.stage.seco <- as.numeric(sim$epi$prev.stage.seco[2600, ])
round(quantile(prev.stage.seco, probs = c(0.025, 0.5, 0.975)), 3)

prev.stage.earlat <- as.numeric(sim$epi$prev.stage.earlat[2600, ])
round(quantile(prev.stage.earlat, probs = c(0.025, 0.5, 0.975)), 3)

prev.stage.latelat <- as.numeric(sim$epi$prev.stage.latelat[2600, ])
round(quantile(prev.stage.latelat, probs = c(0.025, 0.5, 0.975)), 3)

prev.stage.latelatelat <- as.numeric(sim$epi$prev.stage.latelatelat[2600, ])
round(quantile(prev.stage.latelatelat, probs = c(0.025, 0.5, 0.975)), 3)

prev.stage.alllatelat <- as.numeric(sim$epi$prev.stage.alllatelat[2600, ])
round(quantile(prev.stage.alllatelat, probs = c(0.025, 0.5, 0.975)), 3)

prev.stage.tert <- as.numeric(sim$epi$prev.stage.tert[2600, ])
round(quantile(prev.stage.tert, probs = c(0.025, 0.5, 0.975)), 3)


prev.gc <- as.numeric(sim$epi$prev.gc[2600, ])
round(quantile(prev.gc, probs = c(0.025, 0.5, 0.975)), 3)

prev.rgc <- as.numeric(sim$epi$prev.rgc[2600, ])
round(quantile(prev.rgc, probs = c(0.025, 0.5, 0.975)), 3)

prev.ugc <- as.numeric(sim$epi$prev.ugc[2600, ])
round(quantile(prev.ugc, probs = c(0.025, 0.5, 0.975)), 3)

prev.ct <- as.numeric(sim$epi$prev.ct[2600, ])
round(quantile(prev.ct, probs = c(0.025, 0.5, 0.975)), 3)

prev.rct <- as.numeric(sim$epi$prev.rct[2600, ])
round(quantile(prev.rct, probs = c(0.025, 0.5, 0.975)), 3)

prev.uct <- as.numeric(sim$epi$prev.uct[2600, ])
round(quantile(prev.uct, probs = c(0.025, 0.5, 0.975)), 3)

prev.rgcct <- as.numeric(sim$epi$prev.rgcct[2600, ])
round(quantile(prev.rgcct, probs = c(0.025, 0.5, 0.975)), 3)

prev.ugcct <- as.numeric(sim$epi$prev.ugcct[2600, ])
round(quantile(prev.ugcct, probs = c(0.025, 0.5, 0.975)), 3)

# Summary of 500 sims
# Syphilis prevalence plots
par(mfrow = c(1, 1), oma = c(0,0,2,0))
plot(sim, y = "prev.primsecosyph.hivpos", ylab = "Prevalence", ylim = c(0, 0.15), mean.col = "blue", qnts.col = "blue")
plot(sim, y = "prev.primsecosyph.hivneg", ylab = "Prevalence", add = TRUE, mean.col = "red", qnts.col = "red")
abline(h = 0.103, col = "blue", lty = 2)
abline(h = 0.026, col = "red", lty = 2)
title("P and S Syphilis by HIV Status")
legend("topleft", c("HIV+", "HIV-"), col = c("blue", "red"), lty = c(1, 1))

plot(sim, y = "prev.stage.incubprim", ylab = "Prevalence", ylim = c(0.00, 0.40), mean.col = "blue", qnts.col = "blue")
plot(sim, y = "prev.stage.seco", ylab = "Prevalence", add = TRUE, mean.col = "red", qnts.col = "red")
plot(sim, y = "prev.stage.earlat", ylab = "Prevalence", add = TRUE, mean.col = "purple", qnts.col = "purple")
abline(h = 0.20, col = "blue", lty = 2)
abline(h = 0.077, col = "red", lty = 2)
abline(h = 0.2770, col = "purple", lty = 2)
title("P, S, and Early latent Syphilis Prevalence")
legend("topright", c("Primary", "Seco", "Early latent"), col = c("blue", "red", "purple"), lty = c(1, 1, 1))

plot(sim, y = "prev.stage.alllatelat", ylab = "Prevalence", ylim = c(0.00, 0.60), mean.col = "blue", qnts.col = "blue")
plot(sim, y = "prev.stage.tert", ylab = "Prevalence", add = TRUE, mean.col = "red", qnts.col = "red")
abline(h = 0.44, col = "blue", lty = 2)
abline(h = 0.006, col = "red", lty = 2)
title("Late Latent and Tertiary Syphilis Prevalence")
legend("right", c("Late latent", "Tertiary"), col = c("blue", "red"), lty = c(1, 1))

# All syphilis stages
par(mfrow = c(1, 1), oma = c(0, 0, 2, 0))
plot(sim, y = "prev.stage.incubprim", ylab = "Prevalence", ylim = c(0, 1), mean.col = "blue", qnts.col = "blue")
plot(sim, y = "prev.stage.seco", ylab = "Prevalence", add = TRUE, mean.col = "red", qnts.col = "red")
plot(sim, y = "prev.stage.earlat", ylab = "Prevalence", add = TRUE, mean.col = "orange", qnts.col = "orange")
plot(sim, y = "prev.stage.alllatelat", ylab = "Prevalence", add = TRUE, mean.col = "purple", qnts.col = "purple")
plot(sim, y = "prev.stage.tert", ylab = "Prevalence", add = TRUE, mean.col = "black", qnts.col = "black")
title("Stage-specific Syphilis Prevalence")
abline(h = 0.20, col = "blue", lty = 1)
abline(h = 0.077, col = "red", lty = 1)
abline(h = 0.2770, col = "orange", lty = 1)
abline(h = 0.44, col = "purple", lty = 1)
abline(h = 0.006, col = "black", lty = 1)
legend("topright", c("Primary/Incub", "Secondary", "Early Latent", "Late latent", "Tertiary"),
       col = c("blue", "red", "orange", "purple", "black"), lty = c(1, 1, 1, 1, 1))

# Early vs late syphilis
par(mfrow = c(1, 1), oma = c(0, 0, 2, 0))
plot(sim, y = "prev.earlysyph", ylab = "Prevalence", ylim = c(0, 1), mean.col = "blue", qnts.col = "blue")
plot(sim, y = "prev.latesyph", ylab = "Prevalence", add = TRUE, mean.col = "red", qnts.col = "red")
title("Early and Late-stage Syphilis")
abline(h = 0.554, col = "blue", lty = 1)
abline(h = 0.446, col = "red", lty = 1)
legend("topright", c("Early", "Late"), col = c("blue", "red"), lty = c(1, 1))

# STI Prevalence
par(mfrow = c(2,2), oma = c(0,0,2,0))
plot(sim, y = "prev.primsecosyph", ylab = "Prevalence")
title("P and S Syphilis Prevalence")
abline(h = 0.046, col = "red", lty = 2)
plot(sim, y = "prev.ct", ylab = "Prevalence")
title("CT Prevalence")
plot(sim, y = "prev.gc", ylab = "Prevalence")
title("GC Prevalence")
plot(sim, y = "i.prev", ylim = c(0, 0.3), ylab = "Prevalence")
abline(h = 0.26, col = "red",  lty = 2)
title("HIV Prevalence")

# HIV and STI Incidence
par(mfrow = c(2,2), oma = c(0,0,2,0))
plot(sim, y = "ir100")
abline(h = 3.8, col = "red", lty = 2)
title("HIV Incidence")
plot(sim, y = "ir100.gc")
abline(h = 4.2, col = "red", lty = 2)
title("GC Incidence")
plot(sim, y = "ir100.ct")
abline(h = 6.6, col = "red", lty = 2)
title("CT Incidence")
plot(sim, y = "ir100.syph")
abline(h = 2.0, col = "red", lty = 2)
title("Syph Incidence")

par(mfrow = c(1, 2), oma = c(0, 0, 2, 0))
plot(sim, y = "hiv_sum", mean.col = "blue", qnts.col = "blue", qnts.alpha = 0.2)
plot(sim, y = "sti_hiv_sum", mean.col = "green", qnts.col = "green", qnts.alpha = 0.2, add = TRUE)
plot(sim, y = "sti_u_hiv_sum", mean.col = "red", qnts.col = "red", qnts.alpha = 0.2,add = TRUE)
plot(sim, y = "sti_r_hiv_sum", mean.col = "orange", qnts.col = "orange", qnts.alpha = 0.2,add = TRUE)
plot(sim, y = "sti_syph_hiv_sum", mean.col = "purple", qnts.col = "purple", qnts.alpha = 0.2,add = TRUE)
legend("topleft", lty = c(1, 1, 1, 1, 1), col = c("blue", "green", "red", "orange", "purple"), 
       c("HIV Sum", "STI & HIV sum", "USTI & HIV sum", "RSTI & HIV sum", "Syph & HIV sum"))
title("HIV Infections per time step")

par(mfrow = c(1, 1))
plot(sim, y = "sti_paf", mean.col = "blue", qnts.col = "blue", qnts = c(0.5), qnts.alpha = 0.05)
abline(h = c(seq(0.1, 0.9, 0.1)), lty = 2, col = "gray")
plot(sim, y = "sti_u_paf", mean.col = "green", qnts.col = "green", qnts = c(0.5), qnts.alpha = 0.05, add = TRUE)
plot(sim, y = "sti_r_paf", mean.col = "red", qnts.col = "red", qnts = c(0.5), qnts.alpha = 0.05, add = TRUE)
plot(sim, y = "sti_syph_paf", mean.col = "orange", qnts.col = "orange", qnts = c(0.5), qnts.alpha = 0.05, add = TRUE)
legend("topleft", lty = c(1, 1, 1, 1), col = c("blue", "green", "red", "orange"), 
       c("Total STI PAF", "Urethral NG/CT PAF", "Rectal NG/CT PAF", "Syph PAF"))
title("STI Population Attributable Fraction (PAF) and IQR for HIV Infection")

plot(sim, y = "sti_u_paf", mean.col = "blue")
plot(sim, y = "sti_u_sympt_paf", mean.col = "green", add = TRUE)
plot(sim, y = "sti_u_asympt_paf", mean.col = "red", add = TRUE)
legend("topleft", lty = c(1, 1, 1), col = c("blue", "green", "red"), 
       c("U STI PAF", "U Sympt STI PAF", "U Asympt STI PAF"))

plot(sim, y = "sti_r_paf", mean.col = "blue")
plot(sim, y = "sti_r_sympt_paf", mean.col = "green", add = TRUE)
plot(sim, y = "sti_r_asympt_paf", mean.col = "red", add = TRUE)
legend("topleft", lty = c(1, 1, 1), col = c("blue", "green", "red"), 
       c("R STI PAF", "R Sympt STI PAF", "R Asympt STI PAF"))

plot(sim, y = "sti_syph_paf", mean.col = "blue")
plot(sim, y = "sti_syph_sympt_paf", mean.col = "green", add = TRUE)
plot(sim, y = "sti_syph_asympt_paf", mean.col = "red", add = TRUE)
legend("topleft", lty = c(1, 1, 1), col = c("blue", "green", "red"), 
       c("Syph PAF", "Syph Sympt PAF", "Syph Asympt PAF"))

# Evaluate whether testing is leading to treatment
a <- cbind(sim$epi$txCT, sim$epi$CTsympttests, sim$epi$CTasympttests.pos,
           sim$epi$txGC, sim$epi$GCsympttests, sim$epi$GCasympttests.pos,
           sim$epi$txsyph, sim$epi$syphsympttests, sim$epi$syphasympttests.pos)
colnames(a) <- c("txCT", "CTsympttests", "CTposasympttests",
                 "txGC", "GCsympttests", "GCposasympttests",
                 "txsyph", "syphsympttests", "syphposasympttests") 
tail(a)

par(mfrow = c(1,2), oma = c(0,0,2,0))
plot(sim, y = "GCasympttests", mean.col = "blue")
plot(sim, y = "CTasympttests", mean.col = "red", add = TRUE)
plot(sim, y = "syphasympttests", mean.col = "green", add = TRUE)
plot(sim, y = "GCasympttests.pos", mean.col = "blue", add = TRUE, mean.lty = 2)
plot(sim, y = "CTasympttests.pos", mean.col = "red", add = TRUE, mean.lty = 2)
plot(sim, y = "syphasympttests.pos", mean.col = "green", add = TRUE, mean.lty = 2)
abline(v = sim$param$stitest.start)
legend("topleft", lty = c(1, 1, 1, 2, 2, 2), col = c("blue", "red", "green", "blue", "red", "green"), 
       c("GC Tests", "CT Tests", "Syph Tests", "Pos GC Tests", "Pos CT Tests", "Pos Syph Tests"))
title("Asymptomatic Tests")

plot(sim, y = "txGC", mean.col = "blue", mean.lty = 2, ylim = c(0, 100))
plot(sim, y = "txCT", mean.col = "red", mean.lty = 2, add = TRUE)
plot(sim, y = "txsyph", mean.col = "green", mean.lty = 2, add = TRUE)
plot(sim, y = "GCsympttests", mean.col = "blue", add = TRUE)
plot(sim, y = "CTsympttests", mean.col = "red", add = TRUE)
plot(sim, y = "syphsympttests", mean.col = "green", add = TRUE)
abline(v = sim$param$stitest.start)
legend("topleft", lty = c(1, 1, 1, 2, 2, 2), col = c("blue", "red", "green", "blue", "red", "green"), 
       c("GC Tests", "CT Tests", "Syph Tests", "Tx GC", "Tx CT", "TX Syph"))
title("Symptomatic Tests and All Treated")

plot(sim, y = "txGC", mean.col = "blue", mean.lty = 2)
plot(sim, y = "GCsympttests", mean.col = "red", add = TRUE)
plot(sim, y = "GCasympttests.pos", mean.col = "green", add = TRUE)
abline(v = sim$param$stitest.start)
legend("topleft", lty = c(2, 1, 1), col = c("blue", "red", "green"), 
       c("Tx GC", "GC Sympt Tests", "GC Pos Asympt Tests"))
title("GC")

plot(sim, y = "txCT", mean.col = "blue", mean.lty = 2, ylim = c(0, 60))
plot(sim, y = "CTsympttests", mean.col = "red", add = TRUE)
plot(sim, y = "CTasympttests.pos", mean.col = "green", add = TRUE)
abline(v = sim$param$stitest.start)
legend("topleft", lty = c(2, 1, 1), col = c("blue", "red", "green"), 
       c("Tx CT", "CT Sympt Tests", "CT Pos Asympt Tests"))
title("CT")

plot(sim, y = "prev.rct", mean.col = "purple", ylim = c(0, 0.15))
plot(sim, y = "prev.uct", mean.col = "orange", add = TRUE)
plot(sim, y = "prev.syph", mean.col = "green", add = TRUE)
plot(sim, y = "prev.rgc", mean.col = "blue", add = TRUE)
plot(sim, y = "prev.ugc", mean.col = "red", add = TRUE)
abline(v = sim$param$stitest.start)
legend("topleft", lty = c(1, 1, 1, 1, 1), col = c("purple", "orange", "green","blue", "red"), 
       c("rCT", "uCT", "Syph", "rGC", "uGC"))
title("Prevalence")

plot(sim, y = "ir100.rct", mean.col = "purple", ylim = c(0, 40))
plot(sim, y = "ir100.uct", mean.col = "orange", add = TRUE)
plot(sim, y = "ir100.syph", mean.col = "green", add = TRUE)
plot(sim, y = "ir100.rgc", mean.col = "blue", add = TRUE)
plot(sim, y = "ir100.ugc", mean.col = "red", add = TRUE)
abline(v = sim$param$stitest.start)
legend("topleft", lty = c(1, 1, 1, 1, 1), col = c("purple", "orange", "green","blue", "red"), 
       c("rCT", "uCT", "Syph", "rGC", "uGC"))
title("Incidence")

# STI Testing Indications
par(mfrow = c(1,1), oma = c(0,0,2,0))
plot(sim, y = "stiactiveind", mean.col = "purple", ylim = c(0, 1))#, xlim = c(0, sim$control$nsteps - sim$param$stitest.start), ylim = c(0, 1))
plot(sim, y = "newpartner", mean.col = "orange", add = TRUE)
plot(sim, y = "recentpartners", mean.col = "green", add = TRUE)
plot(sim, y = "concurrpart", mean.col = "blue", add = TRUE)
plot(sim, y = "partnersti", mean.col = "brown", add = TRUE)
plot(sim, y = "uai.nmain", mean.col = "black", add = TRUE)
plot(sim, y = "uai.any", mean.col = "gray", add = TRUE)
plot(sim, y = "recentSTI", mean.col = "red", add = TRUE)
abline(h = c(seq(0.1, 0.9, 0.1)), lty = 2, col = "gray")
legend(400, 0.9, lty = c(rep(1, 8)), 
       col = c("purple", "orange", "green", "blue", "brown", "black", "gray", "red"),
       c("Sexually Active", "New Partner", ">1 Recent Partners", "Partner is Concurrent", 
         "Partner had STI", "CAI in Non-main", "Any CAI", "Recent STI"))
title("STI Testing Indications")


###############################
# Summary of mean simulation
# Run after burnin.process
###############################

# Syphilis prevalence plots
par(mfrow = c(1, 1), oma = c(0,0,2,0))
plot(sim, y = "prev.primsecosyph.hivpos", ylab = "Prevalence", ylim = c(0, 0.15), mean.col = "blue", qnts.col = "blue")
plot(sim, y = "prev.primsecosyph.hivneg", ylab = "Prevalence", add = TRUE, mean.col = "red", qnts.col = "red")
abline(h = 0.103, col = "blue", lty = 2)
abline(h = 0.026, col = "red", lty = 2)
title("P and S Syphilis by HIV Status")
legend("topleft", c("HIV+", "HIV-"), col = c("blue", "red"), lty = c(1, 1))

plot(sim, y = "prev.stage.incubprim", ylab = "Prevalence", ylim = c(0.00, 0.40), mean.col = "blue", qnts.col = "blue")
plot(sim, y = "prev.stage.seco", ylab = "Prevalence", add = TRUE, mean.col = "red", qnts.col = "red")
plot(sim, y = "prev.stage.earlat", ylab = "Prevalence", add = TRUE, mean.col = "purple", qnts.col = "purple")
abline(h = 0.20, col = "blue", lty = 2)
abline(h = 0.077, col = "red", lty = 2)
abline(h = 0.2770, col = "purple", lty = 2)
title("P, S, and Early latent Syphilis Prevalence")
legend("topright", c("Primary", "Seco", "Early latent"), col = c("blue", "red", "purple"), lty = c(1, 1, 1))

plot(sim, y = "prev.stage.alllatelat", ylab = "Prevalence", ylim = c(0.00, 0.60), mean.col = "blue", qnts.col = "blue")
plot(sim, y = "prev.stage.tert", ylab = "Prevalence", add = TRUE, mean.col = "red", qnts.col = "red")
abline(h = 0.44, col = "blue", lty = 2)
abline(h = 0.006, col = "red", lty = 2)
title("Late Latent and Tertiary Syphilis Prevalence")
legend("right", c("Late latent", "Tertiary"), col = c("blue", "red"), lty = c(1, 1))

# All syphilis stages
par(mfrow = c(1, 1), oma = c(0, 0, 2, 0))
plot(sim, y = "prev.stage.incubprim", ylab = "Prevalence", ylim = c(0, 1), mean.col = "blue", qnts.col = "blue")
plot(sim, y = "prev.stage.seco", ylab = "Prevalence", add = TRUE, mean.col = "red", qnts.col = "red")
plot(sim, y = "prev.stage.earlat", ylab = "Prevalence", add = TRUE, mean.col = "orange", qnts.col = "orange")
plot(sim, y = "prev.stage.alllatelat", ylab = "Prevalence", add = TRUE, mean.col = "purple", qnts.col = "purple")
plot(sim, y = "prev.stage.tert", ylab = "Prevalence", add = TRUE, mean.col = "black", qnt.col = "black")
title("Stage-specific Syphilis Prevalence")
abline(h = 0.20, col = "blue", lty = 1)
abline(h = 0.077, col = "red", lty = 1)
abline(h = 0.2770, col = "orange", lty = 1)
abline(h = 0.44, col = "purple", lty = 1)
abline(h = 0.006, col = "black", lty = 1)
legend("topright", c("Primary/Incub", "Secondary", "Early Latent", "Late latent", "Tertiary"),
       col = c("blue", "red", "orange", "purple", "black"), lty = c(1, 1, 1, 1, 1))

# Early vs late syphilis
par(mfrow = c(1, 1), oma = c(0, 0, 2, 0))
plot(sim, y = "prev.earlysyph", ylab = "Prevalence", ylim = c(0, 1), mean.col = "blue", qnts.col = "blue")
plot(sim, y = "prev.latesyph", ylab = "Prevalence", add = TRUE, mean.col = "red", qnts.col = "red")
title("Early and Late-stage Syphilis")
abline(h = 0.554, col = "blue", lty = 1)
abline(h = 0.446, col = "red", lty = 1)
legend("topright", c("Early", "Late"), col = c("blue", "red"), lty = c(1, 1))

# STI Prevalence
par(mfrow = c(2,2), oma = c(0,0,2,0))
plot(sim, y = "prev.primsecosyph", ylab = "Prevalence")
title("P and S Syphilis Prevalence")
abline(h = 0.046, col = "red", lty = 2)
plot(sim, y = "prev.ct", ylab = "Prevalence")
title("CT Prevalence")
plot(sim, y = "prev.gc", ylab = "Prevalence")
title("GC Prevalence")
plot(sim, y = "i.prev", ylim = c(0, 0.3), ylab = "Prevalence")
abline(h = 0.26, col = "red",  lty = 2)
title("HIV Prevalence")

# HIV and STI Incidence
par(mfrow = c(2,2), oma = c(0,0,2,0))
plot(sim, y = "ir100")
abline(h = 3.8, col = "red", lty = 2)
title("HIV Incidence")
plot(sim, y = "ir100.gc")
abline(h = 4.2, col = "red", lty = 2)
title("GC Incidence")
plot(sim, y = "ir100.ct")
abline(h = 6.6, col = "red", lty = 2)
title("CT Incidence")
plot(sim, y = "ir100.syph")
abline(h = 0.9, col = "red", lty = 2)
title("Syph Incidence")

df <- as.data.frame(sim)
head(df$prev.rgcct, 100)

plot(sim, y = "prev.rgcct", mean.smooth = FALSE)

df <- tail(as.data.frame(sim), 500)
sum(df$trans.main) / sum(df$incid)
sum(df$trans.casl) / sum(df$incid)
sum(df$trans.inst) / sum(df$incid)

# sum(df$trans.main.gc) / sum(df$incid.gc)
# sum(df$trans.casl.gc) / sum(df$incid.gc)
# sum(df$trans.inst.gc) / sum(df$incid.gc)
# 
# sum(df$trans.main.ct) / sum(df$incid.ct)
# sum(df$trans.casl.ct) / sum(df$incid.ct)
# sum(df$trans.inst.ct) / sum(df$incid.ct)
# 
# sum(df$trans.main.syph) / sum(df$incid.syph)
# sum(df$trans.casl.syph) / sum(df$incid.syph)
# sum(df$trans.inst.syph) / sum(df$incid.syph)

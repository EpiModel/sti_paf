# Syphilis ratio
boxplot(sim$epi$prev.syph.hivpos / sim$epi$prev.syph.hivneg)

# Synergy plots
par(mfrow=c(2,3), oma=c(0,0,2,0))
plot(sim, y = "prev.syph.hivpos", ylab = "Prevalence", ylim=c(0, 0.12))
plot(sim, y = "prev.syph.hivneg", ylab = "Prevalence", add = TRUE)
abline(h = 0.103, col = "red", lty = 2)
abline(h = 0.026, col = "red", lty = 2)
title("Syphilis by HIV Status")
plot(sim, y = "prev.gc.hivpos", ylab = "Prevalence")
plot(sim, y = "prev.gc.hivneg", ylab = "Prevalence", add = TRUE)
title("GC by HIV Status")
plot(sim, y = "prev.ct.hivpos", ylab = "Prevalence")
plot(sim, y = "prev.ct.hivneg", ylab = "Prevalence", add = TRUE)
title("CT by HIV Status")
plot(sim, y = "prev.hiv.syphpos", ylab = "Prevalence")
plot(sim, y = "prev.hiv.syphneg", ylab = "Prevalence", add = TRUE)
title("HIV by Syphilis +/-")
plot(sim, y = "prev.hiv.gcpos", ylab = "Prevalence")
plot(sim, y = "prev.hiv.gcneg", ylab = "Prevalence", add = TRUE)
title("HIV by GC +/-")
plot(sim, y = "prev.hiv.ctpos", ylab = "Prevalence")
plot(sim, y = "prev.hiv.ctneg", ylab = "Prevalence", add = TRUE)
title("HIV by CT +/-")
title("Syph Tprob = XXX, Syph.HIV.RR = , HIV.Syph.RR =", outer = TRUE)



#Prevalence
par(mfrow=c(2,2), oma=c(0,0,2,0))
plot(sim, y = "prev.syph", ylab = "Prevalence")
title("Syphilis Prevalence")
plot(sim, y = "prev.ct", ylab = "Prevalence")
title ("CT Prevalence")
plot(sim, y = "prev.gc", ylab = "Prevalence")
title("GC Prevalence")
plot(sim, y = "i.prev", ylab = "Prevalence")
abline(h=0.26, col = "red", lty = 2)
title("HIV Prevalence")
title("Syph Tprob = XXX, Syph.HIV.RR = , HIV.Syph.RR =", outer = TRUE)


# Incidence
par(mfrow=c(2,2), oma=c(0,0,2,0))
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
title("Syph Tprob = XXX, Relrisk for Syph<->HIV = XXX", outer = TRUE)

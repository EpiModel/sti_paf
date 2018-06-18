rm(list = ls())
library("EpiModelHIV")
library("EpiModelHPC")
library("dplyr")

load("data/followup/sim.n5001.rda")
#25% reduction in STI screening among HIV-
load("data/followup/sim.n5064.rda")
#50% reduction in STI screening among HIV-
load("data/followup/sim.n5065.rda")

load("est/stimod.burnin.rda")
load("est/stidiagmodchange.burnin.rda")


# PAF-related calcs ------------------------------------
# NG only
quantile(colMeans(tail(sim$epi$prev.rgc.hivneg.only, 52)) + colMeans(tail(sim$epi$prev.ugc.hivneg.only, 52)))
quantile(colMeans(tail(sim$epi$prev.rgc.hivpos.only, 52)) + colMeans(tail(sim$epi$prev.ugc.hivpos.only, 52)))

#any NG
quantile(colMeans(tail(sim$epi$prev.gc.hivneg, 52)))
quantile(colMeans(tail(sim$epi$prev.gc.hivpos, 52)))

# CT only
quantile(colMeans(tail(sim$epi$prev.rct.hivneg.only, 52)) + colMeans(tail(sim$epi$prev.uct.hivneg.only, 52)))
quantile(colMeans(tail(sim$epi$prev.rct.hivpos.only, 52)) + colMeans(tail(sim$epi$prev.uct.hivpos.only, 52)))

#any CT
quantile(colMeans(tail(sim$epi$prev.ct.hivneg, 52)))
quantile(colMeans(tail(sim$epi$prev.ct.hivpos, 52)))

# Syph only
quantile(colMeans(tail(sim$epi$prev.primsecosyph.hivneg.only, 52)))
quantile(colMeans(tail(sim$epi$prev.primsecosyph.hivpos.only, 52)))

# Any syphilis
quantile(colMeans(tail(sim$epi$prev.primsecosyph.hivneg, 52)))
quantile(colMeans(tail(sim$epi$prev.primsecosyph.hivpos, 52)))

# Multiple STI
quantile(colMeans(tail(sim$epi$prev.hivnegmultsti, 52)))
quantile(colMeans(tail(sim$epi$prev.hivposmultsti, 52)))

# Number of edges
quantile(colMeans(tail(sim$epi$prop.edges.negneg, 52)))
quantile(colMeans(tail(sim$epi$prop.edges.negpos, 52)))
quantile(colMeans(tail(sim$epi$prop.edges.pospos, 52)))

# Number of acts by partnership serostatuses
quantile(colMeans(tail(sim$epi$num.acts.negneg, 52)))
quantile(colMeans(tail(sim$epi$num.acts.negpos, 52)))
quantile(colMeans(tail(sim$epi$num.acts.pospos, 52)))

# Proportion of UAI acts by partnership serostatuses
quantile(colMeans(tail(sim$epi$prop.uai.negneg, 52)))
quantile(colMeans(tail(sim$epi$prop.uai.negpos, 52)))
quantile(colMeans(tail(sim$epi$prop.uai.pospos, 52)))

# Proportion of acts by partnership serostatuses
quantile(colMeans(tail(sim$epi$prop.acts.negneg, 52)))
quantile(colMeans(tail(sim$epi$prop.acts.negpos, 52)))
quantile(colMeans(tail(sim$epi$prop.acts.pospos, 52)))

# Syph IR by HIV status
quantile(colMeans(tail(sim$epi$ir100.syph.hivneg, 52)))
quantile(colMeans(tail(sim$epi$ir100.syph.hivpos, 52)))

# Prevalence of HIV among IPS syphilis positive
quantile(colMeans(tail(sim$epi$prev.hiv.primsecosyphpos, 52)))

# Prevalence of diagnosed HIV among IPS syphilis positive
quantile(colMeans(tail(sim$epi$prev.dxhiv.dxipssyph, 52)))

# Rectal and urethral STI
quantile(colMeans(tail(sim$epi$prev.rgcct, 52)))
quantile(colMeans(tail(sim$epi$prev.ugcct, 52)))

# Total GC prevalence
quantile(colMeans(tail(sim$epi$prev.rgc.hivneg.only, 52)) +
           colMeans(tail(sim$epi$prev.ugc.hivneg.only, 52)) +
           colMeans(tail(sim$epi$prev.rgc.hivpos.only, 52)) +
           colMeans(tail(sim$epi$prev.ugc.hivpos.only, 52)))

# Total CT prevalence
quantile(colMeans(tail(sim$epi$prev.rct.hivneg.only, 52)) +
           colMeans(tail(sim$epi$prev.uct.hivneg.only, 52)) +
           colMeans(tail(sim$epi$prev.rct.hivpos.only, 52)) +
           colMeans(tail(sim$epi$prev.uct.hivpos.only, 52)))

# Total syph prevalence
quantile(colMeans(tail(sim$epi$prev.syph.hivneg.only, 52)) +
           colMeans(tail(sim$epi$prev.syph.hivpos.only, 52)))

# Total multiple STI prevalence
quantile(colMeans(tail(sim$epi$prev.hivnegmultsti, 52)) +
           colMeans(tail(sim$epi$prev.hivposmultsti, 52)))
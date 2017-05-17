
## P2 functions

gc_pia_vec <- function(sim.base, sim.comp) {

  sim.base <- truncate_sim(sim.base, at = 2600)
  sim.comp <- truncate_sim(sim.comp, at = 2600)

  ir.base.gc <- unname(colMeans(sim.base$epi$ir100.gc)) * 1000

  ir.comp.gc <- unname(colMeans(sim.comp$epi$ir100.gc, na.rm = TRUE)) * 1000
  vec.nia.gc <- round(ir.base.gc - ir.comp.gc, 1)

  vec.pia.gc <- vec.nia.gc/ir.base.gc

  return(vec.pia.gc)
}


epi_stats <- function(sim.base,
                      sim.comp = NULL,
                      at = 520,
                      qnt.low = 0.025,
                      qnt.high = 0.975) {

  # Base scenario -------------------------------------------------------

  # sim.base <- truncate_sim(sim.base, at = 2600)

  # prevalence
  i.prev <- as.numeric(sim.base$epi$i.prev[at, ])
  prev.base <- round(data.frame(median = median(i.prev),
                                ql = quantile(i.prev, qnt.low, names = FALSE),
                                qu = quantile(i.prev, qnt.high, names = FALSE)), 3)
  
  primsecosyph.prev <- as.numeric(sim.base$epi$prev.primsecosyph[at, ])
  primsecosyph.prev.base <- round(data.frame(median = median(primsecosyph.prev ),
                                            ql = quantile(primsecosyph.prev, qnt.low, names = FALSE),
                                            qu = quantile(primsecosyph.prev, qnt.high, names = FALSE)), 3)

  # incidence
  haz <- as.numeric(colMeans(tail(sim.base$epi$ir100, 52)))
  haz.base <- round(data.frame(median = median(haz),
                               ql = quantile(haz, qnt.low, names = FALSE),
                               qu = quantile(haz, qnt.high, names = FALSE)), 2)

  ir.base <- unname(colMeans(sim.base$epi$ir100)) * 1000
  incid.base <- unname(colSums(sim.base$epi$incid))

  haz.gc <- as.numeric(colMeans(tail(sim.base$epi$ir100.gc, 52)))
  haz.base.gc <- round(data.frame(median = median(haz.gc),
                                  ql = quantile(haz.gc, qnt.low, names = FALSE),
                                  qu = quantile(haz.gc, qnt.high, names = FALSE)), 2)
  ir.base.gc <- unname(colMeans(sim.base$epi$ir100.gc)) * 1000
  incid.base.gc <- unname(colSums(sim.base$epi$incid.gc))

  haz.ct <- as.numeric(colMeans(tail(sim.base$epi$ir100.ct, 52)))
  haz.base.ct <- round(data.frame(median = median(haz.ct),
                                  ql = quantile(haz.ct, qnt.low, names = FALSE),
                                  qu = quantile(haz.ct, qnt.high, names = FALSE)), 2)
  ir.base.ct <- unname(colMeans(sim.base$epi$ir100.ct)) * 1000
  incid.base.ct <- unname(colSums(sim.base$epi$incid.ct))
  
  haz.syph <- as.numeric(colMeans(tail(sim.base$epi$ir100.syph, 52)))
  haz.base.syph <- round(data.frame(median = median(haz.syph),
                                  ql = quantile(haz.syph, qnt.low, names = FALSE),
                                  qu = quantile(haz.syph, qnt.high, names = FALSE)), 2)
  ir.base.syph <- unname(colMeans(sim.base$epi$ir100.syph)) * 1000
  incid.base.syph <- unname(colSums(sim.base$epi$incid.syph))


  # Comparison scenario -------------------------------------------------


  if (!is.null(sim.comp)) {

    # sim.comp <- truncate_sim(sim.comp, at = 2600)

    # prevalence
    i.prev <- as.numeric(sim.comp$epi$i.prev[at, ])
    out.prev <- round(data.frame(median = median(i.prev),
                                  ql = quantile(i.prev, qnt.low, names = FALSE),
                                  qu = quantile(i.prev, qnt.high, names = FALSE)), 3)
    
    primsecosyph.prev <- as.numeric(sim.comp$epi$prev.primsecosyph[at, ])
    out.primsecosyph.prev <- round(data.frame(median = median(primsecosyph.prev ),
                                 ql = quantile(primsecosyph.prev, qnt.low, names = FALSE),
                                 qu = quantile(primsecosyph.prev, qnt.high, names = FALSE)), 3)
    

    # incidence
    haz <- as.numeric(colMeans(tail(sim.comp$epi$ir100, 52)))
    out.haz <- round(data.frame(median = median(haz, na.rm = TRUE),
                                 ql = quantile(haz, qnt.low, names = FALSE, na.rm = TRUE),
                                 qu = quantile(haz, qnt.high, names = FALSE, na.rm = TRUE)), 2)
    haz.gc <- as.numeric(colMeans(tail(sim.comp$epi$ir100.gc, 52)))
    out.haz.gc <- round(data.frame(median = median(haz.gc, na.rm = TRUE),
                                    ql = quantile(haz.gc, qnt.low, names = FALSE, na.rm = TRUE),
                                    qu = quantile(haz.gc, qnt.high, names = FALSE, na.rm = TRUE)), 2)
    haz.ct <- as.numeric(colMeans(tail(sim.comp$epi$ir100.ct, 52)))
    out.haz.ct <- round(data.frame(median = median(haz.ct, na.rm = TRUE),
                                    ql = quantile(haz.ct, qnt.low, names = FALSE, na.rm = TRUE),
                                    qu = quantile(haz.ct, qnt.high, names = FALSE, na.rm = TRUE)), 2)
    haz.syph <- as.numeric(colMeans(tail(sim.comp$epi$ir100.syph, 52)))
    out.haz.syph <- round(data.frame(median = median(haz.syph, na.rm = TRUE),
                                   ql = quantile(haz.syph, qnt.low, names = FALSE, na.rm = TRUE),
                                   qu = quantile(haz.syph, qnt.high, names = FALSE, na.rm = TRUE)), 2)


    # HR
    num <- unname(colMeans(tail(sim.comp$epi$ir100, 52)))
    denom <- unname(colMeans(tail(sim.base$epi$ir100, 52)))
    vec.hr <- num/denom
    out.hr <- round(data.frame(median = median(vec.hr, na.rm = TRUE),
                               ql = quantile(vec.hr, qnt.low, names = FALSE, na.rm = TRUE),
                               qu = quantile(vec.hr, qnt.high, names = FALSE, na.rm = TRUE)), 2)

    num.gc <- unname(colMeans(tail(sim.comp$epi$ir100.gc, 52)))
    denom.gc <- unname(colMeans(tail(sim.base$epi$ir100.gc, 52)))
    vec.hr.gc <- num.gc/denom.gc
    vec.hr.gc <- vec.hr.gc[vec.hr.gc < Inf]
    out.hr.gc <- round(data.frame(median = median(vec.hr.gc, na.rm = TRUE),
                                  ql = quantile(vec.hr.gc, qnt.low, names = FALSE, na.rm = TRUE),
                                  qu = quantile(vec.hr.gc, qnt.high, names = FALSE, na.rm = TRUE)), 2)

    num.ct <- unname(colMeans(tail(sim.comp$epi$ir100.ct, 52)))
    denom.ct <- unname(colMeans(tail(sim.base$epi$ir100.ct, 52)))
    vec.hr.ct <- num.ct/denom.ct
    vec.hr.ct <- vec.hr.ct[vec.hr.ct < Inf]
    out.hr.ct <- round(data.frame(median = median(vec.hr.ct, na.rm = TRUE),
                                  ql = quantile(vec.hr.ct, qnt.low, names = FALSE, na.rm = TRUE),
                                  qu = quantile(vec.hr.ct, qnt.high, names = FALSE, na.rm = TRUE)), 2)
    
    num.syph <- unname(colMeans(tail(sim.comp$epi$ir100.syph, 52)))
    denom.syph <- unname(colMeans(tail(sim.base$epi$ir100.syph, 52)))
    vec.hr.syph <- num.syph/denom.syph
    vec.hr.syph <- vec.hr.syph[vec.hr.syph < Inf]
    out.hr.syph <- round(data.frame(median = median(vec.hr.syph, na.rm = TRUE),
                                  ql = quantile(vec.hr.syph, qnt.low, names = FALSE, na.rm = TRUE),
                                  qu = quantile(vec.hr.syph, qnt.high, names = FALSE, na.rm = TRUE)), 2)

    # NIA
    ir.comp <- unname(colMeans(sim.comp$epi$ir100)) * 1000
    vec.nia <- round(ir.base - ir.comp, 1)
    out.nia <- round(data.frame(median = median(vec.nia),
                                ql = quantile(vec.nia, qnt.low, names = FALSE),
                                qu = quantile(vec.nia, qnt.high, names = FALSE)), 0)

    ir.comp.gc <- unname(colMeans(sim.comp$epi$ir100.gc)) * 1000
    vec.nia.gc <- round(ir.base.gc - ir.comp.gc, 1)

    ir.comp.ct <- unname(colMeans(sim.comp$epi$ir100.ct)) * 1000
    vec.nia.ct <- round(ir.base.ct - ir.comp.ct, 1)
    
    ir.comp.syph <- unname(colMeans(sim.comp$epi$ir100.syph)) * 1000
    vec.nia.syph <- round(ir.base.syph - ir.comp.syph, 1)

    # PIA
    vec.pia <- vec.nia/ir.base
    out.pia <- round(data.frame(median = median(vec.pia),
                                ql = quantile(vec.pia, qnt.low, names = FALSE),
                                qu = quantile(vec.pia, qnt.high, names = FALSE)), 3)

    vec.pia.gc <- vec.nia.gc/ir.base.gc
    vec.pia.gc <- vec.pia.gc[vec.pia.gc > -Inf]
    out.pia.gc <- round(data.frame(median = median(vec.pia.gc),
                                   ql = quantile(vec.pia.gc, qnt.low, names = FALSE),
                                   qu = quantile(vec.pia.gc, qnt.high, names = FALSE)), 3)

    vec.pia.ct <- vec.nia.ct/ir.base.ct
    vec.pia.ct <- vec.pia.ct[vec.pia.ct > -Inf]
    out.pia.ct <- round(data.frame(median = median(vec.pia.ct),
                                   ql = quantile(vec.pia.ct, qnt.low, names = FALSE),
                                   qu = quantile(vec.pia.ct, qnt.high, names = FALSE)), 3)
    
    vec.pia.syph <- vec.nia.syph/ir.base.syph
    vec.pia.syph <- vec.pia.syph[vec.pia.syph > -Inf]
    out.pia.syph <- round(data.frame(median = median(vec.pia.syph),
                                   ql = quantile(vec.pia.syph, qnt.low, names = FALSE),
                                   qu = quantile(vec.pia.syph, qnt.high, names = FALSE)), 3)

    # browser()
    # NNT
    total.hiv.tests <- unname(colMeans(tail(sim.comp$epi$totalhivtests, 1)))
    gc.asympt.tests <- unname(colMeans(tail(sim.comp$epi$totalGCasympttests, 1)))
    ct.asympt.tests <- unname(colMeans(tail(sim.comp$epi$totalCTasympttests, 1)))
    syph.asympt.tests <- unname(colMeans(tail(sim.comp$epi$totalsyphasympttests, 1)))
    total.asympt.tests <- unname(colMeans(tail(sim.comp$epi$totalstiasympttests, 1)))
    
    vec.nnt <- total.asympt.tests/(median(incid.base) - unname(colSums(sim.comp$epi$incid)))
    out.nnt <- round(data.frame(median = median(vec.nnt),
                                ql = quantile(vec.nnt, qnt.low, names = FALSE),
                                qu = quantile(vec.nnt, qnt.high, names = FALSE)), 1)


    vec.nnt.gc <- gc.asympt.tests / (median(incid.base.gc) - unname(colSums(sim.comp$epi$incid.gc)))
    out.nnt.gc <- round(data.frame(median = median(vec.nnt.gc),
                                   ql = quantile(vec.nnt.gc, qnt.low, names = FALSE),
                                   qu = quantile(vec.nnt.gc, qnt.high, names = FALSE)), 1)

    vec.nnt.ct <- ct.asympt.tests / (median(incid.base.ct) - unname(colSums(sim.comp$epi$incid.ct)))
    out.nnt.ct <- round(data.frame(median = median(vec.nnt.ct),
                                   ql = quantile(vec.nnt.ct, qnt.low, names = FALSE),
                                   qu = quantile(vec.nnt.ct, qnt.high, names = FALSE)), 1)
    
    vec.nnt.syph <- syph.asympt.tests / (median(incid.base.syph) - unname(colSums(sim.comp$epi$incid.syph)))
    out.nnt.syph <- round(data.frame(median = median(vec.nnt.syph),
                                   ql = quantile(vec.nnt.syph, qnt.low, names = FALSE),
                                   qu = quantile(vec.nnt.syph, qnt.high, names = FALSE)), 1)

    cat("\n\nHIV Prevalence:")
    print(t(out.prev))
    
    cat("\n\nP&S Syph Prevalence:")
    print(t(out.primsecosyph.prev))

    cat("\nHIV Incidence:")
    print(t(out.haz))

    cat("\nHIV HR:")
    print(t(out.hr))

    # cat("\nHIV NIA:")
    # print(t(out.nia))

    cat("\nHIV PIA:")
    print(t(out.pia))

    cat("\nHIV NNT:")
    print(t(out.nnt))

    cat("\nGC Incidence:")
    print(t(out.haz.gc))

    cat("\nGC HR:")
    print(t(out.hr.gc))

    cat("\nGC PIA:")
    print(t(out.pia.gc))

    cat("\nGC NNT:")
    print(t(out.nnt.gc))

    cat("\nCT Incidence:")
    print(t(out.haz.ct))

    cat("\nCT HR:")
    print(t(out.hr.ct))

    cat("\nCT PIA:")
    print(t(out.pia.ct))

    cat("\nCT NNT:")
    print(t(out.nnt.ct))
    
    cat("\nSyph Incidence:")
    print(t(out.haz.syph))
    
    cat("\nSyph HR:")
    print(t(out.hr.syph))
    
    cat("\nSyph PIA:")
    print(t(out.pia.syph))
    
    cat("\nSyph NNT:")
    print(t(out.nnt.syph))

  } else {

    cat("\n\nHIV Prevalence:")
    print(t(prev.base))

    cat("\n\nP&S Syph Prevalence:")
    print(t(primsecosyph.prev.base))
    
    cat("\nHIV Incidence:")
    print(t(haz.base))

    cat("\nGC Incidence:")
    print(t(haz.base.gc))

    cat("\nCT Incidence:")
    print(t(haz.base.ct))
    
    cat("\nSyph Incidence:")
    print(t(haz.base.syph))

  }

}

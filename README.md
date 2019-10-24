## Proportion of Incident Human Immunodeficiency Virus Cases Among Men Who Have Sex With Men Attributable to Gonorrhea and Chlamydia: A Modeling Analysis

This repository holds the source to code to reproduce the analysis featured in our HIV/STI transmission model among men who have sex with men in the United States. This study investigated how prevalent gonorheaand and chlamydia infections could impact the acquisition and transmission risk of HIV infection.

### Citation

> Jones J, Weiss K, Mermin J, Dietz P, Rosenberg ES, Gift T, Chesson H, Sullivan PS, Lyles C, Bernstein K, Jenness SM. Proportion of Incident HIV Cases among Men Who Have Sex with Men Attributable to Gonorrhea and Chlamydia: A Modeling Analysis. Sexually Transmitted Diseases. 2019; 46(6): 357–363.

Additional details may be found in the journal article: https://www.ncbi.nlm.nih.gov/pubmed/31095100


### Abstract

#### BACKGROUND
Sexually transmitted infections (STIs) are associated with an increased risk of human immunodeficiency virus (HIV) acquisition and transmission. We estimated the proportion of HIV incidence among men who have sex with men attributable to infection with the 2 most common bacterial STIs, Neisseria gonorrhoeae (NG) and Chlamydia trachomatis (CT).

#### METHODS
We used a stochastic, agent-based model of a sexual network of MSM with cocirculating HIV, NG, and CT infections. Relative risk (RR) multipliers, specific to anatomic site of infection, modified the risk of HIV transmission and acquisition based on STI status. We estimated the effect of NG and CT on HIV incidence overall and on HIV acquisition and HIV transmission separately. Each scenario was simulated for 10 years. The population attributable fraction (PAF) was determined for each combination of RRs by comparing the incidence in the final year of a scenario to a scenario in which the RRs associated with NG and CT were set to 1.0.

#### RESULTS
Overall, 10.2% (interquartile range [IQR], 7.9-12.4) of HIV infections were attributable to NG/CT infection. Then in sensitivity analyses, the PAF for HIV transmission ranged from 3.1% (IQR, 0.5-5.2) to 20.4% (IQR, 17.8-22.5) and the PAF for HIV acquisition ranged from 2.0% (IQR, -0.7 to 4.3) to 13.8% (IQR, 11.7-16.0).

#### CONCLUSIONS
Despite challenges in estimating the causal impact of NG/CT on HIV risk, modeling is an alternative approach to quantifying plausible ranges of effects given uncertainty in the biological cofactors. Our estimates represent idealized public health interventions in which STI could be maximally prevented, setting targets for real-world STI interventions that seek to reduce HIV incidence.

## Model Code

These models are written and executed in the R statistical software language. To run these files, it is necessary to first install our epidemic modeling software, [EpiModel](http://epimodel.org/), and our extension package specifically for modeling HIV and STI transmission dynamics among MSM, [EpiModelHIV](http://github.com/statnet/EpiModelHIV).

In R:
```
install.packages("EpiModel", dep = TRUE)

# install remotes if necessary, install.packages("remotes")
remotes::install_github("statnet/tergmLite")
remotes::install_github("statnet/EpiModelHPC")
remotes::install_github("statnet/EpiModelHIV", ref = "paf")
```

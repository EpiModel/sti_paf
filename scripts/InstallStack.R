
# Update stiPrEP stack

system("git pull")
devtools::install_github("statnet/EpiModel")
devtools::install_github("statnet/EpiModelHPC")
devtools::install_github("statnet/tergmLite", subdir = "tergmLite")
devtools::install_github("statnet/EpiModelHIV", ref = "port-sti")


## interface with hyak

# upload scripts
system("scp scripts/burnin/*.burn.[Rs]* hyak:/gscratch/csde/sjenness/sti")
system("scp scripts/followup/*.fu.* hyak:/gscratch/csde/sjenness/sti")

# upload inputs
system("scp est/*.rda hyak:/gscratch/csde/sjenness/sti/est")

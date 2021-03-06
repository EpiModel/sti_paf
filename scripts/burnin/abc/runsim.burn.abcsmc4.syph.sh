#!/bin/bash

### User specs
#PBS -N sti-abc-syph
#PBS -l nodes=1:ppn=16,mem=50gb,feature=16core,walltime=05:00:00:00
#PBS -o /suppscr/csde/kweiss2/sti/out
#PBS -e /suppscr/csde/kweiss2/sti/out
#PBS -j oe
#PBS -d /suppscr/csde/kweiss2/sti
#PBS -m ae

### Standard specs
HYAK_NPE=$(wc -l < $PBS_NODEFILE)
HYAK_NNODES=$(uniq $PBS_NODEFILE | wc -l )
HYAK_TPN=$((HYAK_NPE/HYAK_NNODES))
NODEMEM=`grep MemTotal /proc/meminfo | awk '{print $2}'`
NODEFREE=$((NODEMEM-2097152))
MEMPERTASK=$((NODEFREE/HYAK_TPN))
ulimit -v $MEMPERTASK
export MX_RCACHE=0

### Modules
#module load r_3.2.4

### App
R CMD BATCH --vanilla sim.burn.abcsmc4.syph.R out/sim.burn.abcsmc4.syph.n${NSIM}.p${PACC}.Rout

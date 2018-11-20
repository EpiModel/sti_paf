#!/bin/bash

## Job Name
#SBATCH --job-name=slurm-test

## Nodes
#SBATCH --nodes=1

## Tasks per node
#SBATCH --ntasks-per-node=16

## Walltime
#SBATCH --time=2:00:00

## E-mail notification
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kweiss2@emory.edu

## Memory per node
#SBATCH --mem=50G

## Specify the working directory
#SBATCH --workdir=/suppscr/csde/kweiss2/sti

### Modules
. /suppscr/csde/sjenness/spack/share/spack/setup-env.sh
module load gcc-8.2.0-gcc-4.8.5-rhsxipz
module load r-3.5.1-gcc-8.2.0-4suigve

### App
Rscript sim.fu.paf.R

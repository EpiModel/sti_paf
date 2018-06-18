#!/bin/bash

qsub -q batch -t 1 runsim.reestim.casl.abc.sh
qsub -q batch -t 2-250 runsim.reestim.casl.abc.sh


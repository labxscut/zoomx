#!/usr/bin/env bash

#works from sake
#works from laptop 

nohup ./filterMolecule.R --debug >./filterMolecule.R.log 2>&1 &
nohup ./plotCircus.R --debug >./plotCircus.R.log 2>&1 &
nohup ./plotCumCov.R --debug >./plotCumCov.R.log 2>&1 &
nohup ./plotDenCov.R --debug >./plotDenCov.R.log 2>&1 &
nohup ./plotStats.R --debug >./plotStats.R.log 2>&1 &
nohup ./plotMolecule.R --debug >./plotMolecule.R.log 2>&1 &

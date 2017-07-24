#!/bin/bash
#to use this script need to set $cfgfile and $moleFile
moleFile=$1
bamFile=$2
addOpt=$3
log=$moleFile.log.plotMoleFile
cfgFile=$bamFile.gridScan.cfg.json        #output of gridScan configuraiton
maxSize=$(cat $moleFile | wc -l)
for i in `seq 1 $maxSize`; do
  plotMolCMD="plotMolecule.R $moleFile.Plot_$i $cfgFile $addOpt >>$log 2>&1"
  echo $plotMolCMD
  nohup $plotMolCMD &
done

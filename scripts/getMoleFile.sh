#!/bin/bash
if [[ -z $4 ]]; then 
  echo usage: getMolecule.sh \$moleFile \$bamFile \$cisSize \$sampleType
  exit
fi
moleFile=$1
bamFile=$2
cisSize=$3
sampleType=$4
cfgFile=$bamFile.gridScan.cfg.json        #output of gridScan configuraiton
log=$moleFile.getMoleFile.log
maxSize=$(cat $moleFile | wc -l)
for i in `seq 1 $maxSize`; do
  getMolCMD="getMolecule.sh $moleFile $bamFile $i $moleFile.Plot $cisSize $sampleType >>$log 2>&1"
  echo $getMolCMD
  nohup $getMolCMD &
done
